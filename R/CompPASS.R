#' CompPASS
#' Comparative Proteomic Analysis Software Suite (CompPASS) is based on spoke
#' model. This algorithm was developed by Dr. Mathew Sowa for defining the
#' human deubiquitinating enzyme interaction landscape
#' (Sowa, Mathew E., et al., 2009). The implementation of this algorithm was
#' inspired by Dr. Sowa's online tutorial
#' (\url{http://besra.hms.harvard.edu/ipmsmsdbs/cgi-bin/tutorial.cgi}).
#' The output includes Z-score, S-score, D-score and WD-score.
#'
#' @title CompPASS
#' @param datInput A dataframe with column names: idRun, idBait, idPrey,
#' countPrey. Each row represent one unique protein captured in one pull-down
#' experiment
#' @return A data frame consists of unique bait-prey pairs with Z-score,
#' S-score, D-score and WD-score indicating interacting probabilities.
#'
#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
#' @references Sowa, Mathew E., et al. "Defining the human deubiquitinating
#' enzyme interaction landscape." Cell 138.2 (2009): 389-403.
#' \url{https://doi.org/10.1016/j.cell.2009.04.042}
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @importFrom matrixStats rowSds
#' @examples
#' data(TestDatInput)
#' CompPASS(TestDatInput)
CompPASS <- function(datInput){
    colInput <-
        c("idRun", "idBait", "idPrey", "countPrey")
    
    if(!is.data.frame(datInput)){
        stop("Input data should be data.frame")
    }
    
    if(!all(colInput %in% colnames(datInput))){
        missingCol <-
            setdiff(colInput, 
                    colnames(datInput)[match(colInput, colnames(datInput))])
        stop(paste0("Input data missing: ", paste(missingCol, collapse = ", ")))
    }
    
    idBait <- NULL
    idPrey <- NULL
    countPrey <- NULL
    AvePSM <- NULL
    Stsc <- NULL
    Mtsc <- NULL
    f_sum <- NULL
    . <- NULL
    k <-
        length(unique(datInput$idBait))
    statsTbl <-
        datInput %>%
        group_by(`idBait`, `idPrey`) %>%
        summarise(`AvePSM` = mean(`countPrey`)) %>%
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
    stats <-
        spread(statsTbl[, c("idBait", "idPrey", "AvePSM")], `idBait`, `AvePSM`)
    statsTable <-
        as.matrix(stats[, -1])
    rownames(statsTable) <-
        stats$idPrey
    statsTable[is.na(statsTable)] <- 0
    z <-
        t(scale(t(statsTable), center = TRUE, scale = TRUE))
    zTable <-
        data.frame(`idBait` = unlist(lapply(colnames(z), function(x){
            rep(x, nrow(z))
        })),
        `idPrey` = rep(rownames(z), ncol(z)),
        `scoreZ` = c(z),
        stringsAsFactors = FALSE) %>%
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
    w <-
        data.frame(`idPrey` = rownames(statsTable),
                    `Mtsc` = rowMeans(statsTable),
                    `Stsc` = rowSds(statsTable),
                    stringsAsFactors = FALSE) %>%
        mutate(`w` = ifelse(`Stsc`/`Mtsc` <=1, 1, `Stsc`/`Mtsc`))
    f <-
        unique(datInput[, c("idBait", "idPrey")]) %>%
        group_by(`idPrey`) %>%
        summarise(`f_sum` = n())
    p <-
        datInput[, c("idBait", "idPrey")] %>%
        group_by(`idBait`, `idPrey`) %>%
        summarise(`p` = n()) %>%
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
    scoreTbl <-
        statsTbl %>%
        left_join(., f, by = "idPrey") %>%
        left_join(., p[, c("BP", "p")], by = "BP") %>%
        left_join(., w, by = "idPrey") %>%
        mutate(`k` = k) %>%
        mutate(`scoreS` = sqrt((`AvePSM`)*(`k`)/(`f_sum`))) %>%
        mutate(`scoreD` = sqrt((`AvePSM`)*(((`k`)/(`f_sum`))^`p`))) %>%
        mutate(`WD_inner` = (`k` / `f_sum`) * (`Stsc` / `Mtsc`)) %>%
        mutate(`scoreWD` = sqrt((`AvePSM`)*(((`k`)/(`f_sum`)*`w`)^`p`))) %>%
        left_join(., zTable[, c("BP", "scoreZ")], by = "BP")
    output <-
        as.data.frame(scoreTbl[, c("idBait", "idPrey",
                                    "AvePSM", "scoreZ",
                                    "scoreS", "scoreD", "scoreWD")])
    return(output)
}
