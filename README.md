# Statistical Modelling of AP-MS Data (SMAD)

[![Build Status](https://travis-ci.org/zqzneptune/SMAD.svg?branch=master)](https://travis-ci.org/zqzneptune/SMAD)

This R package implements statistical modelling of affinity purification–mass spectrometry (AP-MS) data to compute confidence scores to identify *bona fide* protein-protein interactions (PPI).

## Installation

The development version can be installed through github:
```{r}
 devtools::install_github(repo="zqzneptune/SMAD")
 library(SMAD)
```
## Quick start

### 1. CompPASS

Comparative Proteomic Analysis Software Suite (CompPASS) is based on spoke model. This algorithm was developed by Dr. Mathew Sowa for defining the human deubiquitinating enzyme interaction landscape [(Sowa, Mathew E., et al., 2009)](https://doi.org/10.1016/j.cell.2009.04.042). The implementation of this algorithm was inspired by Dr. Sowa's [online tutorial](http://besra.hms.harvard.edu/ipmsmsdbs/cgi-bin/tutorial.cgi). The output includes Z-score, S-score, D-score and WD-score.

Prepare input data into the dataframe *datInput* with the following format:

|idRun|idBait|idPrey|countPrey|
|-----|:----:|:----:|:-------:|
|Unique ID of one AP-MS run|Bait ID|Prey ID|Prey peptide count|

Then run:

```{r}
CompPASS(datInput)
```

### 2. CompPASS-Plus

On the basis of CompPASS, CompPASS-plus computes entropy and normalized 
WD-score. In its implementation in BioPlex 1.0 [(Huttlin, Edward L., et al., 2015)](https://doi.org/10.1016/j.cell.2015.06.043) and 
BioPlex 2.0 [(Huttlin, Edward L., et al., 2017)](https://www.nature.com/articles/nature22366), a naive 
Bayes classifier that learns to distinguish true interacting proteins from 
non-specific background and false positive identifications was included in the 
compPASS pipline. This function was optimized from the [source code](https://github.com/dnusinow/cRomppass).

Prepare input data into the dataframe *datInput* with the following format:

|idRun|idBait|idPrey|countPrey|
|-----|:----:|:----:|:-------:|
|Unique ID of one AP-MS run|Bait ID|Prey ID|Prey peptide count|

Then run:

```{r}
CompPASSplus(datInput)
```

### 3. HGScore

HGScore Scoring algorithm based on a hypergeometric distribution error model [(Hart et al., 2007)](https://doi.org/10.1186/1471-2105-8-236) with incorporation of NSAF [(Zybailov, Boris, et al., 2006)](https://doi.org/10.1021/pr060161n). This algorithm was first introduced to predict the protein complex network of Drosophila melanogaster [(Guruharsha, K. G., et al., 2011)](https://doi.org/10.1016/j.cell.2011.08.047). This scoring algorithm was based on matrix model. Unlike CompPASS, we need protein length for each prey in the additional column.

Prepare input data into the dataframe *datInput* with the following format:

|idRun|idBait|idPrey|countPrey|lenPrey|
|-----|:----:|:----:|:-------:|:-------:|
|Unique ID of one AP-MS run|Bait ID|Prey ID|Prey peptide count|Prey protein length|


Then run:

```{r}
HG(datInput)
```
## License

MIT @ Qingzhou Zhang
