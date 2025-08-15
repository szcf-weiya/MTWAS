# MTWAS: Multi-tissue Transcriptome-Wide Association Studies

MTWAS is an R package for the paper

> Song, S., Wang, L., Hou, L., & Liu, J. S. (2024). Partitioning and aggregating cross-tissue and tissue-specific genetic effects to identify gene-trait associations. Nature Communications, 15(1), 5769. https://doi.org/10.1038/s41467-024-49924-4


:bulb: We provide pre-trained eQTL weights with MTWAS on 47 **GTEx version 8** tissues, 13 types of immune cells and 2 activation conditions of the **DICE** dataset, and 14 immune cell types from the single-cell RNA-seq **OneK1K** dataset.

:smiley: You could either run MTWAS with the pre-trained weights (:star:easy and recommended), or train your own models with the expression data!

## Table of contents
* [Prerequisites](#white_check_mark-prerequisites)
* [Installation](#hammer_and_wrench-installation)
* [Prepare GWAS summary statistics](#scroll-prepare-gwas-summary-statistics)
* [Run MTWAS with pre-trained GTEx v8 weights](#rocket-run-mtwas-with-pre-trained-gtex-v8-weights)
* [Run MTWAS with pre-trained DICE weights](#rocket-run-mtwas-with-pre-trained-dice-weights)
* [Run MTWAS with pre-trained OneK1K weights](#rocket-run-mtwas-with-pre-trained-onek1k-weights)
* [Train MTWAS with your own datasets](#key-train-your-own-weights)


## :white_check_mark: Prerequisites

The software is developed and tested in Linux and Windows environments.

- R (>=3.6)
- GNU Scientific Library (GSL) (>=2.3)

## :hammer_and_wrench: Installation

```r
devtools::install_github("szcf-weiya/MTWAS")
```

## :scroll: Prepare GWAS summary statistics
Please prepare the GWAS summary statistics in the following format (including the header line):
```
   chr        rsid     ref   alt       z         
    1      rs4040617    G     A     -0.199    
    1      rs4075116    C     T      0.646     
    1      rs9442385    T     G     -0.016    
    ...
```
**chr**: chromosome

**rsid**: SNP rsid

**ref**: reference allele

**alt**: alternative allele

**z**: GWAS z score


## :rocket: Run MTWAS with pre-trained GTEx v8 weights


### Download the pre-trained files:

```bash
wget -O gtex_v8_mtwas_eqtls.tar.gz "https://cloud.tsinghua.edu.cn/d/8bf84211cf3247a0b480/files/?p=%2Fgtex_v8_mtwas_eqtls.tar.gz&dl=1" --no-check-certificate
tar -zxvf gtex_v8_mtwas_eqtls.tar.gz
```

If you are working on UKBB phenotypes, please download from the following weights instead (We restricted eQTLs within UKBB SNPs):

```bash
wget -O gtex_v8_mtwas_eqtls.tar.gz "https://cloud.tsinghua.edu.cn/d/8bf84211cf3247a0b480/files/?p=%2Fgtex_v8_mtwas_eqtls_ukbb.tar.gz&dl=1" --no-check-certificate
tar -zxvf gtex_v8_mtwas_eqtls.tar.gz
```
The result table for UKBB phenotypes (including 84 self-reported cancer and non-cancer illness phenotypes with effective sample sizes larger than 5,000) can also be directly downloaded from the supplementary files of the paper.

### MTWAS analysis:

We use chromosome 22 on whole blood as an example. The list of tissues is detailed in `ct_use.RData`.

```r
library(MTWAS)
data("summary_stats") ## EXAMPLE GWAS summary stats (could be specified by users, format: a data.frame with colnames: chr, rsid, a1, a2, z)
chr = 22
cell_type = 'Whole_Blood'
### remember to change the path to the downloaded folder!!!
## load twas bim files (downloaded)
load(paste0('./gtex_v8_mtwas_eqtls/twas_bim_chr', chr, '.RData'))  
## load twas eqtl files (downloaded)
load(paste0('./gtex_v8_mtwas_eqtls/', cell_type, '/twas_eqtl_chr', chr, '.RData'))
## Run mtwas and derive the gene-trait association test statistics
results = run_mtwas_easy(summary_stats, twas_bim, twas_eqtl, pred_res) 
head(results)
```

The output `results` is a data.frame with the following format:

```
   gene        MTWAS_Z      MTWAS_P     pred_r2      pred_pv         
  PLA2G3        -1.87        0.062        0.01        0.02    
  PANX2         -1.83        0.066        0.01        0.08     
    ...
```

**gene**: gene name

**MTWAS_Z**: gene-trait association Z score derived by MTWAS

**MTWAS_P**: gene-trait association P value derived by MTWAS

**pred_r2**: prediction accuracy of the gene expression evaluated by $r^2$

**pred_pv**: prediction accuracy of the gene expression evaluated by an F-test

Note that we output results of all genes. Users could specify the criteria of the outputs, e.g.,`results[results$pred_pv < 0.05 & results$MTWAS_P < 5e-6, ]`.

## :rocket: Run MTWAS with pre-trained DICE weights



### Download the pre-trained files:

```bash
wget -O dice_mtwas_eqtls.tar.gz "https://cloud.tsinghua.edu.cn/d/8bf84211cf3247a0b480/files/?p=%2Fdice_mtwas_eqtls.tar.gz&dl=1" --no-check-certificate
tar -zxvf dice_mtwas_eqtls.tar.gz
```


### MTWAS analysis:

We use chromosome 22 on B naive cell line as an example. The list of cell types is detailed in `ct_use.RData` (remember to change the path to the downloaded folder).

```r
library(MTWAS)
data("summary_stats") ## EXAMPLE GWAS summary stats (could be specified by users, format: a data.frame with colnames: chr, rsid, a1, a2, z)
chr = 22
cell_type = 'B_NAIVE'
### remember to change the path to the downloaded folder!!!
## load twas bim files (downloaded)
load(paste0('./dice_mtwas_eqtls/twas_bim_chr', chr, '.RData'))  
## load twas eqtl files (downloaded)
load(paste0('./dice_mtwas_eqtls/', cell_type, '/twas_eqtl_chr', chr, '.RData'))
## Run mtwas and derive the gene-trait association test statistics
results = run_mtwas_easy(summary_stats, twas_bim, twas_eqtl, pred_res) 
head(results)
```

The output `results` is a data.frame with the following format:

```
   gene                 MTWAS_Z      MTWAS_P     pred_r2      pred_pv         
  ENSG00000232754         3.08        0.002        0.04        0.99    
  ENSG00000100225         2.46        0.014        0.07        0.01     
    ...
```

**gene**: gene name

**MTWAS_Z**: gene-trait association Z score derived by MTWAS

**MTWAS_P**: gene-trait association P value derived by MTWAS

**pred_r2**: prediction accuracy of the gene expression evaluated by $r^2$

**pred_pv**: prediction accuracy of the gene expression evaluated by an F-test

Note that we output results of all genes. Users could specify the criteria of the outputs, e.g.,`results[results$pred_pv < 0.05 & results$MTWAS_P < 5e-6, ]`.


## :rocket: Run MTWAS with pre-trained OneK1K weights


### Download the pre-trained files:

```bash
wget -O onek1k_mtwas_eqtls.tar.gz "https://cloud.tsinghua.edu.cn/d/8bf84211cf3247a0b480/files/?p=%2Fonek1k_mtwas_eqtls.tar.gz&dl=1" --no-check-certificate
tar -zxvf onek1k_mtwas_eqtls.tar.gz
```


### MTWAS analysis:

We use chromosome 22 on CD4 NC cell line as an example. The list of cell types is detailed in `ct_use.RData` (remember to change the path to the downloaded folder).

```r
library(MTWAS)
data("summary_stats") ## EXAMPLE GWAS summary stats (could be specified by users, format: a data.frame with colnames: chr, rsid, a1, a2, z)
chr = 22
cell_type = 'CD4_NC'
### remember to change the path to the downloaded folder!!!
## load twas bim files (downloaded)
load(paste0('./onek1k_mtwas_eqtls/twas_bim_chr', chr, '.RData'))  
## load twas eqtl files (downloaded)
load(paste0('./onek1k_mtwas_eqtls/', cell_type, '/twas_eqtl_chr', chr, '.RData'))
## Run mtwas and derive the gene-trait association test statistics
results = run_mtwas_easy(summary_stats, twas_bim, twas_eqtl, pred_res) 
head(results)
```

The output `results` is a data.frame with the following format:

```
   gene           MTWAS_Z      MTWAS_P     pred_r2      pred_pv         
  APOL6            1.80        0.072        0.001        0.65    
  RP6-109B7.3     -1.77        0.076        0.038       7.4e-09     
    ...
```

**gene**: gene name

**MTWAS_Z**: gene-trait association Z score derived by MTWAS

**MTWAS_P**: gene-trait association P value derived by MTWAS

**pred_r2**: prediction accuracy of the gene expression evaluated by $r^2$

**pred_pv**: prediction accuracy of the gene expression evaluated by an F-test

Note that we output results of all genes. Users could specify the criteria of the outputs, e.g.,`results[results$pred_pv < 0.05 & results$MTWAS_P < 5e-6, ]`.



## :key: Train your own weights

We also provide functions to train MTWAS with your own datasets!

### Step 1: Data preparation

In order to derive your own MTWAS weights, three types of data are necessary. There is a demo dataset built in our R package:

```r
library(MTWAS)
data('demo')
# demo$dat
# demo$E.info
# demo$E
```

#### (1) Genotype files (dat)

Format: a **list** of plink bfiles, including bim, fam, bed

One could use the R function `read_plink` in package `EBPRS` to read the plink files:

```r
library(EBPRS)
dat <- read_plink('PATH_TO_PLINK_BFILE')
```

#### (2) Gene information (E.info)

Genes that we are interested in to derive the gene-trait associations. 

Format: a **data.frame** in the following format (including the header line):

```
   chr        start        end        gene         
    22      15528192    15529139     OR11H1   
    22      15690026    15721631     POTEH          
    ...
```

#### (3) Gene expression data (E)

A **list** including all the **gene expression data**, the length of the list is the number of tissues. Each element is a sample*gene matrix. 

Each matrix should have **rownames** (overlapped with `dat$fam$V2`), and **colnames** (overlapped with `E.info$gene`)

The orders of the columns and rows of the matrices are not necessary to be the same. The matrices can have NAs.

:exclamation: Please NOTE that the position of `dat$bim$V4` and `E.info` should be the same build (e.g., both are hg19, or hg38, or etc.)


### Step 2: Data imputation and formatting
```r
library(MTWAS)
data('demo')
### substitute the input with your own dataset
twas_dat <- format_twas_dat(E=demo$E, E.info=demo$E.info, dat=demo$dat) 
names(twas_dat)
```

### Step 3: Model training

```r
# load TWAS data
# select cross-tissue eQTLs
ct.eQTL = select.ct.eQTL(twas_dat, verbose = F, ncores = 1)
# select tissue-specific eQTLs
list.eQTL = select.ts.eQTL(twas_dat, ct.eQTL = ct.eQTL, ncores = 1)
```

### Step 4: Extract cross-tissue eQTLs

```r
gene_name = 'IL17RA' ## gene name
gene_index = which(names(ct.eQTL)==gene_name)
## ct-eQTLs for gene IL17RA
print(list.eQTL[[1]][[gene_index]]$`common.snp`) 
```

### Step 5: Extract tissue-specific eQTLs

```r
gene_name = 'CCT8L2' ## gene name
gene_index = which(names(ct.eQTL)==gene_name)
tissue_index = 1 ## tissue specific
## ts-eQTLs for gene CCT8L2 on tissue 1
print(list.eQTL[[tissue_index]][[gene_index]]$`single.snp`)
```

### Step 6: Gene-trait association tests

```r
# load GWAS summary statistics (data.frame, colnames: rsid, a1, a2, chr, z)
data("summary_stats")
# association test
twas.single.trait(summary_stats, twas_dat, list.eQTL)
```

The details of output and data formats can be found in the auto-generated vignette: https://hohoweiya.xyz/MTWAS/articles/mtwas.html

## References
**MTWAS software**
> Song, S., Wang, L., Hou, L., & Liu, J. S. (2024). Partitioning and aggregating cross-tissue and tissue-specific genetic effects to identify gene-trait associations. Nature Communications, 15(1), 5769. https://doi.org/10.1038/s41467-024-49924-4

**The GTEx dataset**
> https://gtexportal.org/home/

**The DICE dataset**
> https://dice-database.org/

> Schmiedel, Benjamin J., et al. "Impact of genetic polymorphisms on human immune cell gene expression." Cell 175.6 (2018): 1701-1715.

**The OneK1K dataset**
> https://onek1k.org/

> Yazar, Seyhan, et al. "Single-cell eQTL mapping identifies cell type–specific genetic control of autoimmune disease." Science 376.6589 (2022): eabf3041.
