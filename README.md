# MTWAS: Multi-tissue Transcriptome-Wide Association Studies

MTWAS is an R package for the paper

> Song, S., Wang, L., Hou, L., & Liu, J. S. (2023). Partitioning and aggregating cross-tissue and tissue-specific
genetic effects in identifying gene-trait associations.


:bulb: We provide pre-trained eQTL weights with MTWAS on 47 **GTEx version 8** tissues, 13 types of immune cells and 2 activation conditions of the **DICE** dataset, and 14 immune cell types from the single-cell RNA-seq **OneK1K** dataset.

:smiley: You could either download the pre-trained MTWAS weights (easy and recommended), or generate your own files!

## Table of contents
* [Prerequisites](#white-check-mark-prerequisites)
* [Install](#hammer_and_wrench-install)
* [LD prepared](#rocket-getting-started)
* [Estimation of heritability and inflation factor](#rocket-estimation-of-heritability-and-inflation-factor)
* [Output](#bulb-output)
* [A Simplified Pipeline](#key-a-simplified-pipeline)


## :white_check_mark: Prerequisites

The software is developed and tested in Linux and Windows environments.

- R (>=3.6)
- GNU Scientific Library (GSL) (>=2.3)

## :hammer_and_wrench: Installation

```r
devtools::install_github("szcf-weiya/MTWAS")
```

## :rocket: Getting Started


## :one: Run MTWAS with pre-trained GTEx v8 weights (:star:Easy and Recommended)

### Download the pre-trained files:

```bash
wget -O gtex_v8_mtwas_eqtls.tar.gz https://cloud.tsinghua.edu.cn/f/5633911d7c39431b8be8/?dl=1 --no-check-certificate
tar -zxvf gtex_v8_mtwas_eqtls.tar.gz
```

### Prepare the GWAS summary statistics in the following format (including the header line):

```
   chr        rsid      a1    a2       z         
    1      rs4040617    G     A     -0.199    
    1      rs4075116    C     T      0.646     
    1      rs9442385    T     G     -0.016    
    ...
```
**chr**: chromosome

**rsid**: SNP rsid

**a1**: reference allele

**a2**: alternative allele

**z**: GWAS z score

### MTWAS analysis:

We use chromosome 22 on whole blood as an example. The list of cell types is detailed in `ct_use.RData`.

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
results = run_mtwas_easy(summary_stats, twas_bim, twas_eqtl) 
head(results)
```



## :two: Train your own weights

We also provide functions to train MTWAS with your own datasets!

### Data preparation

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


### Data imputation and formatting
```r
library(MTWAS)
data('demo')
### substitute the input with your own dataset
twas_dat <- format_twas_dat(E=demo$E, E.info=demo$E.info, dat=demo$dat) 
names(twas_dat)
```

### Model training

```r
# load TWAS data
# select cross-tissue eQTLs
ct.eQTL = select.ct.eQTL(twas_dat, verbose = F, ncores = 1)
# select tissue-specific eQTLs
list.eQTL = select.ts.eQTL(twas_dat, ct.eQTL = ct.eQTL, ncores = 1)
```

### Extract cross-tissue eQTLs

```r
gene_name = 'IL17RA' ## gene name
gene_index = which(names(ct.eQTL)==gene_name)
## ct-eQTLs for gene IL17RA
print(list.eQTL[[1]][[gene_index]]$`common.snp`) 
```

### Extract tissue-specific eQTLs

```r
gene_name = 'CCT8L2' ## gene name
gene_index = which(names(ct.eQTL)==gene_name)
tissue_index = 1 ## tissue specific
## ts-eQTLs for gene CCT8L2 on tissue 1
print(list.eQTL[[tissue_index]][[gene_index]]$`single.snp`)
```

### Gene-trait association tests

```r
# load GWAS summary statistics (data.frame, colnames: rsid, a1, a2, chr, z)
data("summary_stats")
# association test
twas.single.trait(summary_stats, twas_dat, list.eQTL)
```

The details of output and data formats can be found in the auto-generated vignette: https://hohoweiya.xyz/MTWAS/articles/mtwas.html
