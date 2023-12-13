# MTWAS: Multi-tissue Transcriptome-Wide Association Studies

MTWAS is an R package for the paper

> Song, S., Wang, L., Hou, L., & Liu, J. S. (2023). Partitioning and aggregating cross-tissue and tissue-specific
genetic effects in identifying gene-trait associations.

## :white_check_mark: Prerequisites

The software is developed and tested in Linux and Windows environments.

- R (>=3.6)
- GNU Scientific Library (GSL) (>=2.3)

## :hammer_and_wrench: Installation

```r
devtools::install_github("szcf-weiya/MTWAS")
```

## :rocket: Getting Started

 :smiley: You could either download the pre-trained GTEx v8 MTWAS weights, or generate your own files!

## :one: Run MTWAS with pre-trained GTEx v8 weights (Recommended)

### Download the pre-trained files:

```bash
wget -O gtex_v8_mtwas_eqtls.tar.gz https://cloud.tsinghua.edu.cn/f/5633911d7c39431b8be8/?dl=1 --no-check-certificate
tar -zxvf gtex_v8_mtwas_eqtls.tar.gz
```

### MTWAS analysis:

We use chromosome 22 on whole blood as an example. The list of cell types is detailed in `ct_use.RData`.

```r
library(MTWAS)
data("summary_stats") ## EXAMPLE GWAS summary stats (could be specified by users, format: a data.frame with colnames: rsid, a1, a2, chr, z)
chr = 22
cell_type = 'Whole_Blood'
### remember to change the path to the downloaded folder!!!
load(paste0('./gtex_v8_mtwas_eqtls/twas_bim_chr', chr, '.RData'))  ## load twas bim files (downloaded)
load(paste0('./gtex_v8_mtwas_eqtls/', cell_type, '/twas_eqtl_chr', chr, '.RData')) ## load twas eqtl files (downloaded)
results = run_mtwas_easy(summary_stats, twas_bim, twas_eqtl) ## MTWAS summary statistics
head(results)
```



## :two: Train your own weights

```r
# load TWAS data
library(MTWAS)
data("twas_data")
# select cross-tissue eQTLs
ct.eQTL = select.ct.eQTL(twas_data, verbose = F, ncores = 1)
# select tissue-specific eQTLs
list.ts.eQTL = select.ts.eQTL(twas_data, ct.eQTL = ct.eQTL, ncores = 1)
```

### Extract cross-tissue eQTLs

```r
gene_name = 'CCT8L2' ## gene name
gene_index = which(names(ct.eQTL)==gene_name)
print(ct.eQTL[[gene_index]])
```

### Extract tissue-specific eQTLs

```r
gene_name = 'CCT8L2' ## gene name
gene_index = which(names(ct.eQTL)==gene_name)
tissue_index = 1
print(list.ts.eQTL[[tissue_index]][[gene_index]]$single.snp)
```

### Gene-trait association tests

```r
# load GWAS summary statistics (data.frame, colnames: rsid, a1, a2, chr, z)
data("summary_stats")
# association test
twas.single.trait(summary_stats, twas_data, list.ts.eQTL)
```

Data formats and function details: https://hohoweiya.xyz/MTWAS/articles/mtwas.html
