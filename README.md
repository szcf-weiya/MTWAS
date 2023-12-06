# MTWAS: Multi-tissue Transcriptome-Wide Association Studies

MTWAS is an R package for the paper

> Song, S., Wang, L., Hou, L., & Liu, J. S. (2023). Partitioning and aggregating cross-tissue and tissue-specific
genetic effects in identifying gene-trait associations.

## :hammer_and_wrench: Installation

```r
devtools::install_github("szcf-weiya/MTWAS")
```

## :rocket: Get Started

```r
library(MTWAS)
```

### Train weights

```r
# load TWAS data
data("twas_data")
# select cross-tissue eQTLs
ct.eQTL = select.ct.eQTL(twas_data, verbose = F, ncores = 1)
# select tissue-specific eQTLs
list.ts.eQTL = select.ts.eQTL(twas_data, ct.eQTL = ct.eQTL, ncores = 1)
```

#### Extract cross-tissue eQTLs

```r
gene_name = 'CCT8L2' ## gene name
gene_index = which(names(ct.eQTL)==gene_name)
print(ct.eQTL[[gene_index]])
```

#### Extract tissue-specific eQTLs

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
