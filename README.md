# MTWAS

## :hammer_and_wrench: Installation

```r
devtools::install_github("szcf-weiya/MTWAS")
```

## :rocket: Get Started

```r
# load TWAS data
data("twas_data")
# select cross-tissue eQTLs
ct.eQTL = select.ct.eQTL(twas_data, verbose = F, ncores = 1)
# select tissue-specific eQTLs
list.ts.eQTL = select.ts.eQTL(twas_data, ct.eQTL = ct.eQTL, ncores = 1)
# load summary statistics
data("summary_stats")
# association test
twas.single.trait(summary_stats, twas_data, list.ts.eQTL)
```
