# code to prepare demo data (take a subset of 401 genes and 47 cell types)
ngene = 3
ntissue = 5
idx.max.snp = max(sapply(1:ngene, function(i) max(pos[[i]]))) # 509
twas_data = list(
    E = lapply(1:ntissue, function(i) E[[i]][, 1:ngene]),
    dat = list(
        bed = dat$bed[, 1:idx.max.snp],
        bim = dat$bim[1:idx.max.snp, ],
        fam = dat$fam
    ),
    E.info = E.info[1:ngene,], # data.table -> data.frame
    pos = pos[1:ngene]
)
#save(twas_data, file = "data/twas_data.RData")
usethis::use_data(twas_data, overwrite = T) # smaller size
stats <- read.table("../twas_demo/gwas_summary_stats.txt", header = T) # 17301 rows
# take first 1000 rows
summary_stats = stats[1:1000, ]
usethis::use_data(summary_stats, overwrite = T)
