#' Impute the list of expression matrix `E`
#' @param x twas data list
#' @return imputed twas data list and the imputation MSE error
#' @importFrom missForest missForest
#' @export
impute.E = function(x) {
    ntissue = length(x$E)
    mse = vector("list", length = ntissue)
    for (j in 1:ntissue) {
        nomiss = which(apply(x$E[[j]], 2, function(z) {any(!is.na(z))}))
        mat = x$E[[j]][, nomiss]
        colnames(mat) = NULL
        res = missForest(mat, variablewise = T)
        x$E[[j]][, nomiss] = res$ximp
        msej = numeric(length(nomiss))
        msej[nomiss] = res$OOBerror
        mse[[j]] = msej
    }
    list(imp = x, mse.error = mse)
}
