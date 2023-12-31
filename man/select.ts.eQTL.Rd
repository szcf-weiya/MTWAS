% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/train_weights.R
\name{select.ts.eQTL}
\alias{select.ts.eQTL}
\title{Select tissue-specific eQTLs}
\usage{
select.ts.eQTL(x, ct.eQTL, id.tissue = -1, ncores = 1, ...)
}
\arguments{
\item{x}{data}

\item{ct.eQTL}{selected cross-tissue eQTLs from \code{select.ct.eQTL}}

\item{id.tissue}{the index of tissue (cell type) to be considered. Default is -1,
which returns the ts.eQTL for all tissues.}

\item{ncores}{number of cores. If larger than 1, then use parallel computing.}

\item{...}{
  Arguments passed on to \code{\link[=soda_model_imp_single]{soda_model_imp_single}}
  \describe{
    \item{\code{gam}}{tuning parameter \eqn{\gamma} in EBIC criteria. Default is 1.}
    \item{\code{scale}}{whether to scale the expression and the genetic data when regressing the expression on the genetic data. Default is \code{FALSE}.}
    \item{\code{nF0}}{minimum size of the set}
  }}
}
\description{
Select tissue-specific eQTLs
}
