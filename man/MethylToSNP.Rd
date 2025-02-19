% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MethylToSNP.R
\docType{package}
\name{MethylToSNP}
\alias{MethylToSNP}
\title{This method identifies novel SNPs in Illumina DNA methylation array data}
\usage{
MethylToSNP(
  data,
  min.obs = 50,
  gap.ratio = 0.75,
  gap.sum.ratio = 0.5,
  outlier.sd = 3,
  verbose = FALSE,
  show.plot = FALSE,
  SNP = NULL
)
}
\arguments{
\item{data}{A matrix or a data frame or an GenomicRatioSet, GenomicMethylSet, MethylSet, or RatioSet object (see minfi package)}

\item{min.obs}{Minimum number of non-missing data required to test (default: 50).}

\item{gap.ratio}{The ratio of two gaps should be above the threshold.}

\item{gap.sum.ratio}{The ratio of the sum of two gaps relative to the total range of values should be above the threshold.}

\item{outlier.sd}{Do not consider outliers that are more than the specified number of standard deviations from the cluster center.}

\item{verbose}{Show additional information. Useful for debugging.}

\item{show.plot}{plot distribution of betas of a particular probe with information of weights, clusters and outliers. 
TRUE is allowed only for <= 10 probes in data matrix.}

\item{SNP}{Optional SNP annotation}
}
\value{
Detected probes with 3-tier SNP-like methylation pattern along with their reliability scores and SNP annotation
}
\examples{
MethylToSNP(
data.frame(row.names=c('cg0000001', 'cg0000002'), sample1=c(0.1, 0.5), sample2=c(0.5, 0.9), sample3=c(0.4, 0.8)))
}
