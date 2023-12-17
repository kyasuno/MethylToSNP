#' This method identifies novel SNPs in Illumina DNA methylation array data
#'
#' @docType package
#' @name MethylToSNP
NULL


# require(Ckmeans.1d.dp)

#' Identify sites that may have underlying SNPs in methylation array data
#' 
#' @param data A matrix or a data frame or an GenomicRatioSet, GenomicMethylSet, MethylSet, or RatioSet object (see minfi package)
#' @param min.obs Minimum number of non-missing data required to test (default: 50). 
#' @param gap.ratio The ratio of two gaps should be above the threshold.
#' @param gap.sum.ratio The ratio of the sum of two gaps relative to the total range of values should be above the threshold.
#' @param outlier.sd Do not consider outliers that are more than the specified number of standard deviations from the cluster center.
#' @param verbose Show additional information. Useful for debugging.
#' @param show.plot plot distribution of betas of a particular probe with information of weights, clusters and outliers. 
#'   TRUE is allowed only for <= 10 probes in data matrix. 
#' @param SNP Optional SNP annotation
#' @return Detected probes with 3-tier SNP-like methylation pattern along with their reliability scores and SNP annotation
#' @examples
#' MethylToSNP(
#'	data.frame(row.names=c('cg0000001', 'cg0000002'), sample1=c(0.1, 0.5), sample2=c(0.5, 0.9), sample3=c(0.4, 0.8)))
#' @export
MethylToSNP <- function(data, min.obs = 50, gap.ratio = 0.75, gap.sum.ratio = 0.5, outlier.sd = 3.0, verbose=FALSE, show.plot=FALSE, SNP=NULL) {
  
  if ( (!length(dim(data)) == 2) || (dim(data)[2] < 2)){
    stop("[MethylToSNP] There have to be at least 2 samples")
  }
  
  if (dim(data)[2] < 50){
    message("[MethylToSNP] Warning, SNP detection may be unreliable in datasets with less than 50 samples")
  }
  
  if (gap.ratio >= 1 || gap.ratio <= 0){
    stop("[MethylToSNP] You have entered an unacceptable gap.ratio value. This must be a number between zero and one.")
  }
  
  if (gap.sum.ratio <= 0 || gap.sum.ratio >= 1){
    stop("[MethylToSNP] You have entered an unacceptable gap.sum.ratio value. This must be a number between zero and one.")
  }
  
  if (is(data, "GenomicRatioSet") || is(data, "GenomicMethylSet") || is(data, "MethylSet") || is(data, "RatioSet")) {
    if (verbose) 
      message("[MethylToSNP] Calculating beta matrix. minfi is required")
    .checkMinfi()
    data <- getBeta(data)
  }
  
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("[MethylToSNP] Input must be a numeric matrix of methylation beta values")
  }
  
  if (is.null(SNP)){
    message("[MethylToSNP] Optionally, specify SNPs in a data frame with row names corresponding to cg probes (such as SNPs.147CommonSingle in minfiData or minfiDataEPIC package)")
  }
  
  ##########
  # STEP 1.
  # Identify probes that could be potential SNPs in a given set of samples
  
  potentialSNPs <- NULL
  
  # filter probes by the span (width) of data
  n_obs <- apply(data, 1, function(x) sum(!is.na(x)))
  spans <- apply(matrixStats::rowRanges(data, na.rm=TRUE), 1, diff)
  use <- spans >= 0.5 & n_obs >= min.obs
  message(Sys.time(), ": ", sum(use), " (", round(100*sum(use)/length(use), 1), "%)",
          " probes with the span >= 0.5 and ", 
          "at least ", min.obs, " non-missing values will be analyzed.")
  
  if (sum(use) == 0) {
    message(Sys.time(), ": No probe is likely to be a SNP. Return NULL.")
    return(NULL)
  }
  
  data <- data[use, , drop=FALSE]
  spans <- spans[use]
  
  if (show.plot & nrow(data) > 10) {
    stop("As the data has > 10 rows to be tested, show.plot is disabled.")
    show.plot <- FALSE
  }
  
  probes <- rownames(data)
  
  # for each probe:
  for (i in 1: dim(data)[1]) {
    
    # set variables
    probe <- probes[i]
    span <- spans[i]
    x <- data[i, ] |> na.omit()
    
    if (length(x) < 3) {
      next
    }
    
    if (span >= 0.5) {
      
      # inverse density of beta values
      weights <- 1.0 / approxfun(density(x))(x)
      
      # optimal 1D clustering with dynamic programming into 3 clusters
      # no randomization involved
      kmeans <- Ckmeans.1d.dp(x=t(x), y=weights, k=3)
      clusters <- kmeans$cluster
      centers <- kmeans$centers
      
      top <- which(centers == max(centers))
      bottom <- which(centers == min(centers))
      middle <- c(1:3)[-c(top, bottom)]
      
      # disregard outliers (assign them to non-existing cluster #4)
      # corrected
      if (outlier.sd != FALSE) {
        top_outliers <- which(abs(scale(x[clusters == top])) > outlier.sd)
        if (length(top_outliers) > 0) {
          top_outliers <- which(clusters == top)[top_outliers]
          clusters[top_outliers] <- 4
        }
        
        middle_outliers <- which(abs(scale(x[clusters == middle])) > outlier.sd)
        if (length(middle_outliers) > 0) {
          middle_outliers <- which(clusters == middle)[middle_outliers]
          clusters[middle_outliers] <- 4
        }
        
        bottom_outliers <- which(abs(scale(x[clusters == bottom])) > outlier.sd)
        if (length(bottom_outliers) > 0) {
          bottom_outliers <- which(clusters == bottom)[bottom_outliers]
          clusters[bottom_outliers] <- 4
        }
        
        # if we removed outliers and one of the clusters is empty -- should not happen
        if (length(x[clusters == top]) * length(x[clusters == middle]) * length(x[clusters == bottom]) == 0) {
          next
        }
      }
      
      # two gaps:
      # min of the top cluster minus max of the middle cluster
      # and min of the middle cluster minus max of the bottom cluster			
      gaps.sorted <- sort(c(
        min(x[clusters == top]) - max(x[clusters == middle]),
        min(x[clusters == middle]) - max(x[clusters == bottom])),
        decreasing = TRUE)
      
      gaps.largest <- gaps.sorted[1]
      gaps.smallest <- gaps.sorted[2]

      if (show.plot) {
        d <- tibble(
          x=x, 
          w=weights, 
          Cluster=factor(clusters, levels=c(1,2,3,4))
        ) 
        p1 <- d |> 
          ggplot(aes(x=x, y=w, colour=Cluster)) + theme_pubr(legend="right") + 
          geom_point() + 
          scale_colour_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#999999"), drop=FALSE) +
          scale_x_continuous(breaks=seq(0,1,by=0.1), limits=c(0,1)) +
          labs(x="beta", y="Weights", title=probe, 
               caption=paste0("Gap sizes: ", round(gaps.smallest, digits=3), " and ", 
                              round(gaps.largest, digits=3), "\n", 
                              "Gap ratio = ", round(gaps.smallest / gaps.largest, digits=3), 
                              "; and\nProp gap vs. width = ", round(sum(gaps.sorted)/span, digits=3)))
        p2 <- d |> 
          ggplot(aes(x=x, fill=Cluster)) + theme_pubr(legend="right") + 
          geom_histogram(boundary=0, binwidth=0.025) + 
          scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "#999999"), drop=FALSE) +
          scale_x_continuous(breaks=seq(0,1,by=0.1), limits=c(0,1)) +
          labs(x="beta")
        print(
          plot_grid(p1, p2, ncol=1, align="hv", axis="lr")
        )
      }
      
      ###
      # Apply gap thresholds to decide whether to add a cg probe to the list of potential SNPs
      #
      if ((gaps.smallest >= gap.ratio * gaps.largest) && (sum(gaps.sorted) >= gap.sum.ratio * span)) {
        if (verbose) {
          message(probe)
        }
        potentialSNPs <- append(potentialSNPs, probe)
      }
    } else {
      # data span (max - min) is too narrow
    }
    
    ###
    # Progress indicator
    #
    if (verbose){
      if((i %% 1000) == 0){
        message('[MethylToSNP] Processed: ', i, " Identified potential SNPs: ", length(potentialSNPs))
      }
    }
  }
  ###
  # Summary
  #
  if(verbose){
    if (length(potentialSNPs) > 0 ){
      message("[MethylToSNP] Number of potential SNPs found: ",length(potentialSNPs))
    } else{
      warning("[MethylToSNP]  No potential SNPs found. NULL value is returned.")
      return(NULL)
    }
  }
  ##########
  # STEP 2.
  #
  # Calculate SNP confidence
  # 
  data <- data[potentialSNPs, , drop=FALSE]
  
  counts <- apply(data, 1, function(x) {
    c(low = sum(x <= 0.25, na.rm = TRUE),
      high = sum(x >= 0.75, na.rm = TRUE),
      mid = sum(x > 0.25 & x < 0.75, na.rm = TRUE))
  })
  
  L <- (counts["low", ] > 0) * (counts["mid", ] > 0) * (counts["high", ] > 0)
  snp.conf <- round(L * (counts["high", ] + 0.5 * counts["mid", ]) / colSums(counts), digits = 3)
  
  ######
  # Concatenate with SNP information
  #
  results <- tibble(
    row.names = potentialSNPs,
    confidence = snp.conf,
    samples_low = counts["low", ],
    samples_mid = counts["mid", ],
    samples_high = counts["high", ]
  )

  ######
  # Concatenate with SNP information
  #
  # results <- NULL
  if (!missing(SNP)) {
    SNP <- as_tibble(SNP[potentialSNPs, ])
    results <- bind_cols(results, SNP)
  }
  return(results)
}


.checkMinfi <- function() {
  if (!exists('getBeta')) {
    if (!is_installed('minfi')) {
      stop("[MethylToSNP] Bioconductor minfi package is required to ")
    } else {
      library('minfi')
    }        	
  }
}


#' Plot distribution of beta values for array probes that were identified as SNPs
#' 
#' @param x a data frame with probes identified as SNPs by MethylToSNP()
#' @param betas A matrix or a data frame or an GenomicRatioSet, GenomicMethylSet, MethylSet, or RatioSet object (see minfi package)
#' @param horizontal plot orientation (horizontal=TRUE by default)
#' @examples
#' x <- MethylToSNP(betas)
#' x <- plotPotentialSNPs(x, data)
#' @export
plotPotentialSNPs <- function(x, betas, horizontal=TRUE) {
  plotProbes(betas[rownames(x), ], horizontal=horizontal)
}


#' Plot distribution of beta values for array probes
#' 
#' @param betas A matrix or a data frame or an GenomicRatioSet, GenomicMethylSet, MethylSet, or RatioSet object (see minfi package)
#' @param horizontal plot orientation (horizontal=TRUE by default)
#' @examples
#' x <- MethylToSNP(betas)
#' x <- plotPotentialSNPs(x, data)
#' @export
plotProbes <- function(betas, horizontal=TRUE) {
  if (is(betas, "GenomicRatioSet") || is(betas, "GenomicMethylSet") || is(betas, "MethylSet") || is(betas, "RatioSet")) {
    message("[MethylToSNP] Extracting beta values. Minfi required")
    .checkMinfi()
    betas <- getBeta(betas)
  }
  
  old_par <- par(no.readonly = TRUE)
  if (horizontal) {
    par(mar=c(6,5,2,2)) # bottom, left, top and right 
    stripchart(as.data.frame(t(betas)), ylab='Beta value', pch='-', vertical=TRUE, las=2, cex.axis=0.6) # las=horizontal)
  } else {
    par(mar=c(2,8,2,2))
    stripchart(as.data.frame(betas), xlab='Beta value', pch='|', las=2, cex=0.6)
  }
  par(old_par)
}

plotProbesMerged <- function(betas){
  old_par <- par(no.readonly = TRUE)
  par(mar=c(4,8,2,2))
  d <- dim(betas)
  plot(betas[1], ylim=c(0,1), xlab="Samples", ylab="Beta value")
  for (i in 2:d[1]) {
    plot(betas[i], add=TRUE)
  }
  par(old_par)
}

