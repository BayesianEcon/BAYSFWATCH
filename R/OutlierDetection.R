library(MixMatrix)
library(ggplot2)

#' Outlier Detection in 3D Data
#'
#' This function detects outliers in a 3D data array using minimum Bayes Factor.
#' It can estimate necessary matrices using Maximum Likelihood if not provided.
#'
#' @param X3D A 3D array representing the dataset.
#' @param M (Optional) Mean matrix. If NULL, it is estimated using ML.
#' @param Sl (Optional) Row covariance matrix. If NULL, it is estimated using ML.
#' @param V (Optional) Column covariance matrix. If NULL, it is estimated using ML.
#' @param Training A vector indicating indices for training data. Default is `1:60`.
#' @param RollW Logical; whether to use rolling window approach. Default is TRUE.
#' @param vv Define factor vv for the reduction of the prior variance. Default is 0.1.
#' @param qq Type 1 Error. Default is 0.01.
#' @param plot Logical; whether to generate a plot of results. Default is TRUE.
#' @param verbose Logical; whether to print status messages. Default is FALSE.
#'
#' @return A data frame containing the outlier detection results.
#'
#' @details The function first estimates mean and covariance matrices if they are not provided, then tests for outliers
#' using the `OutlierBFJ` function via Minimum Bayes Factor.
#' The resulting dataset contains the following columns:
#'   - Column Outdet is -1 if the observed matrix is considered an outlier.
#'   - Column Hv reports the Minimum Bayes Factor.
#'   - Columns hlb and hub report the indeterminacy lower and upper bounds.
#' If `plot = TRUE`, it generates a plot of the outlier statistics.
#'
#' @examples
#' # Generate a Synthetic Dataset
#' Data<-GenData(5, 3, 6, Rcn = 3, Rrn = 2)
#' # Run Outlier Detection
#' results <- OutlierDetection(Data$X3D,  M = Data$M, Sl = Data$Sl, V =  Data$V)
#' @export
OutlierDetection <- function(X3D,
                             M = NULL, Sl = NULL, V = NULL,
                             Training = 1:60,
                             RollW = TRUE, vv = 0.1,
                             qq = 0.01, plot = TRUE, verbose = FALSE) {

  kk =  dim(X3D)[3]  - max(Training) #look for outlier from right-after the last training observation
  dtvs = 0.05 #Lower bound for inflation under-the-alternative parameter alpha_t. Default is 0.05.
  dte = 0.999 #Upper bound for inflation under-the-alternative parameter alpha_t. Default is 0.999.
  dts = 100   #Steps from lower bound @param dtvs to Upper bound @param dte. Default is 100.
  km = 100    #Number of iterations for kernel-based method. Default is 100.
  hl = 0.5    #Find dt such that P[ht<hl|0] = qq. Default is 0.5.
  step = 300  #Step size for analysis. Default is 300.

  # Estimate matrices if not provided
  estim <- FALSE
  if (is.null(M) || is.null(Sl) || is.null(V)) {
    if (verbose) message("Estimating matrices via ML...")
    Estimation <- MixMatrix::MLmatrixnorm(X3D[,,Training])
    M <- Estimation$mean
    Sl <- Estimation$U  # Assumption for normalization here
    V <- Estimation$V*Estimation$var
    estim <- TRUE
  } else {
    if (verbose) message("Using predefined matrices for M, Sl, and V.")
    estim <- FALSE
  }

  # Compute outliers
  if (verbose) message("Computing outliers...")
  OutObj <- OutlierBFJ(X3D, M, Sl, V, RollW, kk, vv, dtvs, dte, dts, km, hl, qq, step,estim)

  # Convert to data frame and add index
  OutObj <- data.frame(OutObj)
  OutObj$index <- seq_len(nrow(OutObj))

  # Generate plot if requested
  if (plot) {
    if (verbose) message("Generating plot...")
    library(ggplot2)
    p2 <- ggplot2::ggplot(OutObj, aes(x = index, y = Hv))
    p2 <- p2 +  labs(x = "", y = "")
    p2 <- p2 +  geom_line(linewidth = 0.3)
    p2 <- p2 +  geom_line(aes(x = index, y = hlb), linewidth = 0.3, linetype = "dashed", col = "red")
    p2 <- p2 + geom_line(aes(x = index, y = hub), linewidth = 0.3, linetype = "dashed", col = "red")
    p2 <- p2 +   labs(title = "Minimum BF")
    p2 <- p2 +  theme_bw()
    p2 <- p2 +  theme(legend.position = "none",
                      axis.text.x = element_text(size = 14, vjust = -1.1),
                      axis.text.y = element_text(size = 14),
                      axis.ticks.length = unit(-1.2, "mm"),  # Internal tick marks
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())

    print(p2)
  }


  if (verbose) message("Outlier detection completed.")
  return(OutObj)
}
