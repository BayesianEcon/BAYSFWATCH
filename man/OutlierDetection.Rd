% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OutlierDetection.R
\name{OutlierDetection}
\alias{OutlierDetection}
\title{Outlier Detection in 3D Data}
\usage{
OutlierDetection(
  X3D,
  M = NULL,
  Sl = NULL,
  V = NULL,
  Training = 1:60,
  RollW = TRUE,
  vv = 0.1,
  qq = 0.01,
  plot = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{X3D}{A 3D array representing the dataset.}

\item{M}{(Optional) Mean matrix. If NULL, it is estimated using ML.}

\item{Sl}{(Optional) Row covariance matrix. If NULL, it is estimated using ML.}

\item{V}{(Optional) Column covariance matrix. If NULL, it is estimated using ML.}

\item{Training}{A vector indicating indices for training data. Default is \code{1:60}.}

\item{RollW}{Logical; whether to use rolling window approach. Default is TRUE.}

\item{vv}{Define factor vv for the reduction of the prior variance. Default is 0.1.}

\item{qq}{Type 1 Error. Default is 0.01.}

\item{plot}{Logical; whether to generate a plot of results. Default is TRUE.}

\item{verbose}{Logical; whether to print status messages. Default is FALSE.}
}
\value{
A data frame containing the outlier detection results.
}
\description{
This function detects outliers in a 3D data array using minimum Bayes Factor.
It can estimate necessary matrices using Maximum Likelihood if not provided.
}
\details{
The function first estimates mean and covariance matrices if they are not provided, then tests for outliers
using the \code{OutlierBFJ} function via Minimum Bayes Factor.
The resulting dataset contains the following columns:
\itemize{
\item Column Outdet is -1 if the observed matrix is considered an outlier.
\item Column Hv reports the Minimum Bayes Factor.
\item Columns hlb and hub report the indeterminacy lower and upper bounds.
If \code{plot = TRUE}, it generates a plot of the outlier statistics.
}
}
\examples{
# Generate a Synthetic Dataset
Data<-GenData(5, 3, 6, Rcn = 3, Rrn = 2)
# Run Outlier Detection
results <- OutlierDetection(Data$X3D,  M = Data$M, Sl = Data$Sl, V =  Data$V)
}
