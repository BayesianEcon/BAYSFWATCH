#' Generate Simulation Data for Bayesian Factor Analysis
#'
#' This function generates a 3D dataset (5x7x100) for Bayes Factor (BF) outlier detection by simulating matrix-variate normal distributions and injecting outliers based on different scenarios.
#'
#' @param outv A scalar value representing the amount to add as outliers.
#' @param outc An integer (1 to 4) defining the type of outlier injection:
#'   - 1: Adds `outv` to all values in the matrix.
#'   - 2: Adds `outv` to a single element (1,1) in the matrix.
#'   - 3: Adds `outv` to `outn` randomly chosen elements in the matrix.
#'   - 4: Adds `outv` to `Rcn` columns and `Rrn` rows, randomly selected.
#' @param outn The number of outliers to inject (only used in option 3). Must be lower than `p * n = 35`.
#' @param Rcn The number of columns to contain outliers (used in option 4). Default is 1.
#' @param Rrn The number of rows to contain outliers (used in option 4). Default is 1.
#'
#' @return A list containing:
#'   - `X3D0`: The generated 3D dataset without outliers.
#'   - `X3D`: The 3D dataset with outliers.
#'   - `M`: The mean matrix.
#'   - `Sl`: The row covariance matrix.
#'   - `Sp`: The scaled row covariance matrix.
#'   - `V`: The column covariance matrix.
#'
#' @details
#' The function first generates a dataset (a 5x7x100 array) from a matrix-variate normal distribution using `mvtnorm::rmvnorm`.
#' An outlier matrix is injected at the 80th position.
#' Depending on `outc`, the function injects outliers in different ways, modifying elements at specific positions or randomly distributing them.
#' The mean matrix is generated randomly and the factor controlling variance reduction in the prior is set to 0.01.
#'
#' @examples
#' set.seed(123)
#' data <- GenData(outv = 1.5, outc = 3, outn = 10, outp = 80)
#'
#' @export
GenData <- function(outv, outc, outn, Rcn = 1, Rrn = 1){

  #Generate simulation data for BF
  # @param Ts The number of time steps (or matrices) in the generated 3D dataset. Default is 100. #remove
  # @param p The number of rows in each matrix. Default is 5.  #remove
  # @param n The number of columns in each matrix. Default is 7.  #remove
  # @param M0 An integer: 1 generates the mean matrix randomly, 0 sets it to zero. Default is 1.  #remove
  # @param vv0 A factor controlling variance reduction in the prior. If 0, variance scales with `Ts`, otherwise uses `vv0`. Default is 0.1.  #remove
  # @param outp A vector specifying which matrices in the dataset should receive outliers.

  p = 5
  n = 7
  Ts = 100
  M0 = 1
  vv0 = 0.1
  outp = 80

  if (outn > p * n) {
    stop("The number of outliers must be lower than p * n")
  }

  M <- M0 * matrix(rnorm(p * n), p, n)
  tau <- matrix(rnorm(p * p), p, p)
  Sl <- tau %*% t(tau)
  vv <- ifelse(vv0 > 0, vv0, Ts)
  Sp <- Sl / vv
  gamma <- matrix(rnorm(n * n), n, n)
  V <- gamma %*% t(gamma)

  Mv <- as.vector(M)
  Sv <- kronecker(V, Sl)
  X3D <- array(NA, dim = c(p, n, Ts))

  for (i in 1:Ts) {
    X3D[,,i] <- matrix(mvtnorm::rmvnorm(1, mean = Mv, sigma = Sv), p, n)
  }

  X3D0 <- X3D  # Store data without outliers

  for (i in seq_along(outv)) {
    if (outc == 1) {
      X3D[,,outp[i]] <- outv[i] + X3D[,,outp[i]]
    } else if (outc == 2) {
      X3D[1,1,outp[i]] <- outv[i] + X3D[1,1,outp[i]]
    } else if (outc == 3) {
      outMp <- sample(1:(p * n), outn)
      outMs <- matrix(0, p, n)
      for (r in 1:p) {
        indices <- which((r - 1) * n < outMp & outMp <= r * n)
        outMs[r, outMp[indices] - (r - 1) * n] <- 1
      }
      X3D[,,outp[i]] <- outv[i] * outMs + X3D[,,outp[i]]
    } else if (outc == 4) {
      Rr <- sample(1:p, Rrn)
      Rc <- sample(1:n, Rcn)
      outMs <- matrix(0, p, n)
      outMs[Rr, Rc] <- 1
      X3D[,,outp[i]] <- outv[i] * outMs + X3D[,,outp[i]]
    }
  }

  return(list(X3D0 = X3D0, X3D = X3D, M = M, Sl = Sl, Sp = Sp, V = V))
}
