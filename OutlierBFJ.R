library(MixMatrix)

OutlierBFJ <- function(X3D, M, Sl, V, RollW = TRUE, kk = 40, vv = 0.1, dtvs = 0.05, dte = 0.999, dts = 100, km = 100, hl = 0.5, qq = 0.01, step = 3000, estim = FALSE) {
  options(width = 80)

  X0 <- X3D
  p <- dim(X0)[1]
  n <- dim(X0)[2]
  T0 <- dim(X0)[3]

  Sp <- Sl / vv
  dtv <- seq(dtvs, dte, (dte-dtvs)/(dts-1))

  Htkv <- matrix(NA, nrow = kk, ncol = dts)
  a1vv <- matrix(NA, nrow = kk, ncol = dts)

  if (step > 0) {
    etij <- array(NA, dim = c(dts, step, kk, 5))
    etijt <- array(NA, dim = c(dts, kk, 7))
    htijpI <- matrix(NA, nrow = dts, ncol = kk)
    htijvI <- matrix(NA, nrow = dts, ncol = kk)
  }

  Mku <- array(NA, dim = c(p, n, kk))
  eee <- 1e-8

  pb <- txtProgressBar(min = 0, max = kk, style = 3)
  for (k in 1:kk) {

    ###loading bar

    kroll <- ifelse(RollW, k - 1, 0)

    Xk <- X0[, , (1 + kroll):(T0 - kk + (k - 1))]
    T <- dim(Xk)[3]

    if(estim == TRUE){
      Estim<-MixMatrix::MLmatrixnorm(Xk, row.mean = FALSE,  col.mean = FALSE)
      Sl <- Estim$U
      V <- Estim$V*Estim$var
      Sp <- Sl/vv}

    Y <- matrix(aperm(Xk, c(1, 3, 2)), ncol = n)

    if (k > 1) {
      M <- Mku[, , k - 1]
      Mku[, , k] <- M + (Sp %*% solve(Sl + Sp) %*% (X0[, , T0 - kk + k] - M))
    }

    Yt <- X0[, , T0 - kk + k]

    Xbar <- apply(Xk, c(1, 2), mean)  # (p x n) average over T

    A <- Sp
    B <- Sl / T + Sp

    Ms <- M + A %*% solve(B, Xbar - M)
    Ss <- A - A %*% solve(B, A)

    if (k == 1) Mku[, , k] <- Ms

    Sd <- Sl + Ss

    for (i in 1:dts) {

      dti <- dtv[i]
      SAdi <- Sl + Ss / dti
      SHi <- (dti / (1 - dti)) * SAdi %*% solve(Ss) %*% Sd
      Htkv[k, i] <- (det(SAdi) / det(Sd))^(n / 2) * exp(-0.5 * sum(diag(solve(SHi) %*% (Yt - Ms) %*% solve(V) %*% t(Yt - Ms))))
      a1vv[k, i] <- (det(SAdi) / det(Sd))^(n / 2)

      if (step > 0) {
        a1ui <- a1vv[k, i]
        htub <- a1vv[k, i]
        htui <- seq(eee, htub, (htub-eee)/(step+1))[2:(step+1)]
        htvl <- htui
        htvu <- htui

        # fHtcdfm_cpp(htvu[1], a1ui, n, p, Ms, Ms, SHi, Sd, V, km)
        # etij[i, , k, 1] <- 1-fHtcdfm(htvu, a1ui, n, p, Ms, Ms, SHi, Sd, V, km)
        # etij[i, , k, 2] <- fHtcdfm(htvl, a1ui, n, p, Ms, Ms, SHi, SAdi, V, km)
        etij[i, , k, 3] <- fHtcdfm(htvl, a1ui, n, p, Ms, Ms, SHi, Sd, V, km)
        # etij[i, , k, 4] <- 1 - fHtcdfm(htvu, a1ui, n, p, Ms, Ms, SHi, SAdi, V, km)
        # etij[i, , k, 5] <- etij[i, , k, 3] + etij[i, , k, 4]
        min_idx <- which.min(abs(etij[i, , k, 3] - qq))
        htijpI[i, k] <- min_idx
        htijvI[i, k] <- htui[min_idx]
        # etijt[i, k, 1] <- etij[i, min_idx, k, 1]
        # etijt[i, k, 2] <- etij[i, min_idx, k, 2]
        # etijt[i, k, 3] <- etij[i, min_idx, k, 3]
        # etijt[i, k, 4] <- etij[i, min_idx, k, 4]
        # etijt[i, k, 5] <- 1 - fHtcdfm_cpp(max(eee, 2 - htijvI[i, k]), a1ui, n, p, Ms, Ms, SHi, SAdi, V, km)
        # etijt[i, k, 6] <- etij[i, min_idx, k, 3] + etij[i, min_idx, k, 4]
        # etijt[i, k, 7] <- abs(etij[i,  min_idx ,  k, 1] - hl) + abs(etij[i, min_idx, k, 2] - (1 - hl))
      }
    }

    setTxtProgressBar(pb, k)
    flush.console()  # Force refresh on macOS
  }

  jj <- 0.75
  jjv <- floor(jj * dts)
  hlb <- htijvI[jjv, ]
  hub <- 2 - hlb
  Outdet <- ifelse(Htkv[, jjv] < hlb, -1, ifelse(Htkv[, jjv] > hub, 1, 0))
  Hv <- Htkv[, jjv]

  return(list(Outdet = Outdet, Hv = Hv, hlb = hlb, hub = hub))
}
