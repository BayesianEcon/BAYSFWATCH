
fHtpdfm <- function(h, a1, n, p, B, Ms, SH, Sl, V, km) {
  
  eig_SH <- eigen(SH)
  QS <- eig_SH$vectors
  LS <- diag(eig_SH$values)
  SHri <- QS %*% solve(sqrt(LS)) %*% t(QS)
  SHr <- QS %*% sqrt(LS) %*% t(QS)
  Sz <- SHri %*% Sl %*% SHri
  Sz <- 0.5 * (Sz + t(Sz))  # Ensure symmetry
  
  TH <- SHr %*% solve(Sl) %*% (B - Ms) %*% solve(V) %*% t(B - Ms) %*% SHri
  TH <- 0.5 * (TH + t(TH))  # Ensure symmetry
  OM <- TH / 2
  
  eig_Sz <- eigen(Sz)
  Q <- eig_Sz$vectors
  lv <- abs(eig_Sz$values)  # Ensure real values
  bv <- diag(t(Q) %*% OM %*% Q)
  la <- mean(lv)
  
  dk <- numeric(km)
  for (k in 1:km) {
    dk[k] <- (n / (2 * k)) * sum((1 - la / lv)^k) + la * sum((bv / lv) * (1 - la / lv)^(k - 1))
  }
  
  f0 <- 1
  fk <- numeric(km)
  fk[1] <- dk[1]
  for (k in 2:km) {
    fka <- sum((1:(k - 1)) * dk[1:(k - 1)] * fk[(k - 1):1])
    fka <- fka + k * dk[k] * f0
    fk[k] <- fka / k
  }
  
  c0 <- exp(-sum(bv)) * prod((lv / la)^(-n / 2)) * f0
  ck <- numeric(km)
  for (k in 1:km) {
    ck[k] <- exp(-sum(bv)) * prod((lv / la)^(-n / 2)) * fk[k]
  }
  
  Htpdf <- 2 * c0 * dgamma(-2 * log(h / a1), shape = (n * p / 2) + 0, scale = 2 * la) / h
  
  for (k in 1:km) {
    Htpdf <- Htpdf + (2 * ck[k]) * dgamma(-2 * log(h / a1), shape = (n * p / 2) + k, scale = 2 * la) / h
  }
  
  return(Htpdf)
}
