rewardtransformparm <- function(rewards, initprob, subintensemat){
  d <- sum(rewards > 0)
  p <- length(initprob)
  qmat <- matrix(rep(0,p^2), ncol = p)
  for(i in 1:p){
    for(j in (1:p)[-i]){
      qmat[i,j] <- -subintensemat[i,j]/subintensemat[i,i]
    }
  }
  ##
  # If all rewards are stricly postive everything is simpler
  ##
  if(d == p){
    pmat <- qmat
    alphavec <- initprob
  }
  else{
    qplusplus <- qmat[(rewards > 0),(rewards > 0)]
    qpluszero <- qmat[(rewards > 0),(rewards == 0)]
    qzeroplus <- qmat[(rewards == 0),(rewards > 0)]
    qzerozero <- qmat[(rewards == 0 ), (rewards == 0)]
    pmat <- qplusplus + qpluszero %*% solve(diag(1, nrow = p-d)-qzerozero) %*% qzeroplus
    piplus <- initprob[(rewards > 0)]
    pizero <- initprob[(rewards == 0)]
    alphavec <- piplus + pizero %*% solve(diag(1, nrow = p-d)-qzerozero) %*% qzeroplus 
    subintensemat <- as.matrix(subintensemat[(rewards > 0), (rewards >0)])
    rewards <- rewards[rewards > 0]
  }
  pvec <- 1 - rowSums(pmat)
  Tstarmat <- matrix(rep(0,d^2), ncol = d)
  tstarvec <- rep(0,d)
  for(i in 1:d){
    for(j in (1:d)[-i]){
      Tstarmat[i,j] <- -subintensemat[i,i]/rewards[i]*pmat[i,j]
    }
    tstarvec[i] <- -subintensemat[i,i]/rewards[i]*pvec[i]
    Tstarmat[i,i] <- -sum(Tstarmat[i,])-tstarvec[i]
  }
  list("newinitprob" = alphavec, "newsubintensitymatrix" = Tstarmat, "defect" = 1 - sum(alphavec))
}
