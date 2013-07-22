###############################################################
# modified Lincoln-Petersen (Chapman) estimator
# to estimate abundance from single mark / single recap mark recapture dataa
###############################################################
ChapmanMR = function(data, rmInvalid=TRUE){
  if(is.data.frame(data)==F) data = as.data.frame(data)
  names(data) = c('M', 'C', 'R')
  N.hat = with(data, ((M + 1) * (C + 1)) / (R + 1) -1)
  N.hat.SE = with(data, sqrt(((M + 1) * (C + 1) * (M - R) * (C - R)) / ((R + 1)^2 * (R + 2))))
  p.hat = with(data, R / M)
  # using Robson & Regier criteria for valid abundance estimates
  Valid = with(data, (M * C)) > (N.hat * 4)
# Valid = (data$M * data$C) > (data$N.hat * 4) | data$R >= 7
  if(rmInvalid==T){
    N.hat[Valid==F] = NA
    N.hat.SE[Valid==F] = NA
  }
  return(data.frame(N.hat=N.hat, N.hat.SE=N.hat.SE, p.hat=p.hat, Valid=Valid))
}