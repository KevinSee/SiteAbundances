###############################################################
# estimator from Guy & Brown (2007)
# to estimate abundance from 3 pass delpletion data
###############################################################
GuyBrownDepl = function(data){
  x = 2 * (data[,1] + data[,2])
  y = rowSums(data, na.rm=T)
  # estimate abundance
  N.hat = (6*x^2 - 3*x*y - y^2 + y*sqrt(y^2 + 6*x*y - 3*x^2)) / (18*(x- y))
  # estimate probability of capture
  p.hat = (3*x - y -sqrt(y^2 + 6*x*y - 3*x^2)) / (2*x)
  # calculate standard error of N.hat
  N.hat.SE = sqrt((N.hat*(1-p.hat^3)*p.hat^3) / ((1-p.hat^3)^2 - (3*(1-p.hat))^2*p.hat^2))
  valid = ifelse(N.hat<0 | N.hat==Inf | N.hat==-Inf, F, T)
  return(data.frame(N.hat=N.hat, N.hat.SE=N.hat.SE, p.hat=p.hat, Valid=valid))
}
