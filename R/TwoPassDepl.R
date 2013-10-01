##############################################################
# function to estimate abundance from 2 pass delpletion data #
##############################################################
TwoPassDepl = function(data, rmInvalid=FALSE, p.threshold=0.2){
  y1 = data[,1]
  y2 = data[,2]
  N.hat = y1^2 / (y1 - y2)
  # if no fish caught on either pass, make N.hat 0
  N.hat[which(y1==0 & y2==0)] = 0
  p.hat = (y1 - y2) / y1
  N.hat.SE = sqrt((y1^2 * y2^2 * (y1 + y2)) / (y1 - y2)^4)
  # if the 2nd pass caught more fish than the first, make it invalid
  N.hat[which(y1<=y2)] = N.hat.SE[which(y1<=y2)] = p.hat[which(y1<=y2)] = NA
  # mark as invalid any estimate with p.hat < p.threshold
  if(rmInvalid==T) N.hat[p.hat<p.threshold] = N.hat.SE[p.hat<p.threshold] = NA
  return(data.frame(N.hat = N.hat, N.hat.SE=N.hat.SE, p.hat=p.hat))
}
