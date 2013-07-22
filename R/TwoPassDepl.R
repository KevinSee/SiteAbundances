##############################################################
# function to estimate abundance from 2 pass delpletion data #
##############################################################
TwoPassDepl = function(data, rmInvalid=FALSE, p.threshold=0.2){
  N.hat = data[,1]^2 / (data[,1] - data[,2])
  # if no fish caught on either pass, make N.hat 0
  N.hat[which(data[,1]==0 & data[,2]==0)] = 0
  p.hat = (data[,1] - data[,2]) / data[,1]
  # if the 2nd pass caught more fish than the first, make it invalid
  N.hat[which(data[,1]<=data[,2])] = p.hat[which(data[,1]<=data[,2])] = NA
  if(rmInvalid==T) N.hat[p.hat<p.threshold] = NA
  return(data.frame(N.hat = N.hat, p.hat=p.hat))
}
