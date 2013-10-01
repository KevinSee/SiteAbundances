######################################################################
# function to estimate abundance from 2 or more pass delpletion data #
######################################################################
ZippinDepl = function(data, rmInvalid=FALSE, n.pass=2){
# n.pass is number of removal passes
  if(n.pass==2) {
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
    if(rmInvalid==tot.catch) N.hat[p.hat<p.threshold] = N.hat.SE[p.hat<p.threshold] = NA
  }
  if(n.pass > 2) {
    k = n.pass
    # total catch
    tot.catch = rowSums(data)
    zip.sum = rep(0, nrow(data))
    for(i in 1:n.pass) zip.sum = zip.sum + (i-1) * data[,i]
    
    f = function(tot.catch, X, k) {
      n.all = seq(tot.catch, 100*tot.catch)
      test = 2
      while(test>1) {
        for(n in n.all) {
          test = (n+1) / (n - tot.catch + 1) * prod((k*n - X - tot.catch + 1 + (k-1:k)) / (k*n - X + 2 + (k-1:k)))
          if(test <= 1) break
        }
      }
      N = n.all[which(n.all==n)]
      if(tot.catch==0) N=0
      return(N)
    }
    
    N.hat = mapply(f, rowSums(data), zip.sum, MoreArgs=list(k=n.pass))
    
    my.R = zip.sum / tot.catch
    f = function(R, q, k) R - (q / (1-q) + ((k*q^k) / (1-q^k)))
    q.hat = vector('numeric', nrow(data))
    for(i in 1:nrow(data)) q.hat[i] = ifelse(is.na(R[i]), NA, uniroot(f=f, interval=c(0,1), R=R[i], k=n.pass)$root)
    q.hat = sapply(as.list(my.R), function(x){
      if(!is.na(x)) result = uniroot(f=f, interval=c(0,1), R=x, k=n.pass)$root
      if(is.na(x)) result = NA
      return(result)})
    
  }
  return(data.frame(N.hat = N.hat, p.hat=p.hat))
}
