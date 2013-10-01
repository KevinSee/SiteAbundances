##############################################################
# function to estimate abundance from k-pass delpletion data #
##############################################################
KPassDepl = function(data, n.pass=3){
#   k = n.pass
#   # total catch
#   tot.catch = rowSums(data)
#   # intermediate statistic
#   zip.sum = rep(0, nrow(data))
#   for(i in 1:n.pass) zip.sum = zip.sum + (i-1) * data[,i]
#   
#   lik.func = function(params, data, k) {
#     p = params[1]
#     N = params[2]
#     # total catch
#     tot.catch = rowSums(data)
#     # intermediate statistic
#     zip.sum = rep(0, nrow(data))
#     for(i in 1:k) zip.sum = zip.sum + (i-1) * data[i]
#     
#     lik = (factorial(N) * p^tot.catch * (1-p)^(k*N - zip.sum - tot.catch)) / (factorial(N - tot.catch) * prod(factorial(data)))
#     log.lik = -log(lik)
#     return(log.lik)
#   }
#   
#   optim(c(p=0.6, N=rowSums(data)+10), lik.func, data=data, k=3)
#     
#   }
#   
#   return(data.frame(N.hat = N.hat, N.hat.SE=N.hat.SE, p.hat=p.hat))
}
