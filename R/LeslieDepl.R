###############################################################
# Leslie estimator
# to estimate abundance from multiple pass delpletion data
###############################################################
LeslieDepl = function(data, site.spec.p = FALSE, Ricker.correction=TRUE)
  # pass in data frame with catches for each pass across each row
{
  n.pass = ncol(data)
  # manipulate data to be in the form of catch for each pass, and cumulative catch on all previous passes
  long.df = data.frame(site.num = factor(rep(1:nrow(data), each=n.pass)), pass.num = rep(1:n.pass, nrow(data)), C = c(t(data)))
  long.df$K = with(long.df, stack(tapply(C, site.num, cumsum))[,1] - C)
  
  if(Ricker.correction==T) long.df$K = long.df$K + long.df$C/2
  
  if(site.spec.p==F) # fit Leslie model, assuming common catchability
  {	
    if(nlevels(long.df$site.num)>1) {
      mod1 = lm(C ~ -1 + site.num + K, data=long.df)
      p.hat = -coef(mod1)[which(names(coef(mod1))=="K")]
      N.hat = coef(mod1)[1:(which(names(coef(mod1))=="K")-1)] / p.hat
      N.hat.SE = summary(mod1)$sigma / p.hat * sqrt((1/nrow(data)) + ((N.hat - with(long.df, tapply(K, site.num, mean))) / (nrow(data)-1) * with(long.df, tapply(K, site.num, var))))
    }
    if(nlevels(long.df$site.num)==1) {
      mod1 = lm(C ~ K, data=long.df)
      p.hat = -coef(mod1)[which(names(coef(mod1))=="K")]
      N.hat = coef(mod1)[1] / p.hat
      N.hat.SE = summary(mod1)$sigma / p.hat * sqrt((1/nrow(data)) + ((N.hat - with(long.df, tapply(K, site.num, mean))) / (nrow(data)-1) * with(long.df, tapply(K, site.num, var))))
    }
  }
  
#   if(site.spec.p==T) # fit a model with site specific catchability coefficients
#   {
#     mod2 = lm(C ~ -1 + site.num * K, data=long.df)
#     N.hat = coef(mod2)[grep('K', names(coef(mod2)), invert=T)]
#     p.hat = -c(coef(mod2)['K'] + coef(mod2)[grep(':K', names(coef(mod2)))])
#   }
  if(site.spec.p==T) # fit a model with site specific catchability coefficients
  {
    mod2 = lmer(C ~ -1 + site.num + K + (-1 + K | site.num), data=long.df)
    N.hat = fixef(mod2)[grepl('site.num', names(fixef(mod2)))]
    p.hat = -c(fixef(mod2)["K"] + ranef(mod2)$site.num[,1])
    # not working yet
#     N.hat.SE = summary(mod2)@sigma / p.hat * sqrt((1/nrow(data)) + ((N.hat - with(long.df, tapply(K, site.num, mean))) / (nrow(data)-1) * with(long.df, tapply(K, site.num, var))))
    N.hat.SE = NA
  }
  return(data.frame(N.hat=N.hat, N.hat.SE=N.hat.SE, p.hat=p.hat))
} 
