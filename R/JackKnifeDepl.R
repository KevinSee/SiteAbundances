###############################################################
# jacknife estimator from Hankin & Mohr (2001)
# to estimate abundance from 3 pass delpletion data
###############################################################
JackKnifeDepl = function(data, group = NULL, grouped=FALSE)
{
  names(data) = c('p1', 'p2', 'p3')[1:ncol(data)]
  if(grouped) data$grp = group
  # which rows having missing data?
  miss.Dat = sort(unique(which(is.na(data), T)[,1]))
  N.hat = rep(NA, nrow(data))
  # remove missing data rows
  data = na.omit(data)
  if(grouped==F){
    p.hat = 1 - (sum(rowSums(data[,c('p1', 'p2', 'p3')])) - sum(data$p1)) / (sum(rowSums(data[,c('p1', 'p2', 'p3')])) - sum(data$p3))
    if(length(miss.Dat)>0) N.hat[-miss.Dat] = rowSums(data[,c('p1', 'p2')]) + data$p3/ p.hat
    if(length(miss.Dat)==0) N.hat = rowSums(data[,c('p1', 'p2')]) + data$p3/ p.hat
    results = data.frame(N.hat = N.hat, p.hat=p.hat)
  }
  if(grouped==T){
    for(group in unique(data$grp))
    {
      p.hat = 1 - (sum(rowSums(data[data$grp==group,c('p1', 'p2', 'p3')])) - sum(data$p1[data$grp==group])) / (sum(rowSums(simDat[data$grp==group,c('p1', 'p2', 'p3')])) - sum(data$p3[data$grp==group]))
      N.est = data$p1[data$grp==group] + data$p2[data$grp==group] + data$p3[data$grp==group]/p.hat
      if(group == unique(data$grp)[1]) results = data.frame(N.hat = N.hat, p.hat=p.hat)
      else results = rbind(results, data.frame(N.hat = N.est, p.hat=p.hat))
    }
  }
  return(results)
}
