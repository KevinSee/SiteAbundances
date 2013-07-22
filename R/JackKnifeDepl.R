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
  # remove missing data rows
  data = na.omit(data)
  if(grouped==F){
    N.est = p.hat = N.est.SE = rep(NA, nrow(data))
    p.hat = 1 - (sum(rowSums(data[,c('p1', 'p2', 'p3')])) - sum(data$p1)) / (sum(rowSums(data[,c('p1', 'p2', 'p3')])) - sum(data$p3))
    if(length(miss.Dat)>0){
      N.est[-miss.Dat] = rowSums(data[,c('p1', 'p2')]) + data$p3/ p.hat
      N.est.SE[-miss.Dat] = 3*2*data$p3
    }
    if(length(miss.Dat)==0) {
      N.est = rowSums(data[,c('p1', 'p2')]) + data$p3/ p.hat
      N.est.SE = 3*2*data$p3
    }
    results = data.frame(N.hat = N.est, N.hat.SE=N.est.SE, p.hat=p.hat)
  }
  if(grouped==T){
    results = data.frame(N.hat = rep(NA, nrow(data)), N.hat.SE=NA, p.hat=NA, row.names=rownames(data))
    for(group in unique(data$grp))
    {
      if(sum(rowSums(data[data$grp==group,c('p1', 'p2', 'p3')]))==0){
        p.hat=NA
        N.est.SE=NA
        N.est = 0
      } 
      if(sum(rowSums(data[data$grp==group,c('p1', 'p2', 'p3')]))>0){
        p.hat = 1 - (sum(rowSums(data[data$grp==group,c('p1', 'p2', 'p3')])) - sum(data$p1[data$grp==group])) / (sum(rowSums(data[data$grp==group,c('p1', 'p2', 'p3')])) - sum(data$p3[data$grp==group]))
        N.est = data$p1[data$grp==group] + data$p2[data$grp==group] + data$p3[data$grp==group]/p.hat
        N.est.SE = 3*2*data$p3[data$grp==group]
      }
        results[match(rownames(subset(data, grp==group)), rownames(results)),] = data.frame(N.hat = N.est, N.hat.SE=N.est.SE, p.hat=p.hat)
    }
  }
  return(results)
}
