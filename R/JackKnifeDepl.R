###############################################################
# jacknife estimator from Hankin & Mohr (2001)
# to estimate abundance from 3 pass delpletion data
###############################################################
JackKnifeDepl = function(data, group = NULL, grouped=FALSE)
{
  names(data) = c('p1', 'p2', 'p3')[1:ncol(data)]
  if(grouped) data$grp = group
  # set up vectors to capture results
  results = data.frame(N.hat = rep(NA, nrow(data)), N.hat.SE=NA, p.hat=NA, row.names=rownames(data))
  # which rows having missing data?
  miss.data = sort(unique(which(is.na(data), T)[,1]))
  # remove missing data rows
  data = na.omit(data)
  if(grouped==F) {
    if(length(miss.data)>0) {
      results$p.hat[-miss.data] = 1 - (sum(rowSums(data[,c('p1', 'p2', 'p3')])) - sum(data$p1)) / (sum(rowSums(data[,c('p1', 'p2', 'p3')])) - sum(data$p3))
      results$N.hat[-miss.data] = rowSums(data[,c('p1', 'p2')]) + data$p3/ results$p.hat[-miss.data]
      results$N.hat.SE[-miss.data] = 3*2*data$p3
    }
    if(length(miss.data)==0) {
      results$p.hat = 1 - (sum(rowSums(data[,c('p1', 'p2', 'p3')])) - sum(data$p1)) / (sum(rowSums(data[,c('p1', 'p2', 'p3')])) - sum(data$p3))
      results$N.hat = rowSums(data[,c('p1', 'p2')]) + data$p3/ results$p.hat
      results$N.hat.SE = 3*2*data$p3
    }  
  }
  if(grouped==T) {
    for(group in unique(data$grp)) {
      if(sum(rowSums(data[data$grp==group,c('p1', 'p2', 'p3')]))==0) {
        p.hat=NA
        N.est.SE=NA
        N.est = 0
        valid = T
      } 
      if(sum(rowSums(data[data$grp==group,c('p1', 'p2', 'p3')]))>0) {
        p.hat = 1 - (sum(rowSums(data[data$grp==group,c('p1', 'p2', 'p3')])) - sum(data$p1[data$grp==group])) / (sum(rowSums(data[data$grp==group,c('p1', 'p2', 'p3')])) - sum(data$p3[data$grp==group]))
        N.est = data$p1[data$grp==group] + data$p2[data$grp==group] + data$p3[data$grp==group]/p.hat
        N.est.SE = 3*2*data$p3[data$grp==group]
        valid=T
        if(N.est==Inf) {
        	N.est = rowSums(data[data$grp==group,c('p1', 'p2', 'p3')])
        	N.est.SE = NA
        	p.hat = NA
        	valid=F
        }
      }
        results[match(rownames(subset(data, grp==group)), rownames(results)),] = data.frame(N.hat = N.est, N.hat.SE=N.est.SE, p.hat=p.hat, Valid=valid)
    }
  }
  return(results)
}
