####
odd.ratio = function(t, y){
  p1 = mean(y[t==1])
  p0 = mean(y[t==0])
  val = (p1*(1-p0))/(p0*(1-p1))
  return(val)
}

####
## Data Generating Process
N=2000

t = rbinom(N, 1, prob=0.3)
y1 = rbinom(N, 1, prob=0.25)
y0 = rbinom(N, 1, prob=0.25)

y = t*y1 + (1-t)*y0

####
# recall bias
nsim=10000

eta1.seq = seq(0, 0.5, by=0.02)
val.mat = matrix(NA, nrow = nsim, ncol = length(eta1.seq))
for(k in 1:nsim){
  
  
  t = rbinom(N, 1, prob=0.3)
  y1 = rbinom(N, 1, prob=0.25)
  y0 = rbinom(N, 1, prob=0.25)
  
  y = t*y1 + (1-t)*y0
  
  for(j in 1:length(eta1.seq)){
    eta1 = eta1.seq[j]
    eta0 = 0
    
    rand1 = rbinom(N, 1, prob=eta1)
    rand0 = rbinom(N, 1, prob=eta0)
    
    tstar = t + (1-t)*y*rand1 + (1-t)*(1-y)*rand0
    
    val.mat[k,j] = odd.ratio(tstar, y)
  }
  
  if(k%%100==0)cat("..", k)
}

mean.val = apply(val.mat, 2, mean)
sd.val = apply(log(val.mat), 2, sd)
appx.lower.val = exp(log(mean.val)-1.96*sd.val)
appx.upper.val = exp(log(mean.val)+1.96*sd.val)


plot(eta1.seq, mean.val, type = "l", ylim = c(min(lower.int.val), max(upper.int.val)), xlab = expression(eta[1]), ylab = "Odds Ratio")
lines(eta1.seq, appx.lower.val, lty=2)
lines(eta1.seq, appx.upper.val, lty=2)
abline(h=1, col="red")


