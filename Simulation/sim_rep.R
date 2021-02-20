### 
source("../Rfunctions/basic_functions.R")

# library(optmatch)

## Number of repetitions
nsim = 10

rec.OR.mat = matrix(NA, nrow = nsim, ncol = 6)
for(k in 1:nsim){
  
  ## Data - generating 
  n = 800
  
  sub.x = rbind(c(0,0,0,0), c(0,0,0,1), c(0,0,1,0), c(0,0,1,1),
                c(0,1,0,0), c(0,1,0,1), c(0,1,1,0), c(0,1,1,1),
                c(1,0,0,0), c(1,0,0,1), c(1,0,1,0), c(1,0,1,1),
                c(1,1,0,0), c(1,1,0,1), c(1,1,1,0), c(1,1,1,1))
  rep.sub.x = do.call(rbind, replicate(n/16, sub.x, simplify=FALSE))
  x1 = rep.sub.x[,1]
  x2 = rep.sub.x[,2]
  x3 = rep.sub.x[,3]
  x4 = rep.sub.x[,4]
  
  x0 = rep(1, n)
  x = cbind(x0, rep.sub.x)
  
  beta.vec = c(-1, log(2), log(2), log(2), 0)
  gamma.vec0 = c(-2, log(2), log(2), 0, log(2))
  beta.condi = 0.5
  gamma.vec1 = c(-2 + beta.condi, log(2), log(2), 0, log(2))
  
  logit.ps = x %*% beta.vec
  logit.os0 = x %*% gamma.vec0
  logit.os1 = x %*% gamma.vec1
  
  ps = exp(logit.ps)/(1 + exp(logit.ps))
  os0 = exp(logit.os0)/(1 + exp(logit.os0))
  os1 = exp(logit.os1)/(1 + exp(logit.os1))
  
  t = rbinom(n, 1, prob = ps)
  y0 = rbinom(n, 1, prob = os0)
  y1 = rbinom(n, 1, prob = os1)
  y = t*y1 + (1-t)*y0
  
  
  # 
  mean(t); mean(y0); mean(y1); (mean(y1) * (1-mean(y0)))/(mean(y0) * (1-mean(y1)))
  # 0.509, 0.291, 0.395, 1.591
  
  exp(beta.condi)
  # 1.649
  
  ##########
  # Recall bias
  eta0 = 0.1
  eta1 = 0.1
  
  rand1 = rbinom(n, 1, prob=eta1)
  rand0 = rbinom(n, 1, prob=eta0)
  
  tstar = t + (1-t)*y*rand1 + (1-t)*(1-y)*rand0
  
  #################
  #### observable data
  obs.data = cbind(y, tstar, x1, x2, x3, x4)
  ## True
  true.y1 = mean(y1)
  true.y0 = mean(y0)
  true.OR = (true.y1 * (1-true.y0))/(true.y0 * (1-true.y1))
  
  
  ## 0. Crude
  crude.y1 = mean(y[tstar == 1])
  crude.y0 = mean(y[tstar == 0])
  crude.OR = (crude.y1 * (1-crude.y0))/(crude.y0 * (1-crude.y1))
  
  ## 1. ML
  m1 = ML.logistic.over(data = obs.data, eta = c(0.1,0.1))
  # m1$odds.ratio
  
  ## 2. Stratification 
  obs.data = as.data.frame(obs.data)
  
  prog.score.model = glm(y ~ x1 + x2 + x3 + x4, data = obs.data[obs.data$tstar == 0, ], family = "binomial")
  prog.score = as.vector(x %*% prog.score.model$coefficients)
  
  prop.score.model = glm(tstar ~ x1 + x2 + x3 + x4, data = obs.data, family = "binomial", x = TRUE)
  prop.score = as.vector(x %*% prop.score.model$coefficients)
  
  # 2.1. Propensity score stratification 
  
  prob.length = 6
  ps.seq = quantile(prop.score, prob = seq(0, 1, length.out = prob.length))
  
  astar.vec = bstar.vec = cstar.vec = dstar.vec = rep(NA, length(ps.seq)-1)
  for(i in 1:(length(ps.seq)-1)){
    lower = ps.seq[i]
    upper = ps.seq[i+1]
    sbg = obs.data[(prop.score > lower & prop.score <= upper),]
    
    if(i == 1){
      sbg = obs.data[(prop.score >= lower & prop.score <= upper),]
    }
    
    astar.vec[i] =  sum(sbg$tstar == 1 & sbg$y == 1)
    bstar.vec[i] =  sum(sbg$tstar == 1 & sbg$y == 0)
    cstar.vec[i] =  sum(sbg$tstar == 0 & sbg$y == 1)
    dstar.vec[i] =  sum(sbg$tstar == 0 & sbg$y == 0)
  }
  count.prop = cbind(astar.vec, bstar.vec, cstar.vec, dstar.vec)
  
  # 2.2. Prognostic score startification 
  prob.length = 6
  ps.seq = quantile(prog.score, prob = seq(0, 1, length.out = prob.length))
  
  astar.vec = bstar.vec = cstar.vec = dstar.vec = rep(NA, length(ps.seq)-1)
  for(i in 1:(length(ps.seq)-1)){
    lower = ps.seq[i]
    upper = ps.seq[i+1]
    sbg = obs.data[(prog.score > lower & prog.score <= upper),]
    
    if(i == 1){
      sbg = obs.data[(prog.score >= lower & prog.score <= upper),]
    }
    astar.vec[i] =  sum(sbg$tstar == 1 & sbg$y == 1)
    bstar.vec[i] =  sum(sbg$tstar == 1 & sbg$y == 0)
    cstar.vec[i] =  sum(sbg$tstar == 0 & sbg$y == 1)
    dstar.vec[i] =  sum(sbg$tstar == 0 & sbg$y == 0)
  }
  count.prog = cbind(astar.vec, bstar.vec, cstar.vec, dstar.vec)
  
  m21 = SP.inference.over(count.mat = count.prop, eta = c(0.1, 0.1))
  m22 = SP.inference.over(count.mat = count.prog, eta = c(0.1, 0.1))
  
  #### 3. Full matching - This is not shown in the main manuscript, but discussed in supplementary materials
  treatment = obs.data$y
  Xmat = prop.score.model$x[,-1]
  
  # Rank based Mahalanobis distance
  distmat=smahal(treatment,Xmat)
  distmat2 = distmat
  
  ### Create a subject index and name the rows and columns of distance matrix by ### this subject index
  subject.index=seq(1,length(treatment),1)
  rownames(distmat2)=subject.index[treatment==1]
  colnames(distmat2)=subject.index[treatment==0]
  treated=treatment
  
  ## Exactly matched on "female", "marital.prob", "rural".
  maxadd=max(distmat2)+100
  treated.indices=subject.index[treatment==1]
  control.indices=subject.index[treatment==0]
  for(i in 1:nrow(distmat2)){
    distmat2[i,]=distmat2[i,]+maxadd*(obs.data$x1[control.indices] != obs.data$x1[treated.indices[i]])
    distmat2[i,]=distmat2[i,]+maxadd*(obs.data$x2[control.indices] != obs.data$x2[treated.indices[i]])
    distmat2[i,]=distmat2[i,]+maxadd*(obs.data$x3[control.indices] != obs.data$x3[treated.indices[i]])
    distmat2[i,]=distmat2[i,]+maxadd*(obs.data$x4[control.indices] != obs.data$x4[treated.indices[i]])
  }
  
  matchvec = suppressWarnings(fullmatch(distmat2))
  
  
  # Create vectors of the subject indices of the treatment units ordered by
  # their matched set and corresponding control unit
  
  subject.indices=as.numeric(names(matchvec))
  unique.match = unique(matchvec)[-(sum(treated==1)+1)]
  
  num.matched.pair = length(unique.match)
  
  astar.vec = bstar.vec= cstar.vec = dstar.vec = n.vec = rep(NA, num.matched.pair)
  for(i in 1:num.matched.pair){
    temp.indices = which(matchvec == unique.match[i])
    temp.subject.indices = subject.indices[temp.indices]
    matched.treated = temp.subject.indices[which(treatment[temp.subject.indices]==1)]
    matched.control = temp.subject.indices[which(treatment[temp.subject.indices]==0)]
    
    temp.abuse.case = obs.data$tstar[matched.treated]
    temp.abuse.control = obs.data$tstar[matched.control]
    
    astar.vec[i] = sum(temp.abuse.case==1)
    cstar.vec[i] = sum(temp.abuse.case==0)
    
    bstar.vec[i] = sum(temp.abuse.control == 1)
    dstar.vec[i] = sum(temp.abuse.control == 0)
    n.vec[i] = astar.vec[i] + bstar.vec[i] + cstar.vec[i] + dstar.vec[i]
  }
  count.matching = cbind(astar.vec, bstar.vec, cstar.vec, dstar.vec)
  
  m3 = MH.inference.over(count.mat = count.matching, eta = c(0.1, 0.1))
  
  
  rec.OR.mat[k,] = c(true.OR, crude.OR, m1$odds.ratio, m21$est.OR, m22$est.OR, m3$est.OR)

}

round(rec.OR.mat, 2)
apply(rec.OR.mat, 2, mean)
