########################################################
### Application of the Proposed Methods 
########################################################
source("../Rfunctions/basic_functions.R")

## You need to create "wls_anger.csv" first. 
## Use WLS_to_dataset.R to create this csv file. 
dataset = read.csv("../Rcode/wls_anger.csv") 

########################################################
### Data Construction
########################################################

## Outcome variable: anger score
## We use the 90th percentile (18 points) to dichotomize the anger score
dataset$anger90 = (dataset$anger >= 18)

## Exposure variable: childhood abuse
## We use a binary abuse exposure variable 
## Abused if there was abuse by either father or mother ("some" or "a lot" response) 
dataset$abuse = (dataset$abusefa >= 2 | dataset$abusemo >= 2)

## There are 7 covariates - (to adjust for confounders)
## (1) Sex, (2) age at the time of the interview, (3) father's education, 
## (4) mother's education, (5) parental income, (6) farm background, 
## (7) an indicator of parents' marital problems or single parent. 

## Make (1), (3), (4), (6), (7) binary variables
dataset$female = (dataset$sex == 2) # (1) sex
dataset$hsd.fa = (dataset$edufa >= 12) # (3) father's education - high school diploma
dataset$hsd.mo = (dataset$edumo >= 12) # (4) mother's education - high school diploma
dataset$rural = (dataset$farmback == 1) # (6) farm background - rural or not
dataset$marital.prob = (dataset$problems <= 1) # (7) marital problem

# Also, make age discrete with 3 categories. 
dataset$age.cat = rep(2, length(dataset$age))
dataset$age.cat[dataset$age <= 50]=1 
dataset$age.cat[dataset$age >= 57]=3
# dataset$age.cat = as.factor(dataset$age.cat)

# Finally, use a log transformation on (5)
dataset$log.income = log(dataset$income)


## Now, we finally choose the variables to use. 
## !! We set 
## 1st column: outcome
## 2nd column: exposure
## 3rd to the last column: covariates
newdata = dataset[, c("anger90", "abuse", "female", "age", "hsd.fa", "hsd.mo", "marital.prob", "log.income", "rural")]

########################################################
### Data Analysis: Maximum Likelihood method
## Since exposure is assumed to be under-reported, 
## We use "ML.logistic.under" 

## Without accounting for recall bias
ML.est = ML.logistic.under(data = newdata, zeta = c(0,0))
ML.est$est.OR


########################################################
### Data Analysis: Prognostic Stratification method

### Estimate the prognostic score assuming no effect modification
prog.score.model = glm(anger90 ~ abuse + female + age + hsd.fa + hsd.mo + log.income + rural + marital.prob, data = newdata, family = "binomial", x = TRUE)
totX = prog.score.model$x

prog.score = as.vector(totX[,-2] %*% prog.score.model$coefficients[-2])

### Make 5 strata based on prog.score
prob.length = 6
ps.seq = quantile(prog.score, prob = seq(0, 1, length.out = prob.length))

astar.vec = bstar.vec = cstar.vec = dstar.vec = rep(NA, length(ps.seq)-1)
for(i in 1:(length(ps.seq)-1)){
  lower = ps.seq[i]
  upper = ps.seq[i+1]
  sbg = newdata[(prog.score > lower & prog.score <= upper),]
  
  if(i == 1){
    sbg = newdata[(prog.score >= lower & prog.score <= upper),]
  }
  
  astar.vec[i] = sum(sbg$abuse == 1 & sbg$anger90 == 1)
  bstar.vec[i] = sum(sbg$abuse == 1 & sbg$anger90 == 0)
  cstar.vec[i] = sum(sbg$abuse == 0 & sbg$anger90 == 1)
  dstar.vec[i] = sum(sbg$abuse == 0 & sbg$anger90 == 0)
}
n.vec = astar.vec + bstar.vec + cstar.vec + dstar.vec

count.mat = cbind(astar.vec, bstar.vec, cstar.vec, dstar.vec)
SP.est = SP.inference.under(count.mat = count.mat, zeta = c(0,0))
SP.est$est.OR


########################################################
### Table 3
########################################################
## Recovered Marginal COR
round(ML.est$est.OR, 2)
round(SP.est$est.OR, 2)


## Parallel computing
library(parallel)
library(doParallel)
library(foreach)

cl <- parallel::makeCluster(10) # need to specify the number of cores to use. 
doParallel::registerDoParallel(cl)

time.before = Sys.time();
nsim = 10
boot.mat = foreach(k = 1:nsim, .combine = "rbind") %dopar% {
  N = dim(newdata)[1]
  sample.index = sample(1:N, N, replace = T)
  boot.data = newdata[sample.index, ]
  
  ## Bootstrap ML
  boot.ML.est = ML.logistic.under(data = boot.data, zeta = c(0,0))$est.OR
  
  ## Bootstrap Stratification
  boot.prog.score.model = glm(anger90 ~ abuse + female + age + hsd.fa + hsd.mo + log.income + rural + marital.prob, data = boot.data, family = "binomial", x = TRUE)
  boot.totX = boot.prog.score.model$x
  boot.prog.score = as.vector(boot.totX[,-2] %*% boot.prog.score.model$coefficients[-2])
  
  ### Make 5 strata based on boot.prog.score
  prob.length = 6
  ps.seq = quantile(boot.prog.score, prob = seq(0, 1, length.out = prob.length))
  
  astar.vec = bstar.vec = cstar.vec = dstar.vec = rep(NA, length(ps.seq)-1)
  for(i in 1:(length(ps.seq)-1)){
    lower = ps.seq[i]
    upper = ps.seq[i+1]
    sbg = boot.data[(boot.prog.score > lower & boot.prog.score <= upper),]
    
    if(i == 1){
      sbg = boot.data[(boot.prog.score >= lower & boot.prog.score <= upper),]
    }
    
    astar.vec[i] = sum(sbg$abuse == 1 & sbg$anger90 == 1)
    bstar.vec[i] = sum(sbg$abuse == 1 & sbg$anger90 == 0)
    cstar.vec[i] = sum(sbg$abuse == 0 & sbg$anger90 == 1)
    dstar.vec[i] = sum(sbg$abuse == 0 & sbg$anger90 == 0)
  }
  boot.count.mat = cbind(astar.vec, bstar.vec, cstar.vec, dstar.vec)
  boot.SP.est = SP.inference.under(count.mat = boot.count.mat, zeta = c(0,0))$est.OR
  c(boot.ML.est, boot.SP.est)
}
time.after = Sys.time();
time.after - time.before

parallel::stopCluster(cl)

## Bootstrap Estimation of the Variance of 
sd(log(boot.mat[,1]))
sd(log(boot.mat[,2]))

## 95% Confidence Intervals 
# 95% CI for ML
ML.ci95 = c(exp(log(ML.est$est.OR) - 1.96*sd(log(boot.mat[,1]))), exp(log(ML.est$est.OR) + 1.96*sd(log(boot.mat[,1]))))

# 95% CI for ML
SP.ci95 = c(exp(log(SP.est$est.OR) - 1.96*sd(log(boot.mat[,2]))), exp(log(SP.est$est.OR) + 1.96*sd(log(boot.mat[,2]))))

round(ML.ci95, 2); round(SP.ci95, 2)


########################################################
### Table 4 & Figure 2
########################################################
## Table 4 shows the five values of zeta
## To replicate Figure 2, you may want to use a finer vector for zeta. 
## You may want to use "replication.R" to generate the same table and figure. 

zeta.seq = seq(0, 0.5, by = 0.1) 

nsim = 500

sensi.mat = matrix(NA, nrow = length(zeta.seq), ncol = 8)

time.before = Sys.time();

cl <- parallel::makeCluster(10) # need to specify the number of cores to use. 
doParallel::registerDoParallel(cl)

boot.mat = foreach(k = 1:nsim, .combine = "rbind") %dopar% {
  N = dim(newdata)[1]
  sample.index = sample(1:N, N, replace = T)
  boot.data = newdata[sample.index, ]
  
  ## Bootstrap Stratification
  boot.prog.score.model = glm(anger90 ~ abuse + female + age + hsd.fa + hsd.mo + log.income + rural + marital.prob, data = boot.data, family = "binomial", x = TRUE)
  boot.totX = boot.prog.score.model$x
  boot.prog.score = as.vector(boot.totX[,-2] %*% boot.prog.score.model$coefficients[-2])
  
  ### Make 5 strata based on boot.prog.score
  prob.length = 6
  ps.seq = quantile(boot.prog.score, prob = seq(0, 1, length.out = prob.length))
  
  astar.vec = bstar.vec = cstar.vec = dstar.vec = rep(NA, length(ps.seq)-1)
  for(i in 1:(length(ps.seq)-1)){
    lower = ps.seq[i]
    upper = ps.seq[i+1]
    sbg = boot.data[(boot.prog.score > lower & boot.prog.score <= upper),]
    
    if(i == 1){
      sbg = boot.data[(boot.prog.score >= lower & boot.prog.score <= upper),]
    }
    
    astar.vec[i] = sum(sbg$abuse == 1 & sbg$anger90 == 1)
    bstar.vec[i] = sum(sbg$abuse == 1 & sbg$anger90 == 0)
    cstar.vec[i] = sum(sbg$abuse == 0 & sbg$anger90 == 1)
    dstar.vec[i] = sum(sbg$abuse == 0 & sbg$anger90 == 0)
  }
  boot.count.mat = cbind(astar.vec, bstar.vec, cstar.vec, dstar.vec)
  
  boot.ML.est = boot.SP.est = rep(NA, length(zeta.seq))
  for(j in 1:length(zeta.seq)){
    zeta.val = zeta.seq[j]
    
    ## Bootstrap ML
    boot.ML.est[j] = ML.logistic.under(data = boot.data, zeta = c(zeta.val, zeta.val))$est.OR
    boot.SP.est[j] = SP.inference.under(count.mat = boot.count.mat, zeta = c(zeta.val, zeta.val))$est.OR
  }
  c(boot.ML.est, boot.SP.est)
}
parallel::stopCluster(cl)


time.after = Sys.time();
time.after - time.before


boot.mat.ML = boot.mat[,1:length(zeta.seq)]
boot.mat.SP = boot.mat[,(length(zeta.seq)+1):(2*length(zeta.seq))]

ML.sd = apply(log(boot.mat.ML), 2, sd)
S.sd = apply(log(boot.mat.SP), 2, sd)


#### 
time.before = Sys.time();

cl <- parallel::makeCluster(10) # need to specify the number of cores to use. 
doParallel::registerDoParallel(cl)

est.vec = foreach(j = 1:length(zeta.seq), .combine = "rbind") %dopar% {
  zeta.val = zeta.seq[j]
  ML.est.sensi = ML.logistic.under(data = newdata, zeta = c(zeta.val, zeta.val))$est.OR
  SP.est.sensi = SP.inference.under(count.mat = count.mat, zeta = c(zeta.val, zeta.val))$est.OR
  c(ML.est.sensi, SP.est.sensi)
}
parallel::stopCluster(cl)

time.after = Sys.time();
time.after - time.before

#####
## Sensi.mat
sensi.mat = cbind(zeta.seq, est.vec[,1], ML.sd, exp(log(est.vec[,1]) - 1.96*ML.sd), exp(log(est.vec[,1]) + 1.96*ML.sd),
                  est.vec[,2], ML.sd, exp(log(est.vec[,2]) - 1.96*S.sd), exp(log(est.vec[,2]) + 1.96*S.sd))

colnames(sensi.mat) = c("zeta.seq", "ML", "ML.sd", "ML.lci", "ML.uci", "S", "S.sd", "S.lci", "S.uci")

sensi.mat 

########################################################
### Figure 3
########################################################
zeta0.seq = seq(0, 0.5, by = 0.01)
zeta1.seq = seq(0, 0.5, by = 0.01)

mor.mat = matrix(NA, nrow = length(zeta0.seq), ncol = length(zeta1.seq))
for(i in 1:length(zeta0.seq)){
  for(j in 1:length(zeta1.seq)){
    mor.mat[i,j] = SP.inference.under(count.mat = count.mat, zeta = c(zeta0.seq[i], zeta1.seq[j]))$est.OR
  }
}

sensi.mor = as.matrix(mor.mat)
contour(zeta0.seq, zeta1.seq, sensi.mor, 
        xlab = expression(zeta[0]), ylab = expression(zeta[1]),
        nlevels = 20, main = "Marginal COR")
abline(a=0, b=1, lty=2)






