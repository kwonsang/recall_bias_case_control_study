
logit.model.under = function(param, data, zeta = c(0,0)){
  # data
  y = data[,1]
  t = data[,2]
  x = data[,-c(1,2)]
  
  x = as.matrix(x)
  x = cbind(rep(1, length(x[,1])), x)
  
  # parameters
  n = dim(x)[1]
  m = dim(x)[2]
  beta = param[1:m]
  gamma = param[(m+1):(2*m)]
  gamma_t = param[(2*m+1)]
  
  # sensitivity parameter
  zeta0 = zeta[1]
  zeta1 = zeta[2]
  
  t.vec = matrix(1, nrow = dim(x)[1], ncol = 1)
  # outcome models + propensity score model
  m0 = exp((x %*% gamma))/(1 + exp((x %*% gamma)))
  m1 = exp((x %*% gamma) + t.vec*gamma_t)/(1 + exp((x %*% gamma) + t.vec*gamma_t))
  ps = exp((x %*% beta))/(1 + exp((x %*% beta)))
  
  # p_{ab} = pr(Y = a, Tstar = b | X)
  p11 = (1-zeta1)*m1*ps
  p10 = zeta1*m1*ps + m0*(1-ps)
  p01 = (1-zeta0)*(1-m1)*ps
  p00 = zeta0*(1-m1)*ps + (1-m0)*(1-ps)
  
  like.val = log(p11)*y*t + log(p01)*(1-y)*t + log(p10)*y*(1-t) + log(p00)*(1-y)*(1-t)
  sum.like.val = sum(like.val)
  return(-sum.like.val)
}


logit.model.over = function(param, data, eta = c(0,0)){
  # data
  y = data[,1]
  t = data[,2]
  x = data[,-c(1,2)]
  
  x = as.matrix(x)
  x = cbind(rep(1, length(x[,1])), x)
  
  # parameters
  n = dim(x)[1]
  m = dim(x)[2]
  beta = param[1:m]
  gamma = param[(m+1):(2*m)]
  gamma_t = param[(2*m+1)]
  
  # sensitivity parameter
  eta0 = eta[1]
  eta1 = eta[2]
  
  # outcome models + propensity score model
  m0 = exp((x %*% gamma))/(1 + exp((x %*% gamma)))
  m1 = exp((x %*% gamma) + t.vec*gamma_t)/(1 + exp((x %*% gamma) + t.vec*gamma_t))
  ps = exp((x %*% beta))/(1 + exp((x %*% beta)))
  
  # p_{ab} = pr(Y = a, Tstar = b | X)
  p11 = m1*ps + eta1*m0*(1-ps)
  p01 = (1-m1)*ps + eta0*(1-m0)*(1-ps)
  p10 = (1-eta1)*m0*(1-ps)
  p00 = (1-eta0)*(1-m0)*(1-ps)
  
  like.val = log(p11)*y*t + log(p01)*(1-y)*t + log(p10)*y*(1-t) + log(p00)*(1-y)*(1-t)
  sum.like.val = sum(like.val)
  return(-sum.like.val)
}

#################################
## Maximum Likelihood Method
#################################
ML.logistic.under = function(data, zeta = c(0,0)){
  # data
  require(optimx)
  y = data[,1]
  t = data[,2]
  x = data[,-c(1,2)]
  
  x = as.matrix(x)
  init.ps = glm(t ~ x, family="binomial")
  init.m = glm(y ~ x + t, family = "binomial")
  
  x = cbind(rep(1, length(x[,1])), x)
  
  # parameters
  n = dim(x)[1]
  m = dim(x)[2]
  
  init.param = c(init.ps$coefficients, init.m$coefficients)
  
  # opt.res = optim(init.param, logit.model, method = "BFGS", data = data, eta = eta)
  opt.res =suppressWarnings(optimx(init.param, logit.model.under, method = "BFGS", data = data, zeta = zeta))
  
  est.beta = t(opt.res[1, 1:m])
  est.gamma = t(opt.res[1, (m+1):(2*m)])
  est.gamma_t = t(opt.res[1, (2*m+1)])
  
  t.vec = matrix(1, nrow = dim(x)[1], ncol = 1)
  est.e =  exp((x %*% est.beta))/(1 + exp((x %*% est.beta)))
  est.m0 =  exp((x %*% est.gamma))/(1 + exp((x %*% est.gamma)))
  est.m1 =  exp((x %*% est.gamma) + t.vec %*% est.gamma_t)/(1 + exp((x %*% est.gamma) + t.vec %*% est.gamma_t))
  
  est.p0 = mean(est.m0)
  est.p1 = mean(est.m1)
  
  risk.diff = est.p1 - est.p0
  odds.ratio = (est.p1*(1-est.p0))/(est.p0*(1-est.p1))
  
  return(list(potential.outcome = c(est.p1, est.p0), risk.diff = risk.diff, est.OR = odds.ratio, beta = est.beta, gamma = est.gamma, gamma.t = est.gamma_t, ps = est.e, m0 = est.m0, m1 = est.m1, likeli.val = opt.res$value))
}

ML.logistic.over = function(data, eta = c(0,0)){
  # data
  require(optimx)
  y = data[,1]
  t = data[,2]
  x = data[,-c(1,2)]
  
  x = as.matrix(x)
  init.ps = glm(t ~ x, family="binomial")
  init.m = glm(y ~ x + t, family = "binomial")
  
  x = cbind(rep(1, length(x[,1])), x)
  
  # parameters
  n = dim(x)[1]
  m = dim(x)[2]
  
  init.param = c(init.ps$coefficients, init.m$coefficients)
  
  # opt.res = optim(init.param, logit.model, method = "BFGS", data = data, eta = eta)
  opt.res =suppressWarnings(optimx(init.param, logit.model.over, method = "BFGS", data = data, eta = eta))
  
  est.beta = t(opt.res[1, 1:m])
  est.gamma = t(opt.res[1, (m+1):(2*m)])
  est.gamma_t = t(opt.res[1, (2*m+1)])
  
  t.vec = matrix(1, nrow = dim(x)[1], ncol = 1)
  est.e =  exp((x %*% est.beta))/(1 + exp((x %*% est.beta)))
  est.m0 =  exp((x %*% est.gamma))/(1 + exp((x %*% est.gamma)))
  est.m1 =  exp((x %*% est.gamma) + t.vec %*% est.gamma_t)/(1 + exp((x %*% est.gamma) + t.vec %*% est.gamma_t))
  
  est.p0 = mean(est.m0)
  est.p1 = mean(est.m1)
  
  risk.diff = est.p1 - est.p0
  odds.ratio = (est.p1*(1-est.p0))/(est.p0*(1-est.p1))
  
  return(list(potential.outcome = c(est.p1, est.p0), risk.diff = risk.diff, est.OR = odds.ratio, beta = est.beta, gamma = est.gamma, gamma.t = est.gamma_t, ps = est.e, m0 = est.m0, m1 = est.m1, likeli.val = opt.res$value))
}

#################################
## Stratification 
#################################

SP.inference.over = function(count.mat, eta, alpha = 0.05){
  eta0 = eta[1]
  eta1 = eta[2]
  
  ## The following corresponds to the cells in Table 1. 
  a.star = count.mat[,1] # exposed & case
  b.star = count.mat[,2] # exposed & control
  c.star = count.mat[,3] # unexposed & case
  d.star = count.mat[,4] # unexposed & control
  
  ab.star = a.star + b.star; cd.star = c.star + d.star
  n.vec = a.star + b.star + c.star + d.star
  
  a = a.star - (eta1/(1-eta1))*c.star
  b = b.star - (eta0/(1-eta0))*d.star
  c = c.star/(1-eta1)
  d = d.star/(1-eta0)
  
  n1 = a+b
  n0 = c+d
  
  p1 = a/(a+b)
  p0 = c/(c+d)
  
  weighted.p1 = sum(p1*n.vec)/sum(n.vec)
  weighted.p0 = sum(p0*n.vec)/sum(n.vec)
  
  est.mor = (weighted.p1 * (1-weighted.p0))/(weighted.p0 * (1-weighted.p1))
  
  #####
  # variance
  var1 = (1/(weighted.p1*(1-weighted.p1)))^2 * sum((n.vec/sum(n.vec))^2 * (p1*(1-p1)/n1))
  var0 = (1/(weighted.p0*(1-weighted.p0)))^2 * sum((n.vec/sum(n.vec))^2 * (p0*(1-p0)/n0))
  log.var = var0 + var1
  
  #####
  z.score = qnorm(1- alpha/2)
  lower = exp(log(est.mor) - z.score * sqrt(log.var))
  upper = exp(log(est.mor) + z.score * sqrt(log.var))
  conf.int = c(lower, upper)
  
  return(list(est.OR = est.mor, log.var = log.var, conf.int = conf.int))
}

SP.inference.under = function(count.mat, zeta, alpha = 0.05){
  zeta0 = zeta[1]
  zeta1 = zeta[2]
  
  ## The following corresponds to the cells in Table 1. 
  a.star = count.mat[,1] # exposed & case
  b.star = count.mat[,2] # exposed & control
  c.star = count.mat[,3] # unexposed & case
  d.star = count.mat[,4] # unexposed & control
  
  ab.star = a.star + b.star; cd.star = c.star + d.star
  n.vec = a.star + b.star + c.star + d.star
  
  a = a.star/(1-zeta1)
  b = b.star/(1-zeta0)
  c = c.star - (zeta1/(1-zeta1))*a.star
  d = d.star - (zeta0/(1-zeta0))*b.star
  
  n1 = a+b
  n0 = c+d
  
  p1 = a/(a+b)
  p0 = c/(c+d)
  
  weighted.p1 = sum(p1*n.vec)/sum(n.vec)
  weighted.p0 = sum(p0*n.vec)/sum(n.vec)
  
  est.mor = (weighted.p1 * (1-weighted.p0))/(weighted.p0 * (1-weighted.p1))
  
  #####
  # variance
  var1 = (1/(weighted.p1*(1-weighted.p1)))^2 * sum((n.vec/sum(n.vec))^2 * (p1*(1-p1)/n1))
  var0 = (1/(weighted.p0*(1-weighted.p0)))^2 * sum((n.vec/sum(n.vec))^2 * (p0*(1-p0)/n0))
  log.var = var0 + var1
  
  #####
  z.score = qnorm(1- alpha/2)
  lower = exp(log(est.mor) - z.score * sqrt(log.var))
  upper = exp(log(est.mor) + z.score * sqrt(log.var))
  conf.int = c(lower, upper)
  
  return(list(est.OR = est.mor, log.var = log.var, conf.int = conf.int))
}


#################################
## Mantel-Hanszel 
#################################
MH.inference.over = function(count.mat, eta, alpha = 0.05){
  eta0 = eta[1]
  eta1 = eta[2]
  
  ## The following corresponds to the cells in Table 1. 
  a.star = count.mat[,1] # exposed & case
  b.star = count.mat[,2] # exposed & control
  c.star = count.mat[,3] # unexposed & case
  d.star = count.mat[,4] # unexposed & control
  
  ab.star = a.star + b.star; cd.star = c.star + d.star
  n.vec = a.star + b.star + c.star + d.star
  
  a = a.star - (eta1/(1-eta1))*c.star
  b = b.star - (eta0/(1-eta0))*d.star
  c = c.star/(1-eta1)
  d = d.star/(1-eta0)
  
  ## Computation of the RGB variance estimator
  ti = a + b
  ni = a + c
  mi = b + d
  Ni = a + b + c + d
  
  Ri = a*d/Ni
  Si = b*c/Ni
  Nplus = sum(Ni)
  
  Pi = (a + d)/Ni
  Qi = (b + c)/Ni
  Rplus = sum(Ri)
  Splus = sum(Si)
  
  part1 = sum(Pi*Ri)/(2*Rplus^2)
  part2 = sum(Pi*Si + Qi*Ri)/(2*Rplus*Splus)
  part3 = sum(Qi*Si)/(2*Splus^2)
  log.var = part1 + part2 + part3
  
  est.MH = sum(Ri)/sum(Si)
  
  #####
  z.score = qnorm(1- alpha/2)
  lower = exp(log(est.MH) - z.score * sqrt(log.var))
  upper = exp(log(est.MH) + z.score * sqrt(log.var))
  conf.int = c(lower, upper)
  
  return(list(est.OR = est.MH, log.var = log.var, conf.int = conf.int))
}

MH.inference.under = function(count.mat, zeta, alpha = 0.05){
  zeta0 = zeta[1]
  zeta1 = zeta[2]
  
  ## The following corresponds to the cells in Table 1. 
  a.star = count.mat[,1] # exposed & case
  b.star = count.mat[,2] # exposed & control
  c.star = count.mat[,3] # unexposed & case
  d.star = count.mat[,4] # unexposed & control
  
  ab.star = a.star + b.star; cd.star = c.star + d.star
  n.vec = a.star + b.star + c.star + d.star
  
  a = a.star/(1-zeta1)
  b = b.star/(1-zeta0)
  c = c.star - (zeta1/(1-zeta1))*a.star
  d = d.star - (zeta0/(1-zeta0))*b.star
  
  ## Computation of the RGB variance estimator
  ti = a + b
  ni = a + c
  mi = b + d
  Ni = a + b + c + d
  
  Ri = a*d/Ni
  Si = b*c/Ni
  Nplus = sum(Ni)
  
  Pi = (a + d)/Ni
  Qi = (b + c)/Ni
  Rplus = sum(Ri)
  Splus = sum(Si)
  
  part1 = sum(Pi*Ri)/(2*Rplus^2)
  part2 = sum(Pi*Si + Qi*Ri)/(2*Rplus*Splus)
  part3 = sum(Qi*Si)/(2*Splus^2)
  log.var = part1 + part2 + part3
  
  est.MH = sum(Ri)/sum(Si)
  
  #####
  z.score = qnorm(1- alpha/2)
  lower = exp(log(est.MH) - z.score * sqrt(log.var))
  upper = exp(log(est.MH) + z.score * sqrt(log.var))
  conf.int = c(lower, upper)
  
  return(list(est.OR = est.MH, log.var = log.var, conf.int = conf.int))
}

