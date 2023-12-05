########################################################
## Codes for quantile regression for Weibull ###########
## distribution reparameterized by a quantile ##########
## Sanchez et al. (2021) ###############################
########################################################


########################################################
## Required packages
########################################################

require(rootSolve)
require(MASS)
require(maxLik)
require(VGAM) 
require(EnvStats)
require(normalp)
require(cubature)
require(ssym)
library(quantreg)
require(gamlss)
require(gamlss.dist)
require(weibullness)


##############################################################
## PDF
##############################################################

drweibull   <- function(t,Q=0.5,k=1,q=0.5,logar=FALSE)
{
  if(!is.numeric(Q)||!is.numeric(k)||!is.numeric(q))
    if(Q<=0){stop("Q must be positive")}
  if(k<=0) {stop("k must be positive")}
  if(q>=1 || q<=0) {stop("q must be within range (0,1)")}
  
  pdf    <- -k * t^(k-1) * Q^(-k) * log(1-q) * exp( t^k *Q^(-k) * log(1-q)  )
  
  if(logar==TRUE){
    pdf <- log(pdf)
  }
  return(pdf)
}


##############################################################
## CDF
##############################################################


prweibull   <- function(t,Q=0.5,k=1,q=0.5){
  # if(!is.numeric(Q)||!is.numeric(k)||!is.numeric(q))
  #   if(Q<=0){stop("Q must be positive")}
  # if(k<=0) {stop("k must be positive")}
  # if(q>=1 || q<=0) {stop("q must be within range (0,1)")}
  #cdf    <- 1 - (1-q)^((t^k)*Q)
  cdf <- pweibull(t, shape=k, scale = Q*(-log(1-q))^(-1/k), lower.tail = TRUE, log.p = FALSE)
  return(cdf)
}


##############################################################
## Estimation process
##############################################################

quantweibreg.fit = function(y,
                            x,
                            z,
                            quant = 0.5,
                            link = "log"
) 
{
  # x is the design matrix with first column containing only one's
  # z is the design matrix with first column containing only one's
  # y is the vector of values of the response variable  
  
  n           = NROW(x)
  p           = NCOL(x)
  l           = NCOL(z)
  
  linkstr = link
  linkobj = make.link(linkstr)
  linkfun = linkobj$linkfun
  linkinv = linkobj$linkinv
  Q.eta  = linkobj$mu.eta
  
  XY = cbind(y,x[,-1])
  
  if(l>1){
    
  ZY = cbind(y,z[,-1])
  x.u = unique(x[,-1]) # x unico
  z.u = unique(z[,-1]) # z unico
  q.i = c() # vector de cuantiles para cada nivel x_i
  median.i = c() #vector de medianas para cada nivel x.i
  k.i = c()
  lambda.i =c()
  
  for(i in 1:length(x.u)){
    yx.i = which(XY[ ,-1]==x.u[i])
    modelo1 = weibull.mle(y[yx.i], threshold=0)
    lambda.i[i] = modelo1$scale
  }
 
  for(j in 1:length(z.u)){
    yz.i = which(ZY[ ,-1]==z.u[j])
    modelo2 = weibull.mle(y[yz.i], threshold=0)
    k.i[j] = modelo2$shape
  }
  q.i = lambda.i*(-log(1-quant))^(1/k.i)
  regresionQ = lm(linkfun(q.i) ~ x.u)
  regresionk = lm(linkfun(k.i) ~ z.u)
  betai = regresionQ$coefficients
  zetai = regresionk$coefficients
  
  start = list(beta = betai, zeta = zetai)
  if (is.list(start)) 
    start = do.call("c", start)
  
  fr    = function(vp)
  {                                         
    betab = vp[1:p]
    eta1  = as.vector(x%*%betab)
    Q     = linkinv(eta1)
    zetab = vp[-(1:p)]
    eta2  = as.vector(z%*%zetab)
    k     = linkinv(eta2)
    q    = quant 
    vt    = y
    res = sum( log(-log(1-q)) + log(k) + (k-1)*log(vt) - k*log(Q) + 
                 vt^k * Q^(-k) * log(1-q) )
    return(res)
  } 
  
  }else{
    
    ZY = cbind(y,z)
    x.u = unique(x[,-1])
    q.i = c() # vector de cuantiles para cada nivel x_i
    median.i = c() #vector de medianas para cada nivel x.i
    k.i = c()
    lambda.i =c()
    
    for(i in 1:length(x.u)){
      yx.i = which(XY[ ,-1]==x.u[i])
      modelo1 = weibull.mle(y[yx.i], threshold=0)
      lambda.i[i] = modelo1$scale
    }
    
    modelo2 = weibull.mle(y, threshold=0)
    k.i = modelo2$shape
    q.i = lambda.i*(-log(1-quant))^(1/k.i)
    regresionQ = lm(linkfun(q.i) ~ x.u)
    regresionk = lm(linkfun(k.i) ~ 1)
    betai = regresionQ$coefficients
    zetai = regresionk$coefficients
    
    start = list(beta = betai, zeta = zetai)
    if (is.list(start)) 
      start = do.call("c", start)
    
    fr    = function(vp)
    {                                         
      betab = vp[1:p]
      eta1  = as.vector(x%*%betab)
      Q     = linkinv(eta1)
      zetab = vp[-(1:p)]
      eta2  = as.vector(z*zetab)
      k     = linkinv(eta2)
      q    = quant 
      vt    = y
      res = sum( log(-log(1-q)) + log(k) + (k-1)*log(vt) - k*log(Q) + 
                   vt^k * Q^(-k) * log(1-q) )
      return(res)
    } 
  }
  
  
  
  #A = matrix(c(rep(0,p1+p2)),1,p1+p2)
  #B = 0
  
  opt = maxBFGS(fn = fr, start = start)
  
  if (opt$code > 0)  warning("optimization failed to converge")
  
  log.lik.est = opt$maximum
  estimates   = opt$estimate
  
  # Information criteria
  AIC   = - 2 * log.lik.est + 2 * (p+l)
  AICc  = AIC + (2 * (p+l) * ((p+l) + 1)) / (n - (p+l) - 1)
  BIC   = - 2 * log.lik.est + log(n) * (p+l)
  
  beta  = as.vector(estimates[1:p])
  eta1   = as.vector(x%*%beta)
  Q     = linkinv(eta1)
  
  tau   = as.vector(estimates[-(1:p)])
  
  if(l > 1){
    eta2  = as.vector(z%*%tau)
    k     = linkinv(eta2)
  }else{
    eta2  = tau
    k     = linkinv(eta2)
  }
  

  fisher = -opt$hessian
  se     = sqrt(diag(solve(fisher)))
  
  confidence.int <- cbind(estimates - qnorm(0.975)*se, estimates + qnorm(0.975)*se)
  colnames(confidence.int) <- c("2.5%","97.5%")
  
  hess = as.matrix(opt$hessian)
  
  zstatbeta  = beta / se[1:p]
  zstattau   = tau / se[-(1:p)] 
  pvalorbeta = 2 * pnorm(abs(zstatbeta), lower.tail = F)
  pvalortau  = 2 * pnorm(abs(zstattau), lower.tail = F)
  
  names(beta)  = colnames(x)
  names(tau)   = colnames(z)
  
  sebeta = se[1:p]
  setau  = se[-(1:p)]
  tb5 <- miscTools::coefTable(beta,sebeta, df=(n-(p+l)))
  tb6 <- miscTools::coefTable(tau,setau, df=(n-(p+l)))
  
  
  cdfhat      <- prweibull(t=y,Q=Q,k=k,q=quant)
  # 
  GCSresidual <- -log(1-cdfhat)
  RQresidual  <- qnorm(cdfhat)
  #hist(RQresidual)
  # 
  layout(matrix(c(1:2), nrow=1, byrow=TRUE))
  layout.show(2)

  envelopeCS(GCSresidual)
  #par(new=T)
  envelopeRQ(RQresidual)
  
  #RMSE <- sqrt(mean((y-Q)^2))
  
  
  # cat("\n")
  # cat("--------------------------------------------------------------\n")
  # cat("             Quantile Weibull Regression Model                \n")
  # cat("--------------------------------------------------------------\n")
  # cat("--------------------------------------------------------------\n")
  # cat("Maximum Likelihood estimation \n")
  # cat("Log-Likelihood:", log.lik.est, "\n")
  # cat("AIC:", AIC, "AICc:", AICc, "BIC:", BIC, "\n")
  # cat("Number of observations:", n, "\n")
  # cat("Quantile of interest:", quant, "\n")
  # cat("--------------------------------------------------------------\n")
  # cat("Quantile - Coefficients:\n")
  # printCoefmat(tb5, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
  # cat("--------------------------------------------------------------\n")
  # cat("k - Coefficients :\n")
  # printCoefmat(tb6, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
  # cat("--------------------------------------------------------------\n")
  # 
  
  
  rval = list(coefficients = list(beta = beta, tau = tau), 
              se = se, 
              conf.int = confidence.int,
              pvalor = list(beta = pvalorbeta, tau = pvalortau),
              #y = y,
              #x = x,
              Initial.values = start, 
              converged = opt$message, 
              information.criterions = list(aic = AIC,bic = BIC,aicc=AICc), 
              loglik = log.lik.est,
              n = n, p = p, l=l, quantile = quant,
              #Hessian = hess,
              #RMSE     = RMSE#,
              GCSresidual = as.vector(GCSresidual),
              RQresidual = as.vector(RQresidual)#,
              #fitted.values.Q = structure(Q, .Names = names(y)),
              #fitted.values.k = structure(k, .Names = names(y))
              )
  
  
  return(rval)
  
  
}


require(qqplotr)


# ## Otra forma de hacer los gráficos con envelope
# resid = data.frame(RQ = RQresidual)
# ggplot(data = resid, mapping = aes(sample = RQ)) +
#   stat_qq_band() +
#   stat_qq_line() +
#   stat_qq_point() +
#   labs(x = "theoretical quantiles", y = "sample quantiles")
# 
# 
# resid = data.frame(GCS = GCSresidual)
# ggplot(data = resid, mapping = aes(sample = GCS)) +
#   stat_qq_band(distribution="exp", dp =1) +
#   stat_qq_line(distribution="exp", dp =1) +
#   stat_qq_point(distribution="exp", dp =1) +
#   labs(x = "theoretical quantiles", y = "sample quantiles")



envelopeCS <- function(x){
  U	         <- x
  n	         <- length(x)
  d2s 	     <- sort(U)
  xq2 	     <- qexp(ppoints(n))
  Xsim 	     <- matrix(0, 100, n)
  for(i in 1:100){
    u2       <- rexp(n)
    Xsim[i,] <- u2
  }
  Xsim2      <- apply(Xsim, 1, sort)
  d21        <- matrix(0, n, 1)
  d22        <- matrix(0, n, 1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,], 0.025)
    d22[i]  <- quantile(Xsim2[i,], 0.975)
  }
  d2med      <- apply(Xsim2, 1, mean)
  fy         <- range(d2s, d21, d22)
  plot(xq2, d2s, xlab = quote("theoretical quantile"),
       ylab = quote("empirical quantile"), 
       pch = 20, ylim = fy)
  par(new = T)
  plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "")
  par(new = T)
  plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "")
  par(new = T)
  plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "")
}



envelopeRQ <- function(x){
  U	         <- x
  n	         <- length(x)
  d2s 	     <- sort(U)
  xq2 	     <- qnorm(ppoints(n))
  Xsim 	     <- matrix(0, 100, n)
  for(i in 1:100){
    u2       <- rnorm(n)
    Xsim[i,] <- u2
  }
  Xsim2      <- apply(Xsim, 1, sort)
  d21        <- matrix(0, n, 1)
  d22        <- matrix(0, n, 1)
  for(i in 1:n){
    d21[i]  <- quantile(Xsim2[i,], 0.025)
    d22[i]  <- quantile(Xsim2[i,], 0.975)
  }
  d2med      <- apply(Xsim2, 1, mean)
  fy         <- range(d2s, d21, d22)
  plot(xq2, d2s, xlab = quote("theoretical quantile"),
       ylab = quote("empirical quantile"), 
       pch = 20, ylim = fy)
  par(new = T)
  plot(xq2, d21, type = "l", ylim = fy, xlab = "", ylab = "")
  par(new = T)
  plot(xq2, d2med, type = "l", ylim = fy, xlab = "", ylab = "")
  par(new = T)
  plot(xq2, d22, type = "l", ylim = fy, xlab = "", ylab = "")
}


##########################################################
## Simulated data 1 ######################################
##########################################################

x <- runif(5)
z <- runif(5)
#X      <- model.matrix( ~ x)
X      <- cbind(1,x)
Z      <- cbind(1,z)
beta   <- c(0.5,1)
zeta   <- c(0.2, 1)
q      <- 0.5 ## quantile
Q      <- exp(as.vector(X%*%beta))
k      <- exp(as.vector(Z%*%zeta))
lambda <- Q * (-log(1-q))^(-1/k)

Y = c()
X = c()
Z = c()
for(i in 1:5){
  Y      <- c(Y,rweibull(500, shape = k[i], scale = lambda[i]))
  X      <- c(X, rep(x[i], 500))
  Z      <- c(Z, rep(z[i], 500))
}

X = cbind(1,X)
Z = cbind(1,Z)


quant=0.5
link="log"
x=X
z=Z
y=Y

quantweibreg.fit(x=X, z=Z, y=Y, quant=0.5, link="log")


# ##########################################################
# ## Simulated data 2 ######################################
# ##########################################################
# 
# x <- runif(200000)
# #X      <- model.matrix( ~ x)
# X      <- cbind(1,x)
# Z      <- cbind(1,x)
# beta   <- c(0.5,1)
# zeta   <- c(0.2, 1)
# q      <- 0.5 ## quantile
# Q      <- exp(as.vector(X%*%beta))
# k      <- exp(as.vector(Z%*%zeta))
# lambda <- Q * (-log(1-q))^(-1/k)
# Y      <- rweibull(200000, shape = k, scale = lambda)
# quant=0.5
# link="log"
# x=X
# z=Z
# y=Y
# 
# quantweibreg.fit(x=X, z=Z, y=Y, quant=0.5, link="log")