### Nonparametric hierchical Bayes model for survival with competing risks
# Currently implemented only for 2 competing risks 

HBM <-function(data,test.data, n.iter, ntree=500,nsplit=0,sigma.error.f=1,sigma.error.g=1,
                         init.lam=NULL,init.logsig=NULL,burn.in=1,thinning=1){
  ## Inputs:
  # - data: dataframe with training data. Indicator of death cause must be labeled "Cause".
  #         censoring is coded with a 0.
  # - test.data: dataframe with test data.
  # - n.iter: number of MCMC iterations to approximate the posterior distributions
  # - ntree: number of trees in multivariate regression model used to estimate beta
  # - nsplit: whether all possible splits should be explored in ecah tree. (0 means yes)
  # - sigma.error.f: variance parameter in prior distribution of f
  # - sigma.error.g: variance parameter in prior distribution of g
  # - burn.in: number of MCMC iterations discarded as burn-in
  # - thinning: include only 1 in every "thinning" iterations of the MCMC sampler to ensure independent samples

  
  ## Outputs:
  # - beta.train: sequence of draws from the posterior distribution of beta parameters for the training data - location parameter
  # - beta.jump.mat: matrix of jumping probabilities to a new proposed beta value in every iteration of sampler
  # - beta.test: sequence of draws from the posterior distribution of beta parameters for the test data
  # - logsig.mat: sequence of draws from the posterior distribution of log(sigma) parameters - scale parameter
  # - logsig.jump.mat: matrix of jumping probabilities to a new proposed log(sigma) value in every iteration of sampler
  # - lam.mat: sequence of draws from the posterior distribution of lambda parameters - shape parameter
  # - lam.jump.mat: matrix of jumping probabilities to a new proposed lambda value in every iteration of sampler
  # - w.mat: sequence of draws from the posterior distribution of w parameters
  # - w.jump.mat: matrix of jumping probabilities to a new proposed w value in every iteration of sampler
  # - w.mat.test: sequence of draws from the posterior distribution of w parameters for the test data
  # - importance.beta: 3 dimensional matrix with the variable importance (on the median time) statistics for each iteration of the MCMC sampler,
  #                   for each feature and for each cause
  # - importance.w: 3 dimensional matrix with the variable importance (on the event probability) statistics for each iteration of the MCMC sampler,
  #                   for each feature and for each cause
  # - Other parameters are internal used for the purpose of debugging
  
  ### import libraries
  library(randomForestSRC)
  library(flexsurv)
  
  ### some redefinitons
  delta                 <- as.numeric(data$Cause == 0)
  d                     <- sum(delta) #number of censored obervations
  n.samples             <- nrow(data)
  n.samples.test        <- nrow(test.data)
  n.var                 <- ncol(data)-2
  n.cause               <- 2
  cause                 <- data$Cause
  cause.mat             <- matrix(c(as.numeric(data$Cause == 1),as.numeric(data$Cause == 2),
                                    as.numeric(data$Cause == 0)),
                                  ncol=n.cause+1)
  time                  <- as.vector(data$Survival) #event time
  
  ### initialize space allocation
  g.of.x.train          <- array(NA,dim=c(n.iter+1,n.samples,n.cause))
  g.of.x.test           <- array(NA,dim=c(n.iter+1,n.samples.test,n.cause))
  f.of.x.train          <- array(NA,dim=c(n.iter+1,n.samples,n.cause))
  f.of.x.test           <- array(NA,dim=c(n.iter+1,n.samples.test,n.cause))
  
  # initial parameter estimates
  initial               <- flexsurvreg(Surv(data$Survival, data$Cause != 0) ~ 1, dist = "gengamma")
  
  beta.mat              <- array(NA,dim=c(n.iter+1,n.samples,n.cause))
  beta.mat[1,,]         <- rnorm(n.cause*n.samples,0,1) + initial$coeff[[1]]
  beta.mat.test         <- array(NA,dim=c(n.iter+1,n.samples.test,n.cause))
  beta.jump.mat         <- matrix(NA,n.iter+1,n.samples) #matrix for storing beta proposals 
  beta.jump.mat.test    <- matrix(NA,n.iter+1,n.samples.test)
  sigma.jump.beta.mat   <- matrix(NA,(n.iter+1),n.samples)
  
  w.mat                 <- array(NA,dim=c(n.iter+1,n.samples,n.cause))
  w.mat[1,,]            <- runif(n.cause*n.samples,4,5)
  w.mat.test            <- array(NA,dim=c(n.iter+1,n.samples.test,n.cause))
  w.jump.mat            <- matrix(NA,n.iter+1,n.samples)  
  w.jump.mat.test       <- matrix(NA,n.iter+1,n.samples.test)
  sigma.jump.w.mat      <- matrix(NA,(n.iter+1),n.samples)
  
  logsig.mat            <- matrix(NA,(n.iter+1),n.cause) #different sigma=scale for causes
  logsig.mat[1,]        <- rnorm(n.cause,0,0.5) + initial$coeff[[2]]
  logsig.jump.mat       <- rep(NA,(n.iter+1)) 
  sigma.jump.logsig.mat <- rep(NA,(n.iter+1)) 
  
  lam.mat               <- matrix(NA,(n.iter+1),n.cause) #different lambda=shape for causes
  lam.mat[1,]           <- rnorm(n.cause,0,1) + initial$coeff[[3]] 
  lam.jump.mat          <- rep(NA,(n.iter+1))
  sigma.jump.lam.mat    <- rep(NA,(n.iter+1)) 
  
  sigma.jump.beta       <- rep(2.5,n.samples)
  sigma.jump.w          <- rep(0.5,n.samples)
  sigma.jump.logsig     <- rep(0.3,n.cause)
  sigma.jump.lam        <- rep(0.5,n.cause)
  
  importance.beta       <- array(NA,dim=c(n.iter,n.var,n.cause))
  importance.w          <- array(NA,dim=c(n.iter,n.var,n.cause))
  
  #log of posterior (likelihood*prior) of w (assignment variable) to be used in MH ratio
  log.post.w<-function(beta,lam,logsig,delta,time,sigma.error.f,f.of.x,cause.mat,w){
    pi <- exp(w)/sum(exp(w))
    sum(cause.mat[1]*(log(0.001+pi[1]*dgengamma(time,mu=beta[1],sigma=exp(logsig[1]),Q=lam[1])))) + 
      sum(cause.mat[2]*(log(0.001+pi[2]*dgengamma(time,mu=beta[2],sigma=exp(logsig[2]),Q=lam[2]))))+
      sum(cause.mat[3]*(log(1.001-pi[1]*pgengamma(time,mu=beta[1],sigma=exp(logsig[1]),Q=lam[1])-
                              pi[2]*pgengamma(time,mu=beta[2],sigma=exp(logsig[2]),Q=lam[2]))))+
      sum(dnorm(w,f.of.x,sigma.error.f,log=TRUE))
    
  }
  
  #MCMC update for w (assignment variable)
  w.update<-function(sigma.jump.w){
    w.star        <- rnorm(n.cause,w.mat[wei,ii,],sigma.jump.w) 
    log.post.old  <- log.post.w(beta.mat[wei,ii,],lam.mat[wei,],logsig.mat[wei,],delta[ii],time[ii],
                             sigma.error.f,f.of.x.train[wei,ii,],cause.mat[ii,], w.mat[wei,ii,])
    log.post.star <- log.post.w(beta.mat[wei,ii,],lam.mat[wei,],logsig.mat[wei,],delta[ii],time[ii],
                              sigma.error.f,f.of.x.train[wei,ii,],cause.mat[ii,], w.star)
    r             <- ifelse(any(w.star<0),-1,exp(log.post.star-log.post.old) )
    w             <- ifelse(c(runif(1)<r,runif(1)<r),w.star,w.mat[wei,ii,])
    
    if(is.na(w)==TRUE||is.infinite(w)==TRUE||is.nan(w)==TRUE) {w<-c(1,1);print('LOG W NA')}
    p.jump        <- min(r,1) # probability of jumping
    
    list(w=w,p.jump=p.jump)
  }
  
  #log of posterior (likelihood*prior) of beta=location to be used in MH ratio
  log.post.beta<-function(beta,lam,logsig,delta,time,sigma.error.g,g.of.x,cause.mat,w){
    pi <- w/sum(w)
    sum(cause.mat[1]*(log(0.001+pi[1]*dgengamma(time,mu=beta[1],sigma=exp(logsig[1]),Q=lam[1])))) + 
      sum(cause.mat[2]*(log(0.001+pi[2]*dgengamma(time,mu=beta[2],sigma=exp(logsig[2]),Q=lam[2]))))+
      sum(cause.mat[3]*(log(1.001-pi[1]*pgengamma(time,mu=beta[1],sigma=exp(logsig[1]),Q=lam[1])-
                              pi[2]*pgengamma(time,mu=beta[2],sigma=exp(logsig[2]),Q=lam[2]))))+
      sum(dnorm(beta,g.of.x,sigma.error.g,log=TRUE))
    
  }
  
  #MCMC update for beta=location
  beta.update<-function(sigma.jump.beta){
    beta.star     <- rnorm(n.cause,beta.mat[wei,ii,],sigma.jump.beta) 
    log.post.old  <- log.post.beta(beta.mat[wei,ii,],lam.mat[wei,],logsig.mat[wei,],delta[ii],time[ii],
                                sigma.error.g,g.of.x.train[wei,ii,],cause.mat[ii,], w.mat[wei+1,ii,])
    log.post.star <- log.post.beta(beta.star,lam.mat[wei,],logsig.mat[wei,],delta[ii],time[ii],
                                 sigma.error.g,g.of.x.train[wei,ii,],cause.mat[ii,], w.mat[wei+1,ii,])
    r             <- exp(log.post.star-log.post.old) 
    beta          <- ifelse(c(runif(1)<r,runif(1)<r),beta.star,beta.mat[wei,ii,])
    
    if(is.na(beta)==TRUE||is.infinite(beta)==TRUE||is.nan(beta)==TRUE) {beta<-c(1,1);print('BETA NA')}
    p.jump        <- min(r,1) # probability of jumping
    
    list(beta=beta,p.jump=p.jump)
  }
  
  
  #log of posterior (likelihood*prior) of sigma to be used in MH ratio
  log.post.sigma<-function(beta,lam,logsig,delta,time,cause.mat,w){
    #summation <- 0
    pi <- w/rowSums(w)
    sum(cause.mat[,1]*(log(0.001+pi[,1]*dgengamma(time,mu=beta[,1],sigma=exp(logsig[1]),Q=lam[1])))) + 
      sum(cause.mat[,2]*(log(0.001+pi[,2]*dgengamma(time,mu=beta[,2],sigma=exp(logsig[2]),Q=lam[2]))))+
      sum(cause.mat[,3]*(log(1.001-pi[,1]*pgengamma(time,mu=beta[,1],sigma=exp(logsig[1]),Q=lam[1])-
                               pi[,2]*pgengamma(time,mu=beta[,2],sigma=exp(logsig[2]),Q=lam[2]))))+
      sum(dnorm(logsig,mean=initial$coeff[[2]],sd=0.5,log=TRUE))
  }
  
  
  
  #MCMC update for log(sigma=scale)
  logsig.update<-function(sigma.jump.logsig){
    logsig.star   <- rnorm(n.cause,logsig.mat[wei,],sigma.jump.logsig)
    log.post.old  <- log.post.sigma(beta.mat[wei+1,,],lam.mat[wei,],logsig.mat[wei,],delta,time,cause.mat,w.mat[wei+1,,])
    log.post.star <- log.post.sigma(beta.mat[wei+1,,],lam.mat[wei,],logsig.star,delta,time,cause.mat,w.mat[wei+1,,])
    r             <- exp(log.post.star-log.post.old)
    logsig        <- ifelse(c(runif(1)<r,runif(1)<r),logsig.star,logsig.mat[wei,])
    
    if(is.na(logsig)==TRUE) {logsig<-logsig.mat[wei,];print('logsig NA')}
    p.jump        <- min(r,1)
    
    list(logsig=logsig,p.jump=p.jump)
  }
  
  #log of posterior (likelihood*prior) of lam (=shape) to be used in MH ratio
  log.post.lambda<-function(beta,lam,logsig,delta,time,cause.mat,w){
    pi <- w/rowSums(w)
    sum(cause.mat[,1]*(log(0.001+pi[,1]*dgengamma(time,mu=beta[,1],sigma=exp(logsig[1]),Q=lam[1])))) + 
      sum(cause.mat[,2]*(log(0.001+pi[,2]*dgengamma(time,mu=beta[,2],sigma=exp(logsig[2]),Q=lam[2]))))+
      sum(cause.mat[,3]*(log(1.001-pi[,1]*pgengamma(time,mu=beta[,1],sigma=exp(logsig[1]),Q=lam[1])-
                               pi[,2]*pgengamma(time,mu=beta[,2],sigma=exp(logsig[2]),Q=lam[2]))))+
      sum(dnorm(lam,mean=initial$coeff[[3]],sd=2,log=TRUE))
  }
  
  #MCMC update for log(lambda=shape)
  lam.update<-function(sigma.jump.lam){
    lam.star      <- rnorm(n.cause,lam.mat[wei,],sigma.jump.lam)
    log.post.old  <- log.post.lambda(beta.mat[wei+1,,],lam.mat[wei,],logsig.mat[wei+1,],delta,time,cause.mat, w.mat[wei+1,,])
    log.post.star <- log.post.lambda(beta.mat[wei+1,,],lam.star,logsig.mat[wei+1,],delta,time,cause.mat, w.mat[wei+1,,])
    r             <- exp(log.post.star-log.post.old)
    lam           <- ifelse(c(runif(1)<r,runif(1)<r),lam.star,lam.mat[wei,])
    if(is.na(lam)==TRUE||is.nan(lam)==TRUE||is.infinite(lam)==TRUE) {lam<-lam.mat[wei,];print('lam NA')}
    p.jump.lam<-min(r,1)
    
    list(lam=lam,p.jump=p.jump.lam)
  }
  
  #FULL MCMC procedure
  for (wei in 1:n.iter){
    data$lab1               <- beta.mat[wei,,1]; data$lab2 <- beta.mat[wei,,2];
    mrf.beta                <- rfsrc(Multivar(lab1,lab2) ~., data = subset(data, select=-c(Survival,Cause)),ntree=ntree, nsplit=nsplit,importance=TRUE) # specify variables since we are adding stuff in data
    g.of.x.train[wei,,1]    <- mrf.beta$regrOutput$lab1$predicted; g.of.x.train[wei,,2]<-mrf.beta$regrOutput$lab2$predicted; 
    predict                 <- predict(mrf.beta, test.data)
    beta.mat.test[wei,,1]   <- predict$regrOutput$lab1$predicted; beta.mat.test[wei,,2] <- predict$regrOutput$lab2$predicted; 
    importance.beta[wei,,1] <- mrf.beta$regrOutput$lab1$importance; importance.beta[wei,,2] <- mrf.beta$regrOutput$lab2$importance;
    print("Done MRF_beta")
    
    data$lab1               <- w.mat[wei,,1]; data$lab2 <- w.mat[wei,,2]; 
    mrf.w                   <- rfsrc(Multivar(lab1,lab2) ~., data = subset(data, select=-c(Survival,Cause)),ntree=ntree, nsplit=nsplit,importance=TRUE) # specify variables since we are adding stuff in data
    f.of.x.train[wei,,1]    <- mrf.w$regrOutput$lab1$predicted; f.of.x.train[wei,,2]<-mrf.w$regrOutput$lab2$predicted; 
    predict                 <- predict(mrf.w, test.data)
    w.mat.test[wei,,1]      <- predict$regrOutput$lab1$predicted; w.mat.test[wei,,2] <- predict$regrOutput$lab2$predicted; 
    importance.w[wei,,1]    <- mrf.w$regrOutput$lab1$importance;importance.w[wei,,2] <- mrf.w$regrOutput$lab2$importance;
    print("Done MRF_w")
    
    for (ii in 1:n.samples){
      if (wei %% 50 == 0) {sigma.jump.w[ii] <- adaptive_jump(sigma.jump.w[ii],w.jump.mat[,ii],wei)}
      temp                          <- w.update(sigma.jump.w[ii])
      w.mat[wei+1,ii,]              <- temp$w
      w.jump.mat[wei+1,ii]          <- temp$p.jump
      sigma.jump.w.mat[wei+1,ii]    <- sigma.jump.w[ii]
    }
    print('DONE w UPDATE')
    
    for (ii in 1:n.samples){
      if (wei %% 50 == 0) {sigma.jump.beta[ii] <- adaptive_jump(sigma.jump.beta[ii],beta.jump.mat[,ii],wei)}
      temp                           <- beta.update(sigma.jump.beta[ii])
      beta.mat[wei+1,ii,]            <- temp$beta
      beta.jump.mat[wei+1,ii]        <- temp$p.jump
      sigma.jump.beta.mat[wei+1,ii]  <- sigma.jump.beta[ii]
    }
    print('DONE BETA UPDATE')
    
    if (wei %% 50 == 0) {sigma.jump.logsig <- adaptive_jump(sigma.jump.logsig,logsig.jump.mat,wei)}
    temp                          <- logsig.update(sigma.jump.logsig)
    logsig.mat[wei+1,]            <- temp$logsig
    logsig.jump.mat[wei+1]        <- temp$p.jump
    sigma.jump.logsig.mat[wei+1]  <- sigma.jump.logsig
    
    print('DONE LOGSIG UPDATE')
    
    # lambda update
    if (wei %% 50 == 0) {sigma.jump.lam <- adaptive_jump(sigma.jump.lam,lam.jump.mat,wei)}
    temp                       <- lam.update(sigma.jump.lam)
    lam.mat[wei+1,]            <- temp$lam
    lam.jump.mat[wei+1]        <- temp$p.jump
    sigma.jump.lam.mat[wei+1]  <- sigma.jump.lam
    
    print('DONE LAMBDA UPDATE')
    print(sprintf("FULL_ITERATION_%i",wei))
  }
  
  
  #post-processing
  draws<-seq(burn.in,n.iter,thinning)
  
  
  list(parameters=list(beta.train=beta.mat[draws,,],beta.jump.mat=beta.jump.mat[draws,],beta.test=beta.mat.test[draws,,],
                       logsig.mat=logsig.mat[draws,],logsig.jump.mat=logsig.jump.mat[draws],lam.mat=lam.mat[draws,],
                       lam.jump.mat=lam.jump.mat[draws],w.mat=w.mat[draws,,],w.jump.mat=w.jump.mat[draws,],
                       w.mat.test=w.mat.test[draws,,],sigma.jump.beta.mat=sigma.jump.beta.mat[draws,],
                       sigma.jump.logsig.mat=sigma.jump.logsig.mat[draws],sigma.jump.lam.mat=sigma.jump.lam.mat[draws],
                       sigma.jump.w.mat=sigma.jump.w.mat[draws,]),
       importance = list(importance.beta=importance.beta[draws,,],importance.w=importance.w[draws,,]))
  
}
