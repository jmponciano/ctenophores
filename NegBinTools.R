##############  Functions needed to run the models, just run these #########
# Negative Binomial likelihood function (with time)

library(MASS)

negbinom.nll <- function(lguess, elapsed.days,counts.perday){
  
  alpha <- lguess$alpha
  k     <- lguess$k
  
  tvec <- elapsed.days
  x.t  <- counts.perday
  
  n   <- length(counts.perday)
  
  negll.vec <- rep(0,n)
  
  P <- alpha/(alpha+tvec)
  
  for(j in 1:n){
    negll.vec[j] <- -dnbinom(x=x.t[j], size=k, prob=P[j], log=TRUE)
  }
  
  tot.nll <- sum(negll.vec)  
  return(tot.nll)
  
}


negbinom.fit <- function(formula, full.df=mydataf, exp.unit=mydataf$exp.unit,response="embryos"){
  
  # formula = the specific model you are running in standard formula notation
  # full.df = this is the name of the data frame of interest
  # exp.unit = this is the name of the column inside "full.df" which explicitly
  #            keeps track and names each one of the experimental units.
  #            for example...
  
  
  designmat <- model.matrix(formula, data=full.df)
  nbetas    <- ncol(designmat)
  inits <- rep(0.33, nbetas)
  samp.size <- nrow(full.df)
  
  tot.nll <- function(guesses=c(inits,0.33),designmat=designmat,full.df=full.df,response=response){
    
    
    unique.reps <- unique(exp.unit)
    n.exp.units <- length(unique.reps)
    nbetas    <- ncol(designmat)
    log.alpha.inv <- designmat%*%guesses[1:nbetas] 
    full.df$log.alpha.inv <- log.alpha.inv
    k <- exp(guesses[(nbetas+1)])
    
    nll.vec <- rep(0,n.exp.units)
    for(i in 1:n.exp.units){
      
      subdat <- full.df[exp.unit==unique.reps[i],]
      subdat<- subdat[order(subdat$day),]
      sub.alpha <- 1/exp(subdat$log.alpha.inv)
      
      subdays0 <- c(0,subdat$day)
      nsubdays <- length(subdays0)
      elapsed.t<- subdays0[2:nsubdays]-subdays0[1:(nsubdays-1)]
      subdays <- elapsed.t
      if(response=="embryos"){
        subcounts <- subdat$embryos}else if(response=="cydippids")
        {subcounts <- subdat$normal.cyd}
      
      parmlist <- list(alpha=sub.alpha, k=k)      
      nll.vec[i] <- negbinom.nll(lguess=parmlist, elapsed.days=subdays, 
                                 counts.perday=subcounts)
    }
    #print(nll.vec)
    sum.nll <- sum(nll.vec)
    if((sum.nll==0)|is.na(sum.nll)){return(.Machine$double.xmax)}else{
      return(sum.nll)
    }
  }
  
  opt.rslts <- optim(par=c(inits,0.33),fn=tot.nll, method="BFGS", 
                     designmat=designmat,full.df=full.df,response=response,
                     hessian=TRUE)
  
  betas.hat <- opt.rslts$par[1:nbetas]
  k.hat <- exp(opt.rslts$par[(nbetas+1)])
  
  max.llike <- -opt.rslts$value
  
  BIC <- -2*max.llike + (nbetas+1)*log(samp.size)
  
  Fishinv <- ginv(opt.rslts$hessian)
  vars <- diag(Fishinv)
  half.cis <- qnorm(p=0.975)*sqrt(vars)
  
  
  out.list <- list(betas.hat=betas.hat, k.hat=k.hat, max.llike=max.llike, BIC=BIC, my.hess=opt.rslts$hessian, betas.SEs = half.cis)
  return(out.list)
}

# Means and variances
MV.negbinom <- function(lguess, elapsed.days,counts.perday=NULL,new=FALSE){
  
  alpha <- lguess$alpha
  k     <- lguess$k
  tvec  <- elapsed.days
  n     <- length(tvec)
  E.vec <- rep(0,n)
  P <- alpha/(alpha+tvec)
  Q <- 1-P
  E.vec <- (k*Q)/P
  V.vec <- (k*Q)/(P^2)
  
  if(new==FALSE){
    x.t  <- counts.perday
    MVO.mat <- cbind(x.t, E.vec, V.vec)
    colnames(MVO.mat) <- c("Observed", "Expected", "Var.Expected")
    return(MVO.mat)
  }else{
    MV.mat <- cbind(E.vec, V.vec)
    colnames(MV.mat) <- c("Expected", "Var.Expected")
    return(MV.mat)
  }
}

negbinomfit.pred <- function(betas, k, formula, exp.unit,full.df,response="embryos"){
  
  designmat <- model.matrix(formula, data=full.df)
  unique.reps <- unique(exp.unit)
  n.exp.units <- length(unique.reps)
  nbetas    <- ncol(designmat) # has to match 'length(betas)'
  
  log.alpha.inv <- designmat%*%betas 
  full.df$log.alpha.inv <- log.alpha.inv
  
  for(i in 1:n.exp.units){
    
    subdat <- full.df[exp.unit==unique.reps[i],]
    subdat<- subdat[order(subdat$day),]
    sub.alpha <- 1/exp(subdat$log.alpha.inv)
    
    subdays0 <- c(0,subdat$day)
    nsubdays <- length(subdays0)
    elapsed.t<- subdays0[2:nsubdays]-subdays0[1:(nsubdays-1)]
    subdays <- elapsed.t
    subdat$elapsed.days <- elapsed.t
    
    if(response=="embryos"){
      subcounts <- subdat$embryos}else if(response=="cydippids")
      {subcounts <- subdat$normal.cyd}
    
    parmlist <- list(alpha=sub.alpha, k=k)      
    
    MVpreds <- MV.negbinom(lguess=parmlist, elapsed.days=subdays, 
                           counts.perday=subcounts)
    Expected <- MVpreds[,2]
    Var.Expected <- MVpreds[,3]
    Observed <- MVpreds[,1]
    
    subdat$ExpUniNum <- rep(i,length(Observed))
    subdat$tsdays <- subdat$day-1
    subdat$Observed <- Observed
    subdat$Expected <- Expected
    subdat$Var.Expected <- Var.Expected

    if(i==1){
      data.out <- subdat
    }else if(i>1){
      data.out <- rbind(data.out,subdat)
    }
    
  }
  return(data.out) #preds.list
}


negbinom.pred <- function(betas, k, formula, exp.unit,full.df){
  
  designmat <- model.matrix(formula, data=full.df)
  unique.reps <- unique(exp.unit)
  n.exp.units <- length(unique.reps)
  nbetas    <- ncol(designmat) # has to match 'length(betas)'
  
  log.alpha.inv <- designmat%*%betas 
  full.df$log.alpha.inv <- log.alpha.inv
  
  preds.list <- list()
  
  for(i in 1:n.exp.units){
    
    subdat <- full.df[exp.unit==unique.reps[i],]
    subdat<- subdat[order(subdat$day),]
    sub.alpha <- 1/exp(subdat$log.alpha.inv)
    
    subdays0 <- c(0,subdat$day)
    nsubdays <- length(subdays0)
    elapsed.t<- subdays0[2:nsubdays]-subdays0[1:(nsubdays-1)]
    subdays <- elapsed.t
    subdat$elapsed.days <- elapsed.t
    

    parmlist <- list(alpha=sub.alpha, k=k)      
    
    ###### Return here
    MVpreds <- MV.negbinom(lguess=parmlist, elapsed.days=subdays, new=TRUE)
    Expected <- MVpreds[,1]
    Var.Expected <- MVpreds[,2]

    subdat$ExpUniNum <- rep(i,length(Expected))
    subdat$tsdays <- subdat$day-1
    subdat$Expected <- Expected
    subdat$Var.Expected <- Var.Expected

    if(i==1){
      data.out <- subdat
    }else if(i>1){
      data.out <- rbind(data.out,subdat)
    }
    
  }
  return(data.out) 
}


one.sim <- function(lguess,elapsed.days){
  
  # lguess= list of parameter values: [1] predicted alphas using betas hats, 
  # and [2] khat
  
  alpha <- lguess$alpha
  k     <- lguess$k
  tvec <- elapsed.days
  
  n   <- length(tvec)
  
  xt  <- rep(0,n)
  
  P <- alpha/(alpha+tvec)
  
  for(j in 1:n){
    xt[j] <- rnbinom(n=1, size=k, prob=P[j])
  }
  
  nt  <- cumsum(xt)
  outmat <- cbind(nt,xt, elapsed.days)
  colnames(outmat) <- c("Nt", "Xt" , "tau")
  
  return(outmat)  
  
}


one.boot.df <- function(full.df){
  
  orig.expunits <- as.factor(full.df$exp.unit)
  levs4boot <- levels(orig.expunits)
  nlevs.expun <- length(levs4boot)
  npboot.reps <- sample(levs4boot, size=nlevs.expun, replace=TRUE)
  boot.df <- full.df[full.df$exp.unit==npboot.reps[1],]
  beu.count <- rep(1,nrow(boot.df))
  for(i in 2:nlevs.expun){
    ith.boot.dat <- full.df[full.df$exp.unit==npboot.reps[i],]
    boot.df <- rbind(boot.df,ith.boot.dat)
    beu.count <- c(beu.count,rep(i,nrow(ith.boot.dat)))
  }
  beu.countf <- as.factor(beu.count)
  boot.df$beu.countf <- beu.countf
  return(boot.df)
}


negbinom.fit.boot <- function(formula, boot.df=boot.df, response="embryos"){
  
  # formula = the specific model you are running in standard formula notation
  # full.df = this is the name of the data frame of interest
  # exp.unit = this is the name of the column inside "full.df" which explicitly
  #            keeps track and names each one of the experimental units.
  #            for example...
  
  
  designmat <- model.matrix(formula, data=boot.df)
  nbetas    <- ncol(designmat)
  inits <- rep(0.33, nbetas)
  samp.size <- nrow(boot.df)
  
  tot.nll <- function(guesses=c(inits,0.33),designmat=designmat,boot.df=boot.df,response=response){
    
    b.unique.reps <- unique(boot.df$beu.countf)
    bn.exp.units <- length(b.unique.reps)
    nbetas    <- ncol(designmat)
    log.alpha.inv <- designmat%*%guesses[1:nbetas] 
    boot.df$log.alpha.inv <- log.alpha.inv
    k <- exp(guesses[(nbetas+1)])
    
    nll.vec <- rep(0,bn.exp.units)
    for(i in 1:bn.exp.units){
      
      subdat <- boot.df[boot.df$beu.countf==b.unique.reps[i],]
      subdat<- subdat[order(subdat$day),]
      sub.alpha <- 1/exp(subdat$log.alpha.inv)
      
      subdays0 <- c(0,subdat$day)
      nsubdays <- length(subdays0)
      elapsed.t<- subdays0[2:nsubdays]-subdays0[1:(nsubdays-1)]
      subdays <- elapsed.t
      if(response=="embryos"){
        subcounts <- subdat$embryos}else if(response=="cydippids")
        {subcounts <- subdat$normal.cyd}
      
      parmlist <- list(alpha=sub.alpha, k=k)      
      nll.vec[i] <- negbinom.nll(lguess=parmlist, elapsed.days=subdays, 
                                 counts.perday=subcounts)
    }
    #print(nll.vec)
    sum.nll <- sum(nll.vec)
    if((sum.nll==0)|is.na(sum.nll)){return(.Machine$double.xmax)}else{
      return(sum.nll)
    }
  }
  
  opt.rslts <- optim(par=c(inits,0.33),fn=tot.nll, method="BFGS",designmat=designmat,boot.df=boot.df,response=response)
  
  betas.hat <- opt.rslts$par[1:nbetas]
  k.hat <- exp(opt.rslts$par[(nbetas+1)])
  
  max.llike <- -opt.rslts$value
  
  BIC <- -2*max.llike + (nbetas+1)*log(samp.size)
  
  out.list <- list(betas.hat=betas.hat, k.hat=k.hat, max.llike=max.llike, BIC=BIC)
  return(out.list)
}


#####################  End of functions needed to run the code ########``