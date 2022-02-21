# Set working directory to Folder containing the "Data" and the "Rcode" folders:
#data.dir <- 
#Rcode.dir <-

##############  Functions needed to run the models, run these first #########
# Negative Binomial likelihood function (with time)

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
  
  opt.rslts <- optim(par=c(inits,0.33),fn=tot.nll, method="BFGS",designmat=designmat,full.df=full.df,response=response)
  
  betas.hat <- opt.rslts$par[1:nbetas]
  k.hat <- exp(opt.rslts$par[(nbetas+1)])
  
  max.llike <- -opt.rslts$value
  
  BIC <- -2*max.llike + (nbetas+1)*log(samp.size)
  
  out.list <- list(betas.hat=betas.hat, k.hat=k.hat, max.llike=max.llike, BIC=BIC)
  return(out.list)
}

# Means and variances
MV.negbinom <- function(lguess, elapsed.days,counts.perday){
  
  alpha <- lguess$alpha
  k     <- lguess$k
  
  tvec <- elapsed.days
  x.t  <- counts.perday
  
  n   <- length(counts.perday)
  
  E.vec <- rep(0,n)
  
  P <- alpha/(alpha+tvec)
  Q <- 1-P
  
  E.vec <- (k*Q)/P
  V.vec <- (k*Q)/(P^2)
  
  MVO.mat <- cbind(x.t, E.vec, V.vec)
  colnames(MVO.mat) <- c("Observed", "Expected", "Var.Expected")
  return(MVO.mat)
  
}

negbinomfit.pred <- function(betas, k, formula, exp.unit,full.df,response="embryos"){
  
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
    
    
    preds.list[[i]] <- MVpreds
    
    if(i==1){
      data.out <- subdat
    }else if(i>1){
      data.out <- rbind(data.out,subdat)
    }
    
  }
  return(data.out) #preds.list
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


#####################  End of functions needed to run the code ########


########## FIGURE 2 ##########

## effect of size on reproductive output ##

dissogeny_size <- read.csv("~/Ml_spawning_cydippids/dissogeny_size.csv")

#Fig 2 A - histogram of observations of single animals (mixed sizes)
#N = 56 parental cydippids from 4 bio reps checked daily
library(ggplot2)
ggplot(dissogeny_size, aes(embryos/number)) + geom_histogram(binwidth = 1, fill = 'skyblue', color = 'black', bins = 20) + 
  xlab("number of embryos produced") + ylab("number of observed occurrences") +
  stat_bin(binwidth = 1, geom="text", size = 2, aes(label=..count..), vjust=-1.5) + ylim(0, 225) +
  theme(text = element_text(size = 16))

#drop days with 0 production as in P. globosa paper to see that clutch size/body size relationship is similar
#Fig 2 B - clutch size vs body size
#N = 32 parentals from 4 bio reps, producing 98 clutches
clutch_no_0s <- read.csv("~/Ml_spawning_cydippids/clutch_no_0s.csv")
ggplot(clutch_no_0s, aes(size, (embryos/number))) + geom_boxplot(aes(group=size)) + geom_smooth(method = lm) +
  xlab("parental size (mm)") + ylab("clutch size") + 
  theme(text = element_text(size = 16))

# now to deploy our custom negative binomial, time-explicit model:

size.data <- read.csv(paste0(data.dir,"dissogeny_size_singletons.csv"))
names(size.data)

null.model <- negbinom.fit(formula=~1, full.df=size.data, exp.unit = size.data$exp.unit) # 918.2136
model1 <- negbinom.fit(formula=~size, full.df=size.data, exp.unit = size.data$exp.unit) # 912.2451 not a great predictor
model2 <- negbinom.fit(formula=~number, full.df=size.data, exp.unit = size.data$exp.unit) # 923.7704 all number = 1
model3 <- negbinom.fit(formula=~bio.rep, full.df=size.data, exp.unit = size.data$exp.unit) # 899.5539
model4 <- negbinom.fit(formula=~bio.rep+day, full.df=size.data, exp.unit = size.data$exp.unit) # 901.5416
model5 <- negbinom.fit(formula=~day, full.df=size.data, exp.unit = size.data$exp.unit) # 894.5517 *
model6 <- negbinom.fit(formula=~bio.rep+size, full.df=size.data, exp.unit = size.data$exp.unit) # 901.8521
model7 <- negbinom.fit(formula=~day+number, full.df=size.data, exp.unit = size.data$exp.unit) # 900.1085

boot.dat <- one.boot.df(full.df=size.data)
boot.model <- negbinom.fit.boot(formula=~day, boot.df=boot.dat)

#Bootstrap run of models 3 to 7:
B <- 100
bic.bootmat <- matrix(0,nrow=B,ncol=5)
cis.3 <- matrix(0,nrow=B,ncol=4)
cis.4 <- matrix(0,nrow=B,ncol=5)
cis.5 <- matrix(0,nrow=B,ncol=2)
cis.6 <- matrix(0,nrow=B,ncol=5)
cis.7 <- matrix(0,nrow=B,ncol=3)
dat.list <- list()
for(i in 1:B){
  
  boot.dat <- one.boot.df(full.df=size.data)
  
  b3 <- negbinom.fit.boot(formula=~bio.rep, boot.df=boot.dat)
  b4 <- negbinom.fit.boot(formula=~bio.rep+day, boot.df=boot.dat)
  b5 <- negbinom.fit.boot(formula=~day, boot.df=boot.dat)
  b6 <- negbinom.fit.boot(formula=~bio.rep+size, boot.df=boot.dat)
  b7 <- negbinom.fit.boot(formula=~day+number, boot.df=boot.dat)  

  bic.bootmat[i,] <- c(b3$BIC, b4$BIC, b5$BIC, b6$BIC, b7$BIC)

  cis.3[i,] <- b3$betas.hat
  cis.4[i,] <- b4$betas.hat  
  cis.5[i,] <- b5$betas.hat      
  cis.6[i,] <- b6$betas.hat
  cis.7[i,] <- b7$betas.hat  
  
  dat.list[[i]] <- boot.dat
}

fitted.mods <- c("mod3", "mod4", "mod5", "mod6", "mod7")
best.index <- apply(bic.bootmat,1,FUN=function(x){which(x==min(x), arr.ind=TRUE)})
boot.support <- rep("NA",B)
for(i in 1:B){
  boot.support[i] <- fitted.mods[best.index[i]]
}
supports.run1 <- table(boot.support)/B
print(supports.run1)

delta.bic.mat <- matrix(0,nrow=B,ncol=5)
delta.bic.rank <- matrix(0, nrow=B, ncol=5)
delta.bic1 <- rep(0,B)
delta.bic2 <- rep(0,B)
for(i in 1:B){
  x <- bic.bootmat[i,]
  best <- which(x==min(x), arr.ind=TRUE)
  ith.delta.bics <- x - x[best]
  delta.bic.rank[i,] <- rank(ith.delta.bics)
  delta.bic.mat[i,] <- ith.delta.bics
  first.delta <- ith.delta.bics[-best]
  second.best <- which.min(first.delta)
  
  delta.bic1[i] <- min(first.delta)
  delta.bic2[i] <- min(first.delta[-second.best])
}

first.second.deltas <- cbind(delta.bic1,delta.bic2)

mod5.isbest <- which(boot.support=="mod5",arr.ind=TRUE)
mod5.support <- first.second.deltas[mod5.isbest,]
apply(mod5.support,2,mean)

all.where5notbest <- delta.bic.rank[delta.bic.rank[,3]!=1,]

# When Model 5 was not best, which rank did it have?
table(all.where5notbest[,3])

#print(supports.run1)
#boot.support
#mod3 mod4 mod6 mod7 
#0.03 0.18 0.27 0.52   

#####
# Do goodness of fit tests for both model5 (~day), which is the best BIC (894), as well as model6 (~bio.rep+size) (901) ?
# (i.e. because size is shown as graph do we want to evaluate it? it's also only marginally worse than day for cydippids (below))
#####

# Example of how to run the predictions:  the negbinomfit.pred outputs the original data set, re-shuffled to
# match the experimental units order, along with the Expected (predicted) mean counts per experimental unit
# and their associated variance.  The output is then a new data.frame from which variables can be extracted
# to do plots. 

#grouping by size

#Fig 2 C
#model to see output grouped by size
modelX <- negbinom.fit(formula=~size, full.df=size.data, exp.unit = size.data$exp.unit) # 912.2451 *

trial.pred <- negbinomfit.pred(betas=modelX$betas.hat, k=modelX$k.hat, formula=~size, 
                               exp.unit=size.data$exp.unit, full.df=size.data, 
                               response="embryos")

# Trying a time series(staircase) plot for embryos ~size:

max.days <- max(trial.pred$tsdays)

data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

#countsbysizeNrep <- table(trial.pred$size,trial.pred$ExpUniNum)

mycols0 <- c("#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a")

size1 <- data4plot$size[1]
col1 <- mycols0[size1]
par(oma=c(1,1,1,1), mar=c(4,4,3,2), mfrow=c(1,1))
plot(days4plot, cum.expected, type="s", lwd=2, col=col1, bty="l", ylim=c(0,150), xlim=c(0,max.days+1),
     xlab="Days", ylab="Predicted mean cumulative embryos", cex.lab=1.25)


nexpus <- max(trial.pred$ExpUniNum)
for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  
  ith.size <- data4plot$size[1]
  ith.col  <- mycols0[ith.size]
  points(days4plot, cum.expected, type="s", lwd=2, col=ith.col, bty="l",lty=1)
  
}
legend("topleft", legend=1:8, col=mycols0, lty=rep(1,length(mycols0), lwd=2), 
       lwd=rep(1,length(mycols0)), cex=1.1, bty="n")


# ### recolor to color only the 2-4 mm size range - not shown in paper ###
# #model to see output grouped by size
# modelX <- negbinom.fit(formula=~size, full.df=size.data, exp.unit = size.data$exp.unit) # 912.2451 *
# 
# trial.pred <- negbinomfit.pred(betas=modelX$betas.hat, k=modelX$k.hat, formula=~size, 
#                                exp.unit=size.data$exp.unit, full.df=size.data, 
#                                response="embryos")
# 
# # time series(staircase) plot for embryos ~size:
# 
# max.days <- max(trial.pred$tsdays)
# 
# data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
# cum.expected <- cumsum(data4plot$Expected)
# days4plot <- data4plot$tsdays
# 
# #countsbysizeNrep <- table(trial.pred$size,trial.pred$ExpUniNum)
# 
# mycols0 <- c("#999999", "#fcc5c0", "#fa9fb5", "#f768a1", "#999999", "#999999", "#999999", "#999999")
# 
# size1 <- data4plot$size[1]
# col1 <- mycols0[size1]
# par(oma=c(1,1,1,1), mar=c(4,4,3,2))
# plot(days4plot, cum.expected, type="s", lwd=2, col=col1, bty="l", ylim=c(0,150), xlim=c(0,max.days+1),
#      xlab="Days", ylab="Cumulative number of embryos", cex.lab=1.25)
# 
# 
# nexpus <- max(trial.pred$ExpUniNum)
# for(i in 2:nexpus){
#   
#   data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
#   cum.expected <- cumsum(data4plot$Expected)
#   days4plot <- data4plot$tsdays
#   
#   ith.size <- data4plot$size[1]
#   ith.col  <- mycols0[ith.size]
#   points(days4plot, cum.expected, type="s", lwd=2, col=ith.col, bty="l",lty=1)
#   
# }
# legend("topleft", legend=1:8, col=mycols0, lty=rep(1,length(mycols0), lwd=2), 
#        lwd=rep(1,length(mycols0)), cex=1.1, bty="n")

#Fig 2 D
# Trying a time series(staircase) plot for normal 24h cydippids ~size:
where.nas2 <- which(is.na(size.data$normal.cyd), arr.ind=TRUE)
size.data.nonas <- size.data[-where.nas2,]

null.model <- negbinom.fit(formula=~1, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 467.377
model1 <- negbinom.fit(formula=~size, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 454.0638
model2 <- negbinom.fit(formula=~number, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 472.8107
model3 <- negbinom.fit(formula=~bio.rep, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 452.0115
model4 <- negbinom.fit(formula=~bio.rep+day, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 451.5707
model5 <- negbinom.fit(formula=~day, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 444.4695 *
model6 <- negbinom.fit(formula=~tech.rep, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 625.2783

model7 <- negbinom.fit(formula=~day+number, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 449.9032
model8 <- negbinom.fit(formula=~bio.rep+size, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 447.3678

#####

# this breaks down by size and is shown
modelX <- negbinom.fit(formula=~size, full.df=size.data.nonas, exp.unit = size.data.nonas$exp.unit, response="cydippids") # 454.0638 *

trial.pred <- negbinomfit.pred(betas=modelX$betas.hat, k=modelX$k.hat, formula=~size, 
                               exp.unit=size.data.nonas$exp.unit, full.df=size.data.nonas, 
                               response="cydippids")

# Trying a time series (staircase) plot:

max.days <- max(trial.pred$tsdays)

data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

mycols0 <- c("#fde0dd", "#fcc5c0", "#fa9fb5", "#f768a1", "#dd3497", "#ae017e", "#7a0177", "#49006a")

size1 <- data4plot$size[1]
col1 <- mycols0[size1]
par(oma=c(1,1,1,1), mar=c(4,4,3,2))
plot(days4plot, cum.expected, type="s", lwd=1.5, col=col1, bty="l", ylim=c(0,150), xlim=c(0,max.days+1),
     xlab="Days", ylab="Predicted mean cumulative 24h cydippids", cex.lab=1.25)


nexpus <- max(trial.pred$ExpUniNum)
for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  
  ith.size <- data4plot$size[1]
  ith.col  <- mycols0[ith.size]
  points(days4plot, cum.expected, type="s", lwd=1.5, col=ith.col, bty="l",lty=1)
  
}
legend("topleft", legend=1:8, col=mycols0, lty=rep(1,length(mycols0)), 
       lwd=rep(1,length(mycols0)), cex=1.1, bty="n")


########## FIGURE 3 ##########

## effect of temperature on reproductive output ###########################

#temperature and embryo output
temperature.fname <- paste0(data.dir,"temp-fecundity-revised.csv") 
temp.dataf <-read.csv(temperature.fname)

names(temp.dataf)

#t.test((embryos/number) ~ temperature, data = temp.dataf)
#t.test((normal.cyd/number) ~ temperature, data = temp.dataf)

# full.df=mydataf, exp.unit=mydataf$exp.unit
null.model <- negbinom.fit(formula=~1, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit)  # 552.014
model1 <- negbinom.fit(formula=~day, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit) # 546.3797
model2 <- negbinom.fit(formula=~temperature, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit) # 538.3136
model3 <- negbinom.fit(formula=~number, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit) # 550.5957
model4 <- negbinom.fit(formula=~bio.rep, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit) # 556.1047
model5 <- negbinom.fit(formula=~bio.rep+temperature, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit) # 537.3869
# this aligns with my expectations based on 1st pass analysis that temperature is positively correlated with fecundity
# now to check including accounting for variable numbers of parentals

model1.1 <- negbinom.fit(formula=~number+day, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit) # 544.7987
model1.2 <- negbinom.fit(formula=~number+temperature, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit) # 532.1247 *
model1.4 <- negbinom.fit(formula=~number+bio.rep, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit) # 556.2008

#####

#Figure 3 A - embryo output by temperature
model1.2 <- negbinom.fit(formula=~number+temperature, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit) # 532.1247 *

trial.pred <- negbinomfit.pred(betas=model1.2$betas.hat, k=model1.2$k.hat, formula=~number+temperature, 
                               exp.unit=temp.dataf$exp.unit, full.df=temp.dataf, 
                               response="embryos")

max.days <- max(trial.pred$tsdays)
pretty.bioreps <- as.factor(trial.pred$bio.rep)
levels(pretty.bioreps) <- c("R1", "R2", "R3", "R4")
trial.pred$pretty.bioreps <- pretty.bioreps

data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

par(oma=c(1,1,1,1), mar=c(4,4,3,2))
plot(days4plot, cum.expected, type="s", lwd=1.5, col="#56B4E9", bty="l", ylim=c(0,660), xlim=c(0,max.days+1),
     xlab="Days", ylab="Predicted mean cumulative embryos", cex.lab=1.25)
text(x=max(days4plot)+0.35, y=max(cum.expected), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)

nexpus <- max(trial.pred$ExpUniNum)
for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  
  ith.temp <- data4plot$temperature[1]
  if(ith.temp==20){ith.col<-"#56B4E9";my.lty =1}else if(ith.temp==28){ith.col<-"#E69F00";my.lty =1}
  points(days4plot, cum.expected, type="s", lwd=1.5, col=ith.col, bty="l",lty=my.lty)
  text(x=max(days4plot)+0.35, y=max(cum.expected), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)
  
}
legend("topleft", legend=c("28 C.", "20 C."), col=c("#E69F00", "#56B4E9"), lty=c(1,1), lwd=c(1.5,1.5), cex=1.1, bty="n")

# cydippid output by temperature

null.model <- negbinom.fit(formula=~1, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit,response="cydippids") # 488.5174
model1 <- negbinom.fit(formula=~day, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit,response="cydippids") # 482.9868
model2 <- negbinom.fit(formula=~temperature, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit,response="cydippids") # 468.4962 
model3 <- negbinom.fit(formula=~number, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit,response="cydippids") # 489.5525
model4 <- negbinom.fit(formula=~bio.rep, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit, response="cydippids") # 494.4993
model5 <- negbinom.fit(formula=~bio.rep+temperature, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit, response="cydippids") # 468.101
model6 <- negbinom.fit(formula=~bio.rep+number, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit, response="cydippids") # 495.9482
model7 <- negbinom.fit(formula=~day+number, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit, response="cydippids") # 483.7533

model2.1 <- negbinom.fit(formula=~number+temperature, full.df=temp.dataf, exp.unit=temp.dataf$exp.unit,response="cydippids") # 464.9463 *

#####

# Figure 3 B:
trial.pred <- negbinomfit.pred(betas=model2.1$betas.hat, k=model2.1$k.hat, formula=~number+temperature, 
                               exp.unit=temp.dataf$exp.unit, full.df=temp.dataf, 
                               response="cydippids")

max.days <- max(trial.pred$tsdays)
pretty.bioreps <- as.factor(trial.pred$bio.rep)
levels(pretty.bioreps) <- c("R1", "R2", "R3", "R4")
trial.pred$pretty.bioreps <- pretty.bioreps

data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

par(oma=c(1,1,1,1), mar=c(4,4,3,2))
plot(days4plot, cum.expected, type="s", lwd=1.5, col="#56B4E9", bty="l", ylim=c(0,500), xlim=c(0,max.days+1),
     xlab="Days", ylab="Predicted mean cumulative 24h cydippids", cex.lab=1.25)
text(x=max(days4plot)+0.25, y=max(cum.expected), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)

nexpus <- max(trial.pred$ExpUniNum)
for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  
  ith.temp <- data4plot$temperature[1]
  if(ith.temp==20){ith.col<-"#56B4E9";my.lty =1}else if(ith.temp==28){ith.col<-"#E69F00";my.lty =1}
  points(days4plot, cum.expected, type="s", lwd=1.5, col=ith.col, bty="l",lty=my.lty)
  text(x=max(days4plot)+0.25, y=max(cum.expected), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)
  
}
legend("topleft", legend=c("28 C.", "20 C."), col=c("#E69F00", "#56B4E9"), lty=c(1,1), lwd=c(1.5,1.5), cex=1.1, bty="n")

# Figure 3 D:
#copied here from MeanEffectComparisonCalculator.R so all figs can be produced here
source("NegBinTools.R")
library(gam)

# Model 1 setting:  ~number
mod1.test <- negbinom.fit(formula=~number, full.df=temp.dataf, 
                          exp.unit=temp.dataf$exp.unit)

# Prediction for model 1:
mod1.pred <-negbinomfit.pred(betas=mod1.test$betas.hat, k=mod1.test$k.hat, formula=~number, 
                             exp.unit=temp.dataf$exp.unit, full.df=temp.dataf)

# Model 2 setting: ~number + temperature
mod2.test <- negbinom.fit(formula=~number+temperature, full.df=temp.dataf, 
                          exp.unit=temp.dataf$exp.unit) 

# Prediction for model 2:
mod2.pred <- negbinomfit.pred(betas=mod2.test$betas.hat, k=mod2.test$k.hat, 
                              formula=~number+temperature, exp.unit=temp.dataf$exp.unit, 
                              full.df=temp.dataf)


temp.fecund.preds <- data.frame(number =mod1.pred[,4], model1.pred= mod1.pred[,13], 
                                model2.pred=mod2.pred[,13], bio.rep = mod1.pred$bio.rep,
                                temp= mod2.pred[,2])

temp.f20 <- temp.fecund.preds[temp.fecund.preds$temp==20,]
temp.f28 <- temp.fecund.preds[temp.fecund.preds$temp==28,]

lm20 <- lm(model2.pred~number+ I(number^2), data=temp.f20)
lm28 <- lm(model2.pred~number+ I(number^2), data=temp.f28)

df.4pred <- data.frame(number=7:20)
num.vec <- df.4pred$number
pred.lm20 <- predict(lm20,newdata=df.4pred, se.fit=TRUE)$fit/df.4pred$number
pred.lm28 <- predict(lm28,newdata=df.4pred, se.fit=TRUE)$fit/df.4pred$number

beta.hats <- c(-2.16306112,  0.08242668,  0.17022280) # these numbers are from the object 'mod2.test'
beta.ses <- c(1.79162612, 0.04780375, 0.06267349) # these numbers are from the object 'mod2.test'

propag.relses <- 1.5*mean(beta.ses/abs(beta.hats)) # Propagating the uncertainty from the standard error of the
# from the mles into the predicted trend.

#set color transparency
transp.col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}
## END
two.cols <- c("#56B4E9", "#E69F00") 
two.cols.light <- c(transp.col(color="#9ad2f2",percent=50),transp.col(color="#ffc034",percent=50))

#setting the alpha level
alpha <- 0.05
directses.betas <- qnorm(1-(alpha/2))*sqrt(diag(solve(mod2.test$my.hess)))[1:3]
propag.relses <- mean(directses.betas/abs(beta.hats))

# below produces an error: Error in rgb(255, 192, 52) : color intensity 255, not in [0,1]
# two.cols.light <- c(rgb(154, 210, 242, max=255,alpha=),rgb(255, 192, 52))

par(oma=c(1,1,1,1), mar=c(4,4,1,1))
plot(num.vec, pred.lm20,type="l", lwd=2, 
     ylab="Per capita fecundity", xlab="Number of parentals", xlim=c(6.5,21), ylim=c(0,7),
     bty="l")
polygon(x=c(num.vec,rev(num.vec)), y=c(pred.lm20-propag.relses*pred.lm20,
                                       rev(pred.lm20+propag.relses*pred.lm20)), col=two.cols.light[1], border=NA)
points(num.vec,pred.lm20, type="l", lwd=2, col=two.cols[1])

points(num.vec, pred.lm28,type="l", lwd=2, 
       ylab="Per capita fecundity", xlab="number of parentals", xlim=c(5,21), ylim=c(3,100),
       bty="l")
polygon(x=c(num.vec,rev(num.vec)), y=c(pred.lm28-propag.relses*pred.lm28,
                                       rev(pred.lm28+propag.relses*pred.lm28)), col=two.cols.light[2], border=NA)
points(num.vec,pred.lm28, type="l", lwd=2, col=two.cols[2])
legend("topleft", legend=c("20 degrees", "28 degrees"), col=two.cols, lty=1,lwd=2, bty="n")

points(num.vec, pred.lm20,type="l", lwd=2, 
       ylab="Per capita fecundity", xlab="Number of parentals", xlim=c(6.5,21), ylim=c(0,7),
       bty="l")
polygon(x=c(num.vec,rev(num.vec)), y=c(pred.lm20-propag.relses*pred.lm20,
                                       rev(pred.lm20+propag.relses*pred.lm20)), col=two.cols.light[1], border=NA)
points(num.vec,pred.lm20, type="l", lwd=2, col=two.cols[1])

# note: the transparent color cannot be saved as an EPS

########## FIGURE 4 ##########

## effect of culture density on reproductive output ##

#density 1 - dedicated experiment: 50 ml cultures with 1, 2, 4, or 8 individuals ######################################################
#embryos ~density
dens1.data <- read.csv(paste0(data.dir,"dissogeny_density_1.csv"))
names(dens1.data)

null.model <- negbinom.fit(formula=~1, full.df=dens1.data, exp.unit = dens1.data$exp.unit) # 987.667
model1 <- negbinom.fit(formula=~number, full.df=dens1.data, exp.unit = dens1.data$exp.unit) # 993.1294
model2 <- negbinom.fit(formula=~day, full.df=dens1.data, exp.unit = dens1.data$exp.unit) # 993.046
model3 <- negbinom.fit(formula=~bio.rep, full.df=dens1.data, exp.unit = dens1.data$exp.unit) # 905.7654 best BIC for single parameter models
#in contrast to temp, effect of bio rep bigger than any experimental factor; 
#     this aligns with my expectations based on 1st pass analysis that density is not a major factor/trigger

model4 <- negbinom.fit(formula=~bio.rep+day, full.df=dens1.data, exp.unit = dens1.data$exp.unit) # 906.7068
model5 <- negbinom.fit(formula=~number+day, full.df=dens1.data, exp.unit = dens1.data$exp.unit) # 998.5531

model1.2 <- negbinom.fit(formula=~day+number, full.df=dens1.data, exp.unit = dens1.data$exp.unit) # 998.5531
model1.3 <- negbinom.fit(formula=~bio.rep+number, full.df=dens1.data, exp.unit = dens1.data$exp.unit) # 903.8298 * slightly lower than ~bio.rep


#####
# Do goodness of fit tests for model1.3 (bio.rep+number) - best BIC (903)
#####

## Fig 4 A - embryos ~density; linear scale hard to read so log-transformed
trial.pred <- negbinomfit.pred(betas=model1.3$betas.hat, k=model1.3$k.hat, formula=~bio.rep+number, 
                               exp.unit=dens1.data$exp.unit, full.df=dens1.data, 
                               response="embryos")

max.days <- max(trial.pred$tsdays)
pretty.bioreps <- as.factor(trial.pred$bio.rep)
levels(pretty.bioreps) <- c("R1", "R2", "R3", "R4")
trial.pred$pretty.bioreps <- pretty.bioreps

data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

par(oma=c(1,1,1,1), mar=c(4,4,3,2))
plot(days4plot, log(cum.expected), type="s", lwd=1.5, col = "#a1dab4", bty="l", ylim=log(c(0.1,100)), xlim=c(1,max.days+1),
     xlab="Days", ylab="Predicted mean cumulative embryos (log)", cex.lab=1.25)

text(x=max(days4plot)+0.25, y=log(max(cum.expected)), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)

nexpus <- max(trial.pred$ExpUniNum)
for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  
  ith.dens <- data4plot$number[1]
  if(ith.dens==1){ith.col<-"#a1dab4";my.lty=1}else if(ith.dens==2){ith.col<-"#41b6c4";my.lty=5}else if(ith.dens==4){ith.col<-"#2c7fb8";my.lty=2}else if(ith.dens==8){ith.col<-"#253494";my.lty=3}
  points(days4plot, log(cum.expected), type="s", lwd=1.5, col=ith.col, bty="l",lty=my.lty)
  text(x=max(days4plot)+0.25, y=log(max(cum.expected)), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)
  
}
legend("bottomright", legend=c("1 per 50 ml", "2 per 50 ml", "4 per 50 ml", "8 per 50 ml"), col=c("#a1dab4", "#41b6c4", "#2c7fb8", "#253494"), lty=c(1,5,2,3), lwd=c(1.5,1.5), cex=1, bty="n")

#cydippids ~density
dens1.data <- read.csv(paste0(data.dir,"dissogeny_density_1.csv")) 
names(dens1.data)

where.nas <- which(is.na(dens1.data$normal.cyd), arr.ind=TRUE)
dens1.data.nonas <- dens1.data[-where.nas,]

null.model <- negbinom.fit(formula=~1, full.df=dens1.data.nonas, exp.unit = dens1.data.nonas$exp.unit, response="cydippids") # 388.9906
model1 <- negbinom.fit(formula=~number, full.df=dens1.data.nonas, exp.unit = dens1.data.nonas$exp.unit, response="cydippids") # 393.4938
model2 <- negbinom.fit(formula=~day, full.df=dens1.data.nonas, exp.unit = dens1.data.nonas$exp.unit, response="cydippids") # 387.9596
model3 <- negbinom.fit(formula=~bio.rep, full.df=dens1.data.nonas, exp.unit = dens1.data.nonas$exp.unit, response="cydippids") # 360.8277 *
model4 <- negbinom.fit(formula=~bio.rep+day, full.df=dens1.data.nonas, exp.unit = dens1.data.nonas$exp.unit, response="cydippids") # 364.9057
model5 <- negbinom.fit(formula=~bio.rep+number, full.df=dens1.data.nonas, exp.unit = dens1.data.nonas$exp.unit, response="cydippids") # 352.5247
model6 <- negbinom.fit(formula=~day+number, full.df=dens1.data.nonas, exp.unit = dens1.data.nonas$exp.unit, response="cydippids") # 391.9752

#####


# Fig 4B - density (1-8) per 50 ml effect on normal offspring (log scale)

model3.1 <- negbinom.fit(formula=~bio.rep+number, full.df=dens1.data.nonas, exp.unit = dens1.data.nonas$exp.unit, response="cydippids") # 352.5247 *

trial.pred <- negbinomfit.pred(betas=model3.1$betas.hat, k=model3.1$k.hat, formula=~bio.rep+number, 
                               exp.unit=dens1.data.nonas$exp.unit, full.df=dens1.data.nonas, 
                               response="cydippids")

max.days <- max(trial.pred$tsdays)
pretty.bioreps <- as.factor(trial.pred$bio.rep)
levels(pretty.bioreps) <- c("R1", "R2", "R3", "R4")
trial.pred$pretty.bioreps <- pretty.bioreps

data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

par(oma=c(1,1,1,1), mar=c(4,4,3,2))
plot(days4plot, log(cum.expected), type="s", lwd=1.5, col = "#a1dab4", bty="l", ylim=log(c(0.1,100)), xlim=c(1,max.days+1),
     xlab="Days", ylab="Predicted mean cumulative 24 h cydippids (log)", cex.lab=1.25)

text(x=max(days4plot)+0.25, y=log(max(cum.expected)), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)

nexpus <- max(trial.pred$ExpUniNum)
for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  
  ith.dens <- data4plot$number[1]
  if(ith.dens==1){ith.col<-"#a1dab4";my.lty=1}else if(ith.dens==2){ith.col<-"#41b6c4";my.lty=5}else if(ith.dens==4){ith.col<-"#2c7fb8";my.lty=2}else if(ith.dens==8){ith.col<-"#253494";my.lty=3}
  points(days4plot, log(cum.expected), type="s", lwd=1.5, col=ith.col, bty="l",lty=my.lty)
  text(x=max(days4plot)+0.25, y=log(max(cum.expected)), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)
  
}
legend("topleft", legend=c("1 per 50 ml", "2 per 50 ml", "4 per 50 ml", "8 per 50 ml"), col=c("#a1dab4", "#41b6c4", "#2c7fb8", "#253494"), lty=c(1,5,2,3), lwd=c(1.5,1.5), cex=1, bty="n")


# #what about total number of animals in a culture vessel (not adjusted for total volume)?:
# 
# dens2.data <- read.csv(paste0(data.dir,"dissogeny_density_2.csv")) #########################################################################
# names(dens2.data)
# 
# null.model <- negbinom.fit(formula=~1, full.df=dens2.data, exp.unit = dens2.data$exp.unit) # 1775.546
# model1 <- negbinom.fit(formula=~number, full.df=dens2.data, exp.unit = dens2.data$exp.unit) # 1658.436
# model2 <- negbinom.fit(formula=~number+day, full.df=dens2.data, exp.unit = dens2.data$exp.unit) # 1659.695 
# model3 <- negbinom.fit(formula=~bio.rep, full.df=dens2.data, exp.unit = dens2.data$exp.unit) # 1665.123
# model4 <- negbinom.fit(formula=~day, full.df=dens2.data, exp.unit = dens2.data$exp.unit) # 1758.209
# model5 <- negbinom.fit(formula=~bio.rep+number, full.df=dens2.data, exp.unit = dens2.data$exp.unit) # 1550.587 *
# model6 <- negbinom.fit(formula=~bio.rep+day, full.df=dens2.data, exp.unit = dens2.data$exp.unit) # 1670.72

#number alone and bio.rep alone have similar BIC but biorep+number is much better than any other model
# data labels make it pretty clear that bio.rep dominates
# not broken down by number at larger numbers

#####

# # effect of total number of parentals on embryos
# model5 <- negbinom.fit(formula=~bio.rep+number, full.df=dens2.data, exp.unit = dens2.data$exp.unit) # 1550.587 *
# 
# trial.pred <- negbinomfit.pred(betas=model5$betas.hat, k=model5$k.hat, formula=~bio.rep+number, 
#                                exp.unit=dens2.data$exp.unit, full.df=dens2.data, 
#                                response="embryos")
# 
# max.days <- max(trial.pred$tsdays)
# pretty.bioreps <- as.factor(trial.pred$bio.rep)
# levels(pretty.bioreps) <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")
# trial.pred$pretty.bioreps <- pretty.bioreps
# 
# data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
# cum.expected <- cumsum(data4plot$Expected)
# days4plot <- data4plot$tsdays
# 
# par(oma=c(1,1,1,1), mar=c(4,4,3,2))
# plot(days4plot, log(cum.expected), type="s", lwd=1.5, col = "#fd8d3c", bty="l", ylim=log(c(0.1,1000)), xlim=c(1,max.days+1),
#      xlab="Days", ylab="Predicted mean cumulative embryos (log)", cex.lab=1.25)
# 
# text(x=max(days4plot)+0.25, y=log(max(cum.expected)), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)
# 
# nexpus <- max(trial.pred$ExpUniNum)
# for(i in 2:nexpus){
#   
#   data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
#   cum.expected <- cumsum(data4plot$Expected)
#   days4plot <- data4plot$tsdays
#   
#   ith.dens <- data4plot$number[1]
#   if(ith.dens==1){ith.col<-"#fd8d3c";my.lty =1}else if(ith.dens==2){ith.col<-"#fc4e2a";my.lty =1}else if(ith.dens==4){ith.col<-"#e31a1c";my.lty =1}else if(ith.dens==8){ith.col<-"#bd0026";my.lty =1}else if(ith.dens>8){ith.col<-"#800026";my.lty =1}
#   points(days4plot, log(cum.expected), type="s", lwd=1.5, col=ith.col, bty="l",lty=my.lty)
#   text(x=max(days4plot)+0.25, y=log(max(cum.expected)), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)
#   
# }
# legend("bottomright", legend=c("1", "2", "4", "8", ">8"), col=c("#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"), lty=c(1,1,1,1,1), lwd=c(1.5,1.5), cex=1.1, bty="n")


# #  effect of number of total parentals on cydippids
# dens2.data <- read.csv(paste0(data.dir,"dissogeny_density_2.csv"))  
# names(dens2.data)
# 
# where.nas2 <- which(is.na(dens2.data$normal.cyd), arr.ind=TRUE)
# dens2.data.nonas <- dens2.data[-where.nas2,]
# 
# null.model <- negbinom.fit(formula=~1, full.df=dens2.data.nonas, exp.unit = dens2.data.nonas$exp.unit, response="cydippids") # 954.8747
# model1 <- negbinom.fit(formula=~number, full.df=dens2.data.nonas, exp.unit = dens2.data.nonas$exp.unit, response="cydippids") # 884.3817
# model2 <- negbinom.fit(formula=~bio.rep, full.df=dens2.data.nonas, exp.unit = dens2.data.nonas$exp.unit, response="cydippids") # 840.4439
# model3 <- negbinom.fit(formula=~day, full.df=dens2.data.nonas, exp.unit = dens2.data.nonas$exp.unit, response="cydippids") # 947.8156
# model4 <- negbinom.fit(formula=~bio.rep+day, full.df=dens2.data.nonas, exp.unit = dens2.data.nonas$exp.unit, response="cydippids") # 842.1732
# model5 <- negbinom.fit(formula=~number+day, full.df=dens2.data.nonas, exp.unit = dens2.data.nonas$exp.unit, response="cydippids") # 868.6837
# model6 <- negbinom.fit(formula=~bio.rep+number, full.df=dens2.data.nonas, exp.unit = dens2.data.nonas$exp.unit, response="cydippids") # 804.0545
# 
# #####

# # Fig - pred mean cum cyd ~density (all culture sizes)
# model2.1 <- negbinom.fit(formula=~bio.rep+number, full.df=dens2.data.nonas, exp.unit = dens2.data.nonas$exp.unit, response="cydippids")  # 804.0545 *
# 
# trial.pred <- negbinomfit.pred(betas=model2.1$betas.hat, k=model2.1$k.hat, formula=~bio.rep+number, 
#                                exp.unit=dens2.data.nonas$exp.unit, full.df=dens2.data.nonas, 
#                                response="cydippids")
# 
# max.days <- max(trial.pred$tsdays)
# pretty.bioreps <- as.factor(trial.pred$bio.rep)
# levels(pretty.bioreps) <- c("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8")
# trial.pred$pretty.bioreps <- pretty.bioreps
# 
# data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
# cum.expected <- cumsum(data4plot$Expected)
# days4plot <- data4plot$tsdays
# 
# # staircase plot:
# par(oma=c(1,1,1,1), mar=c(4,4,3,2))
# plot(days4plot, log(cum.expected), type="s", lwd=1.5, col = "#fd8d3c", bty="l", ylim=log(c(0.1,1000)), xlim=c(1,max.days+1),
#      xlab="Days", ylab="Predicted mean cumulative 24h cydippids (log)", cex.lab=1.25)
# 
# text(x=max(days4plot)+0.25, y=log(max(cum.expected)), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)
# 
# nexpus <- max(trial.pred$ExpUniNum)
# for(i in 2:nexpus){
#   
#   data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
#   cum.expected <- cumsum(data4plot$Expected)
#   days4plot <- data4plot$tsdays
#   
#   ith.dens <- data4plot$number[1]
#   if(ith.dens==1){ith.col<-"#fd8d3c";my.lty =1}else if(ith.dens==2){ith.col<-"#fc4e2a";my.lty =1}else if(ith.dens==4){ith.col<-"#e31a1c";my.lty =1}else if(ith.dens==8){ith.col<-"#bd0026";my.lty =1}else if(ith.dens>8){ith.col<-"#800026";my.lty =1}
#   points(days4plot, log(cum.expected), type="s", lwd=1.5, col=ith.col, bty="l",lty=my.lty)
#   text(x=max(days4plot)+0.25, y=log(max(cum.expected)), labels=as.character(data4plot$pretty.bioreps[1]), cex=0.7)
#   
# }
# legend("bottomright", legend=c("1", "2", "4", "8", ">8"), col=c("#fd8d3c", "#fc4e2a", "#e31a1c", "#bd0026", "#800026"), lty=c(1,1,1,1,1), lwd=c(1.5,1.5), cex=1.1, bty="n")

# Fig 4 D - effect of culture density on fecundity, separated by biological replicate
source("NegBinTools.R")
library(gam)

dens1.data <- read.csv(paste0(data.dir,"dissogeny_density_1.csv"))
names(dens1.data)

mod1.test <-  negbinom.fit(formula=~bio.rep, full.df=dens1.data, 
                           exp.unit = dens1.data$exp.unit)

mod1.pred <- negbinomfit.pred(betas=mod1.test$betas.hat, k=mod1.test$k.hat, 
                              formula=~bio.rep, 
                              exp.unit=dens1.data$exp.unit, full.df=dens1.data, 
                              response="embryos")

mod2.test <- negbinom.fit(formula=~bio.rep+number, full.df=dens1.data, 
                          exp.unit = dens1.data$exp.unit)

mod2.pred <- negbinomfit.pred(betas=mod2.test$betas.hat, k=mod2.test$k.hat, 
                              formula=~bio.rep+number, 
                              exp.unit=dens1.data$exp.unit, full.df=dens1.data, 
                              response="embryos")

##### Figure comparing the smoothed predictions 

dens1.2modpreds <- data.frame(number=mod1.pred$number, model1.pred=mod1.pred[,13], model2.pred=mod2.pred[,13], 
                              bio.rep=mod1.pred$bio.rep)

betas.mod2 <- c(mod2.test$betas.hat, mod2.test$k.hat)
betasses.mod2 <- mod2.test$betas.SEs

prop.error <- mean(betasses.mod2/betas.mod2)

br.1df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R1",]
br.2df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R2",]
br.3df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R3",]
br.4df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R4",]

lm1 <- lm(model2.pred~number+I(number^2), data=br.1df)
lm2 <- lm(model2.pred~number+I(number^2), data=br.2df)
lm3 <- lm(model2.pred~number+I(number^2), data=br.3df)
lm4 <- lm(model2.pred~number+I(number^2), data=br.4df)

df.4pred <- data.frame(number=1:8)
num.vec <- df.4pred$number
pred1 <- predict(lm1, newdata=df.4pred)/df.4pred$number
pred2 <- predict(lm2, newdata=df.4pred)/df.4pred$number
pred3 <- predict(lm3, newdata=df.4pred)/df.4pred$number
pred4 <- predict(lm4, newdata=df.4pred)/df.4pred$number

# Sensible cols
my.cols  <- c("#a1dab4", "#41b6c4", "#2c7fb8", "#253494")
light.cols <- c("#c6e8d1", "#7bccd6", "#59a4d7","#384ccd")

par(oma=c(1,1,1,1), mar=c(4,4,1,1))
# plot for R1
plot(num.vec, pred1, type="l", lwd=2, ylab="Fecundity per capita", xlab="Number of parentals",
     xlim=c(0,9), ylim=c(0,8), bty="l")
polygon(x=c(num.vec,rev(num.vec)), y=c(pred1-prop.error*pred1, rev(pred1+prop.error*pred1)), 
        col=light.cols[1], border=NA)
points(num.vec, pred1, type="l", lwd=2, col=my.cols[1])

# plot for R2
points(num.vec, pred2, type="l", lwd=2)
polygon(x=c(num.vec,rev(num.vec)), y=c(pred2-prop.error*pred2, rev(pred2+prop.error*pred2)), 
        col=light.cols[2], border=NA)
points(num.vec, pred2, type="l", lwd=2, col=my.cols[2])

# plot for R3
points(num.vec, pred3, type="l", lwd=2)
polygon(x=c(num.vec,rev(num.vec)), y=c(pred3-prop.error*pred3, rev(pred3+prop.error*pred3)), 
        col=light.cols[3], border=NA)
points(num.vec, pred3, type="l", lwd=2, col=my.cols[3])

# plot for R4
points(num.vec, pred4, type="l", lwd=2)
polygon(x=c(num.vec,rev(num.vec)), y=c(pred4-prop.error*pred4, rev(pred4+prop.error*pred4)), 
        col=light.cols[4], border=NA)
points(num.vec, pred4, type="l", lwd=2, col=my.cols[4])

legend("topleft", legend=c("R1", "R2", "R3", "R4"), col=my.cols, lty=1, lwd=2, bty="n")


########## FIGURE 5 ##########

## effect of diet on reproductive output ##

#embryos ~diet
feed2.data <- read.csv(paste0(data.dir,"feed2.csv")) #################################################################################
names(feed2.data)

null.model <- negbinom.fit(formula=~1, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 502.2316
model1 <- negbinom.fit(formula=~food, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 501.7295 
model2 <- negbinom.fit(formula=~number, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 502.7569
model3 <- negbinom.fit(formula=~bio.rep, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 497.2842 *
model4 <- negbinom.fit(formula=~day, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 506.3784
model5 <- negbinom.fit(formula=~parental.dpf, full.df=feed2.data, exp.unit = feed2.data$exp.unit) #  505.047
model6 <- negbinom.fit(formula=~DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 496.4327 *
model7 <- negbinom.fit(formula=~ARA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 505.7573
model8 <- negbinom.fit(formula=~EPA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 506.3734
model9 <- negbinom.fit(formula=~protein, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 505.308
model10 <- negbinom.fit(formula=~lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 499.2041 *
model11 <- negbinom.fit(formula=~carbohydrate, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 502.7014

model12 <- negbinom.fit(formula=~bio.rep+number+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 484.7735
model13 <- negbinom.fit(formula=~bio.rep+number+lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 491.4732
model14 <- negbinom.fit(formula=~bio.rep+number, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 501.4043
model15 <- negbinom.fit(formula=~bio.rep+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 480.6468
model16 <- negbinom.fit(formula=~bio.rep+lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 487.3263
model16 <- negbinom.fit(formula=~number+lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 496.8494
model17 <- negbinom.fit(formula=~number+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 493.2619
model18 <- negbinom.fit(formula=~bio.rep+ARA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 499.2634
model19 <- negbinom.fit(formula=~bio.rep+EPA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 501.1771
model20 <- negbinom.fit(formula=~lipid+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 500.5327
model21 <- negbinom.fit(formula=~lipid+DHA+bio.rep, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 484.784

#####

model6.1 <- negbinom.fit(formula=~number+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 493.2619
model10.1 <- negbinom.fit(formula=~lipid+number, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 496.8494

#best model:
model15 <- negbinom.fit(formula=~bio.rep+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 480.6468

#########################
# Fig 5A - DHA staircase plot (log):
trial.pred <- negbinomfit.pred(betas=model15$betas.hat, k=model15$k.hat, formula=~bio.rep+DHA, 
                               exp.unit=feed2.data$exp.unit, full.df=feed2.data, 
                               response="embryos")

max.days <- max(trial.pred$tsdays)
data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

mycols0 <- c("#785EF0", "#DC267F","#648FFF","#FE6100")
mycols <- rep(mycols0, each=2)

mylty <- c(lty=1,lty=1,lty=1,lty=2)
myltys <- rep(mylty, each=2)

par(oma=c(1,1,1,1), mar=c(4,4,3,2))
plot(days4plot, log(cum.expected), type="s", lwd=2, col=mycols[1], bty="l", log(c(0.1,1000)), xlim=c(0,max.days+1), 
     xlab="Days", ylab="Predicted mean cumulative embryos (log)", cex.lab=1.25)

nexpus <- max(trial.pred$ExpUniNum)

for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  points(days4plot, log(cum.expected), type="s", lwd=2, col=mycols[i], bty="l",lty=myltys[i])
  
}
dha.vals <- as.character(unique(trial.pred$DHA))
my.legend <- paste0("DHA mg/g = ",dha.vals[c(3,1,2,4)])
legend("bottomright", legend=my.legend, col=mycols0[c(3,1,2,4)], lty=c(1,1,1,2), 
       lwd=c(1.5,1.5,1.5,1.5), cex=1.1, bty="n")


####################
# Fig 5C - total % lipid staircase plot for embryos (log):
model16 <- negbinom.fit(formula=~bio.rep+lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit) # 487.3263

trial.pred <- negbinomfit.pred(betas=model16$betas.hat, k=model16$k.hat, formula=~bio.rep+lipid, 
                               exp.unit=feed2.data$exp.unit, full.df=feed2.data, 
                               response="embryos")

max.days <- max(trial.pred$tsdays)
data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

mycols0 <- c("#785EF0", "#DC267F","#648FFF","#FE6100")
mycols <- rep(mycols0, each=2)

mylty <- c(lty=1,lty=1,lty=1,lty=1)
myltys <- rep(mylty, each=2)

par(oma=c(1,1,1,1), mar=c(4,4,3,2))
plot(days4plot, log(cum.expected), type="s", lwd=2, col=mycols[1], bty="l", ylim=log(c(0.1,1000)), xlim=c(0,max.days+1),
     xlab="Days", ylab="Predicted mean cumulative embryos (log)", cex.lab=1.25)

nexpus <- max(trial.pred$ExpUniNum)

for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  points(days4plot, log(cum.expected), type="s", lwd=2, col=mycols[i], bty="l",lty=myltys[i])
  
}
lipid.vals <- as.character(unique(trial.pred$lipid))
my.legend <- paste0("% lipid = ",lipid.vals[c(3,1,2,4)])
legend("bottomright", legend=my.legend, col=mycols0[c(3,1,2,4)], lty=c(1,1,1,1), 
       lwd=c(1.5,1.5,1.5,1.5), cex=1.1, bty="n")

#cydippids ~diet
feed2.data <- read.csv(paste0(data.dir,"feed2.csv"))
names(feed2.data)

null.model <- negbinom.fit(formula=~1, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 352.6024
feed2.data$food.type <- as.factor(feed2.data$food.type)
model1 <- negbinom.fit(formula=~food.type, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 336.6597 *
#feed type is correlated with reproductive output
model2 <- negbinom.fit(formula=~number, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 354.7691

feed2.data$bio.rep <- as.factor(feed2.data$bio.rep)
model3 <- negbinom.fit(formula=~bio.rep, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 351.4429
model4 <- negbinom.fit(formula=~day, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") #  356.7601
model5 <- negbinom.fit(formula=~parental.dpf, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 355.5525
model6 <- negbinom.fit(formula=~DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 337.9548 
model7 <- negbinom.fit(formula=~ARA, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 356.644
model8 <- negbinom.fit(formula=~EPA, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 356.4496
model9 <- negbinom.fit(formula=~protein, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 356.3046
model10 <- negbinom.fit(formula=~lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 332.6167 
model11 <- negbinom.fit(formula=~carbohydrate, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 353.329

model12 <- negbinom.fit(formula=~bio.rep+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 326.6361 *

model13 <- negbinom.fit(formula=~bio.rep+number+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 330.7605
model14 <- negbinom.fit(formula=~bio.rep+number+lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 329.5366
model15 <- negbinom.fit(formula=~bio.rep+lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 325.7526 *
model16 <- negbinom.fit(formula=~bio.rep+number+lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 329.5366
model17 <- negbinom.fit(formula=~number+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 339.5698
model18 <- negbinom.fit(formula=~number+lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 335.286

#####
# Fig 5B - cydippids ~DHA staircase plot:
feed2.data <- read.csv(paste0(data.dir,"feed2.csv"))
names(feed2.data)
feed2.data$bio.rep <- as.factor(feed2.data$bio.rep)

model12 <- negbinom.fit(formula=~bio.rep+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 326.6361

trial.pred <- negbinomfit.pred(betas=model12$betas.hat, k=model12$k.hat, formula=~bio.rep+DHA, 
                               exp.unit=feed2.data$exp.unit, full.df=feed2.data, 
                               response="cydippids")

max.days <- max(trial.pred$tsdays)
data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

mycols0 <- c("#785EF0", "#DC267F","#648FFF","#FE6100")
mycols <- rep(mycols0, each=2)

mylty <- c(lty=1,lty=1,lty=1,lty=1)
myltys <- rep(mylty, each=2)

par(oma=c(1,1,1,1), mar=c(4,4,3,2))
plot(days4plot, log(cum.expected), type="s", lwd=2, col=mycols[1], bty="l", log(c(0.1,1000)), xlim=c(0,max.days+1), 
     xlab="Days", ylab="Predicted mean cumulative cydippids (log)", cex.lab=1.25)

nexpus <- max(trial.pred$ExpUniNum)

for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  points(days4plot, log(cum.expected), type="s", lwd=2, col=mycols[i], bty="l",lty=myltys[i])
  
}
dha.vals <- as.character(unique(trial.pred$DHA))
my.legend <- paste0("DHA mg/g = ",dha.vals[c(3,1,2,4)])
legend("bottomright", legend=my.legend, col=mycols0[c(3,1,2,4)], lty=c(1,1,1,1), 
       lwd=c(1.5,1.5,1.5,1.5), cex=1.1, bty="n")


# Fig 5D - cydippids ~lipid staircase plot:

model15 <- negbinom.fit(formula=~bio.rep+lipid, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 325.7526 

trial.pred <- negbinomfit.pred(betas=model15$betas.hat, k=model15$k.hat, formula=~bio.rep+lipid, 
                               exp.unit=feed2.data$exp.unit, full.df=feed2.data, 
                               response="cydippids")

max.days <- max(trial.pred$tsdays)
data4plot <- trial.pred[trial.pred$ExpUniNum==1,]
cum.expected <- cumsum(data4plot$Expected)
days4plot <- data4plot$tsdays

mycols0 <- c("#785EF0", "#DC267F","#648FFF","#FE6100")
mycols <- rep(mycols0, each=2)

mylty <- c(lty=1,lty=1,lty=1,lty=1)
myltys <- rep(mylty, each=2)

par(oma=c(1,1,1,1), mar=c(4,4,3,2))
plot(days4plot, log(cum.expected), type="s", lwd=2, col=mycols[1], bty="l", ylim=log(c(0.01,300)), xlim=c(0,max.days+1), 
    xlab="Days", ylab="Predicted mean cumulative cydippids (log)", cex.lab=1.25)

nexpus <- max(trial.pred$ExpUniNum)

for(i in 2:nexpus){
  
  data4plot <- trial.pred[trial.pred$ExpUniNum==i,]
  cum.expected <- cumsum(data4plot$Expected)
  days4plot <- data4plot$tsdays
  points(days4plot, log(cum.expected), type="s", lwd=2, col=mycols[i], bty="l",lty=myltys[i])
  
}
lipid.vals <- as.character(unique(trial.pred$lipid))
my.legend <- paste0("% lipid = ",lipid.vals[c(3,1,2,4)])
legend("bottomright", legend=my.legend, col=mycols0[c(3,1,2,4)], lty=c(1,1,1,1), 
       lwd=c(1.5,1.5,1.5,1.5), cex=1.1, bty="n")



############## END OF MAIN PAPER FIGURES ##############


########## SUPPLEMENTAL FIGURES ##########

############  Example of one simulated trajectory for Figure 1 in the appendix

###########   Grab the betas hats for model12, get the first experimental unit conditions
###########   predict the alpha for that first experimental unit and use these alpha values to 
###########   simulate data for the desired vector of observed elapsed days (5)

############  Plot real data of embryos from single individuals that grow well past reported pause size
###########   Overlay information about measured size on each day 

source("NegBinTools.R")

revisions.data <- data.frame(read.csv(paste0(data.dir,"single-Ml-spawn-growth_v3.csv")))

#cydippids ~diet
feed2.data <- read.csv(paste0(data.dir,"feed2.csv"))
model12 <- negbinom.fit(formula=~bio.rep+DHA, full.df=feed2.data, exp.unit = feed2.data$exp.unit, response="cydippids") # 326.6361

trial.pred <- negbinomfit.pred(betas=model12$betas.hat, k=model12$k.hat, formula=~bio.rep+DHA, 
                               exp.unit=feed2.data$exp.unit, full.df=feed2.data, 
                               response="cydippids")
designmat <- model.matrix(~bio.rep+DHA, data=feed2.data)
unique.reps <- unique(feed2.data$exp.unit)
n.exp.units <- length(unique.reps)
nbetas    <- ncol(designmat) 

log.alpha.inv <- designmat%*%model12$betas.hat 
feed2.data$log.alpha.inv <- log.alpha.inv
alphas.vec <- 1/exp(unique(log.alpha.inv)) # We shoud have n.exp.units alphas
elapsed.days <- rep(1,10)

params.list <- list(alpha=alphas.vec[1],k=model12$k.hat)

first.sim <- one.sim(lguess=params.list, elapsed.days=elapsed.days)

nsims <- 5



second.sims <- matrix(0, nrow=length(elapsed.days), ncol=nsims)
for(i in 1:nsims){
  second.sims[,i] <- one.sim(lguess=params.list, elapsed.days=elapsed.days)[,1]
}
max.n <- max(second.sims)


############## Beginning of Supplemental Figure 1 ##############

par(mfrow=c(2,2))

# Supp Fig S1 A
plot(c(0,cumsum(elapsed.days)), c(0,first.sim[,1]), type="s", lwd=2, 
     col="grey", bty="l", ylim=c(0,max.n), xlab="Days", ylab="Cumulative number of embryos")
text(x=0.5,y=80,labels="A.", cex=1.25)

# Supp Fig S1 B

for(i in 1:nsims){
  if(i==1){
    plot(c(0,cumsum(elapsed.days)), c(0,second.sims[,i]), type="s", lwd=2, 
         col="grey", bty="l",xlab="Days", ylim=c(0,max.n), 
         ylab="Cumulative number of embryos")  
  }else{
    points(c(0,cumsum(elapsed.days)), c(0,second.sims[,i]), type="s", lwd=2,col="grey")  
  }
}
text(x=0.5,y=80,labels="B.", cex=1.25)

exp.units <- unique(revisions.data$exp.unit)


# Supp Fig S1 C

revisions.data <- data.frame(read.csv(paste0(data.dir,"single-Ml-spawn-growth_v3.csv")))

# want to show mutliple trajectories that represent a range of fecundity levels
# want to avoid overlap as much as possible for clearer viewing
# 3 seems to be reasonable max to visualize but data for 14 animals available in supp file
# only 5 (aay) goes >250 but plotting that one improves separation

cols <- c("#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#0c2c84")

# Fifth experimental data set
subdat1 <- revisions.data[revisions.data$exp.unit==exp.units[5],] #combos that work: 2,4,5 or 1,2,5 or 1,8,9 or 1,2,9 or 9,10,12
plot(subdat1$day-1,cumsum(subdat1$embryos), type="s", lwd=2,
     col="grey", bty="l", ylim=c(0,250), xlim=c(0,13),
     xlab="Days", ylab="Cumulative number of embryos")

# color vectors for subdat1:
cols.vec1 <- rep(0, nrow(subdat1))
for(i in 1:nrow(subdat1)){
  
  cols.vec1[i] <- cols[subdat1$size[i]]
  
}
points(subdat1$day-1,cumsum(subdat1$embryos),pch=16,col=cols.vec1 )

# Fourth experimental data set
subdat4 <- revisions.data[revisions.data$exp.unit==exp.units[4],] #combos that work: 2,4,5 or 1,2,5 or 1,8,9 or 1,2,9 or 9,10,12
points(subdat4$day-1,cumsum(subdat4$embryos), type="s", lwd=2,
       col="grey", bty="l", ylim=c(0,250), xlim=c(0,13),
       xlab="Days", ylab="Cumulative number of embryos")

# color vectors for subdat1:
cols.vec4 <- rep(0, nrow(subdat4))
for(i in 1:nrow(subdat4)){
  
  cols.vec4[i] <- cols[subdat4$size[i]]
  
}
points(subdat4$day-1,cumsum(subdat4$embryos),pch=16,col=cols.vec4)


# Second experimental data set
subdat2 <- revisions.data[revisions.data$exp.unit==exp.units[2],] #combos that work: 2,4,5 or 1,2,5 or 1,8,9 or 1,2,9 or 9,10,12
points(subdat2$day-1,cumsum(subdat2$embryos), type="s", lwd=2,
       col="grey", bty="l", ylim=c(0,250), xlim=c(0,13),
       xlab="Days", ylab="Cumulative number of embryos")

# color vectors for subdat1:
cols.vec2 <- rep(0, nrow(subdat2))
for(i in 1:nrow(subdat2)){
  
  if(subdat2$size[i]==2.5){subdat2$size[i]<-2}
  cols.vec2[i] <- cols[subdat2$size[i]]
  
}
points(subdat2$day-1,cumsum(subdat2$embryos),pch=16,col=cols.vec2)

plot(0,0, axes=FALSE,type="n",xlab="", ylab="")
legend("topleft", legend=paste0(rep("size = ",length(cols)), 1:length(cols) ), col=cols,
       bty="n", lty=rep(1,length(cols)), lwd=rep(4,length(cols)), cex=1)

######### END Supp Fig 1 #########
