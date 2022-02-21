
source("NegBinTools.R")

library(gam)

data.dir <- paste0("/Users/josemiguel/Dropbox/AllisonEdgar/","Data/")
#######################  October 29 test with made up dataset ########

# Comparison of mean jump for two different models:
# First upoading the data: 
test.data <-data.frame(read.csv("tester.csv"))


# Model 1 setting:  ~number
mod1.test <- negbinom.fit(formula=~treatment+bio.rep, full.df=test.data, 
                          exp.unit=test.data$exp.unit)


mod2.test <- negbinom.fit(formula=~number, full.df=test.data, 
                          exp.unit=test.data$exp.unit)

mod3.test <- negbinom.fit(formula=~bio.rep, full.df=test.data, 
                          exp.unit=test.data$exp.unit)



mod4.test <- negbinom.fit(formula=~treatment, full.df=test.data, 
                          exp.unit=test.data$exp.unit)



# Prediction for model 1:

#mod2.pred <-negbinomfit.pred(betas=mod2.test$betas.hat, k=mod2.test$k.hat, 
#                             formula=~number, exp.unit=temp.dataf$exp.unit, 
#                             full.df=test.data)

#mod2.pred$Expected/mod2.pred$number


# Comparison of mean jump for two different models:
# First upoading the data: 
temperature.fname <- paste0(data.dir,"temp-fecundity-revised.csv") 
temp.dataf <-read.csv(temperature.fname)


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

### To compare the average number of embryos in two models do:

mod1.pred[,13]/mod1.pred[,4] # which gives the mean number of embryos 
                                     # per parental always around 2.05 regardless
                                     # of temperature
#> mod1.pred[1:14,13]/mod1.pred[1:14,4]
#[1] 2.058991 2.058991 2.058991 2.058991 2.058991 2.044011 2.037098 2.081283 2.081283 2.081283 2.081283 2.081283 2.081283 2.037098

# And for model 2, where a temperature effect is included, the mean number of 
# embryos per parental at lower temperature (20) is around 0.88 whereas the
# mean number of embryos per parental at high temperature (28) is around 3.53:

#mod2.pred[1:14,13]/mod2.pred[1:14,4]
#[1] 0.8815909 0.8815909 0.8815909 0.8815909 0.8815909 0.8625784 0.8472857 3.5290316 3.5290316 3.5290316 3.5290316 3.5290316 3.5290316 3.3070781


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


# ### Initial figure submitted with the first manuscript version (December 2021):  
# ### actually not true... this has been modified to produce something else at some point...
# two.cols <- c("#56B4E9", "#E69F00") 
# two.cols.light <- c("#9ad2f2", "#ffc034")
# # Figure for paper
# par(oma=c(1,1,1,1), mar=c(4,4,1,1))
# plot(num.vec, pred.lm20,type="l", lwd=2, 
#      ylab="Per capita fecundity", xlab="Number of parentals", xlim=c(6.5,21), ylim=c(0,7),
#      bty="l")
# polygon(x=c(num.vec,rev(num.vec)), y=c(pred.lm20-propag.relses*pred.lm20,
#                                        rev(pred.lm20+propag.relses*pred.lm20)), col=two.cols.light[1], border=NA)
# points(num.vec,pred.lm20, type="l", lwd=2, col=two.cols[1])
# 
# points(num.vec, pred.lm28,type="l", lwd=2, 
#      ylab="Per capita fecundity", xlab="number of parentals", xlim=c(5,21), ylim=c(3,100),
#      bty="l")
# polygon(x=c(num.vec,rev(num.vec)), y=c(pred.lm28-propag.relses*pred.lm28,
#                                        rev(pred.lm28+propag.relses*pred.lm28)), col=two.cols.light[2], border=NA)
# points(num.vec,pred.lm28, type="l", lwd=2, col=two.cols[2])
# legend("topleft", legend=c("20 degrees", "28 degrees"), col=two.cols, lty=1,lwd=2, bty="n")
# # End figure for paper v1 (initial submission)


#ses.fit20 <- predict(lm20,newdata=df.4pred, se.fit=TRUE)$se.fit
#ses.fit28 <- predict(lm28,newdata=df.4pred, se.fit=TRUE)$se.fit

#### Transparency rgb function:
## Transparent colors
## Mark Gardener 2015
## www.dataanalytics.org.uk



########  Updated Figure for the paper, done during revisions (Feb 14 2022)

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
two.cols.light <- c(rgb(154, 210, 242, max=255,alpha=),rgb(255, 192, 52))

# Figure for paper
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







# 
# plot(temp.f28$number, temp.f28$model2.pred, pch=16, col="blue", 
#      ylab="fecundity", xlab="number of parentals", xlim=c(5,21), ylim=c(3,90),
#      bty="l")
# points(temp.f20$number, temp.f20$model2.pred, pch=16, col="red", 
#      ylab="fecundity", xlab="number of parentals")





br <- unique(temp.fecund.preds$bio.rep)
nbr <- length(br)
par(mfrow=c(2,2), oma=c(1,1,1,1), mar= c(4,4,3,1))
for(i in 1:nbr){
  
  subdat <- temp.fecund.preds[temp.fecund.preds$bio.rep==br[i],]
  plot(subdat$number,subdat$model2.pred/subdat$number, pch=16, col="red",
       ylab="per capita fecundity", xlab="number of parentals",
       main=paste0("per capita fecundity for bio. rep ",br[i]),
       bty="l")
  points(subdat$number,subdat$model1.pred/subdat$number, pch=16, col="blue")
  
}
  

####################  Next comparison: ~bio.rep vs ~bio.rep+number 

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








# We are going to compare the predicted average number of embryos without
# taking into account number of parentals vs. the average number of embryos
# when the parentals are 1,2,4, and 8

num.vec <- unique(dens1.data$number) # 1,2,4,8
means.Exp <- rep(0,4)
list.preds.dens1 <- list()

for(i in 1:4){
  
  subdat.i   <- mod2.pred[mod2.pred$number==num.vec[i],]
  Exp4allbioreps <- unique(subdat.i$Expected)/num.vec[i]
  means.Exp[i] <- mean(Exp4allbioreps)
  
  percapfecmat <- cbind(rep(num.vec[i],length(Exp4allbioreps)),Exp4allbioreps)
  
  list.preds.dens1[[i]] <- percapfecmat
  
  if(i==1){percap.fec.dens1 <- percapfecmat}else{
    
    percap.fec.dens1 <- rbind(percap.fec.dens1,percapfecmat)
  }
  
}


par(oma=c(1,1,1,1), mar=c(5,5,2,1))
plot(percap.fec.dens1[,1], percap.fec.dens1[,2], pch=16, bty="l", xlab="Number of parentals", ylab="Per capita Fecundity", cex.lab=1.55)

library(ggplot2)
library(gam)


qplot(x=percap.fec.dens1[,1],y= percap.fec.dens1[,2], geom='auto', ylim=c(0,0.75), xlab="number of parentals", 
            ylab="Per capita Fecundity")

fit3 <- gam(percap.fec.dens1[,2]~s(percap.fec.dens1[,1], df=6))
pred3 <- predict(fit3, se.fit=TRUE, level=0.90)


plot(percap.fec.dens1[,1],pred3$fit, type="l", col="red", lwd=2)
points(percap.fec.dens1[,1],pred3$fit+ pred3$se.fit, type="l", col="gray", lwd=2)
points(percap.fec.dens1[,1],pred3$fit- pred3$se.fit, type="l", col="gray", lwd=2)




data.dir <- paste0("/Users/josemiguel/Dropbox/AllisonEdgar/","Data/")
#### Case 3:
  
#### The data was:
  
dens2.data <- read.csv(paste0(data.dir,"dissogeny_density_2.csv"))

#### And tested models were

mod1.test <-  negbinom.fit(formula=~bio.rep, full.df=dens2.data, 
                           exp.unit = dens2.data$exp.unit)

mod1.pred <- negbinomfit.pred(betas=mod1.test$betas.hat, k=mod1.test$k.hat, 
                              formula=~bio.rep, exp.unit=dens2.data$exp.unit, 
                              full.df=dens2.data, response="embryos")

#### vs.

mod2.test <- negbinom.fit(formula=~bio.rep+number, full.df=dens2.data, 
                          exp.unit = dens2.data$exp.unit)

mod2.pred <- negbinomfit.pred(betas=mod2.test$betas.hat, k=mod2.test$k.hat, 
                              formula=~bio.rep+number, 
                              exp.unit=dens2.data$exp.unit, full.df=dens2.data, 
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
br.4df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R4",] #"R3.5613" "R4.5624" "R1RG"    "R2RG"

br.5df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R3.5613",]
br.6df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R4.5624",]
br.7df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R1RG",]
br.8df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R2RG",] 


lm1 <- lm(model2.pred~number+I(number^2), data=br.1df)
lm2 <- lm(model2.pred~number+I(number^2), data=br.2df)
lm3 <- lm(model2.pred~number+I(number^2), data=br.3df)
lm4 <- lm(model2.pred~number+I(number^2), data=br.4df)

lm5 <- lm(model2.pred~number+I(number^2), data=br.5df)
lm6 <- lm(model2.pred~number+I(number^2), data=br.6df)
lm7 <- lm(model2.pred~number+I(number^2), data=br.7df)
lm8 <- lm(model2.pred~number+I(number^2), data=br.8df)


df.4pred <- data.frame(number=1:20)
num.vec <- df.4pred$number
pred1 <- predict(lm1, newdata=df.4pred)
pred2 <- predict(lm2, newdata=df.4pred)
pred3 <- predict(lm3, newdata=df.4pred)
pred4 <- predict(lm4, newdata=df.4pred)

pred5 <- predict(lm5, newdata=df.4pred)
pred6 <- predict(lm6, newdata=df.4pred)
pred7 <- predict(lm7, newdata=df.4pred)
pred8 <- predict(lm8, newdata=df.4pred)



# Sensible cols
my.cols  <- c("#a1dab4", "#41b6c4", "#2c7fb8", "#253494")
light.cols <- c("#c6e8d1", "#7bccd6", "#59a4d7","#384ccd")

par(oma=c(1,1,1,1), mar=c(4,4,1,1))
# plot for R1
plot(num.vec, pred1, type="l", lwd=2, ylab="Fecundity", xlab="Number of parentals",
     xlim=c(0,21), ylim=c(0,155), bty="l")
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








########################### Case 4: #################################
feed2.data <- read.csv(paste0(data.dir,"feed2.csv"))


mod1.test <-  negbinom.fit(formula=~bio.rep+DHA, full.df=feed2.data, 
                           exp.unit = feed2.data$exp.unit)



#### vs.

mod2.test <- negbinom.fit(formula=~bio.rep+number+DHA, full.df=feed2.data, 
                          exp.unit = feed2.data$exp.unit)

mod2.pred <- negbinomfit.pred(betas=mod2.test$betas.hat, k=mod2.test$k.hat, 
                              formula=~bio.rep+number+DHA,exp.unit=feed2.data$exp.unit, 
                              full.df=feed2.data, response="embryos")

####
dens1.2modpreds <- data.frame(number=mod2.pred$number, model2.pred=mod2.pred[,22], 
                              bio.rep=mod2.pred$bio.rep, DHA=mod2.pred[,11], percap.pred=mod2.pred[,22]/mod2.pred$number)

betas.mod2 <- c(mod2.test$betas.hat, mod2.test$k.hat)
betasses.mod2 <- mod2.test$betas.SEs

prop.error <- mean(betasses.mod2/abs(betas.mod2))

br.1df <- dens1.2modpreds[(dens1.2modpreds$bio.rep=="R1"),]
br.2df <- dens1.2modpreds[dens1.2modpreds$bio.rep=="R2",]

lm1 <- lm(percap.pred~DHA+I(DHA^2), data=br.1df)
lm2 <- lm(model2.pred~DHA+I(DHA^2), data=br.2df)

df.4pred <- data.frame(DHA=seq(from=0.3, to=42,by=0.1))
num.vec <- df.4pred$DHA
pred1 <- predict(lm1, newdata=df.4pred)
pred2 <- predict(lm2, newdata=df.4pred)

# Sensible cols
my.cols  <- c("#a1dab4", "#41b6c4", "#2c7fb8", "#253494")
light.cols <- c("#c6e8d1", "#7bccd6", "#59a4d7","#384ccd")

par(mfrow=c(1,2),oma=c(1,1,1,1), mar=c(4,4,1,1))
plot(num.vec, pred2, type="l", lwd=2, ylab="Fecundity per capita", xlab="DHA",
     xlim=c(0,42), ylim=c(0.5,70), bty="l", main="R1", cex.main=1.5)
polygon(x=c(num.vec,rev(num.vec)), y=c(pred2-prop.error*pred2, rev(pred2+prop.error*pred2)), 
        col=light.cols[2], border=NA)
points(num.vec, pred2, type="l", lwd=2, col=my.cols[2])


plot(num.vec, pred1, type="l", lwd=2, ylab="Fecundity per capita", xlab="DHA",
     xlim=c(0,42), ylim=c(0.5,70), bty="l", main="R2", cex.main=1.5)
polygon(x=c(num.vec,rev(num.vec)), y=c(pred1-prop.error*pred1, rev(pred1+prop.error*pred1)), 
        col=light.cols[1], border=NA)
points(num.vec, pred1, type="l", lwd=2, col=my.cols[1])



