library(forecast)
library(gmwm)
# Set seed for reproducibility
set.seed(123)


source("WKestimator.R")
source("tbatsHighFreq.R")
source("validation.R")

c <- 0.8 
d <- 0.7282 
a <- -0.9   # MA part negated to be consistent with +signs in the filter defn.
b <- 0.0242 # MA part negated to be consistent with +signs in the filter defn.

m <- 17989.4920 # mean of the NY energy data


mod1 <- list(sigma2 = 10282, coef = c(c,-a,d,-b), mean = m)

MSE_PA <- 1 + (1 - a + c)^2 + (1 - a + c + (-a*c + c^2))^2 + (1 - a + c + (-a*c +c^2) 
                                                              + (-a*c^2 + c^3))^2

denom <- (-a*((1+c)^2)*(1+4*c^2 + c^4) + (1 + a^2)*(c + 2*c^2 + 4*c^3 + 2*c^4 + c^5))
num <- (-a*((1 + c)^2)*(6 + 8*c^2 + 6*c^4) + (1 + a^2)*
          (4 + 6*c + 8*c^2 + 8*c^3 + 8*c^4 + 6*c^5 + 4*c^6))
MSE_AP <- 2*(denom^2)/(num - sqrt(num^2 - 4*(denom^2)))


# generate
modelARIMA <- SARIMA(ar= c, i = 0, ma = -a, sar = d, si = 0, sma = -b, s = 4, sigma2 = mod1$sigma2)
N <- 4*24*(365+60)

yt <- gen_gts(N, modelARIMA) + m



Yt <- ts(yt, freq = 4)
n <- 365*24*4
N <- (365+31)*24*4

Yactual <- ts(colSums(matrix(Yt,4,length(Yt)/4)))

#-----------
# WK
#-----------
Nwk <- (36500+31)*24*4
nwk <- 36500*24*4
yt_wk <- gen_gts(Nwk, modelARIMA) + m
Yactual_wk <- ts(colSums(matrix(yt_wk,4,length(yt_wk)/4)))
startQ <- Nwk - (3650+31)*24*4 + 1
startH <- Nwk/4 - (3650+31)*24 + 1
errPAwk <- WKestimatorPA(mod1, yt_wk[startQ:Nwk]) 
errAPwk <- WKestimatorAP(mod1, Yactual_wk[startH:(Nwk/4)])

fin_wkPAm <- Acc(errPAwk[1:(3650*24)], Yactual_wk[startH:(nwk/4)]) 
fin_wkAPm <- Acc(errAPwk[1:(3650*24)], Yactual_wk[startH:(nwk/4)])
fin_wkPAh <- Acc(errPAwk[1:(31*24) + (3650*24)], Yactual_wk[(nwk/4 + 1):(Nwk/4)])
fin_wkAPh <- Acc(errAPwk[1:(30*24) + (3650*24)], Yactual_wk[(nwk/4 + 1):(Nwk/4 - 24)])


#-----------
# TBATS
#------------
errAPTBATS_sim2 <- list()
seasAP <- 24 
for (i in c(1:31)){
  print(i)
  errAPTBATS_sim2[[i]] <- tbatsHighFreq(Yactual[(1:(366*24)) + ((i-1)*24)], 1, seasAP, 24, 365*24)
}
save.image(file = "simTBATS_AP_2.RData")

seasPA <- 4
errPATBATS_sim2 <- list()
for (i in c(1:31)){
  print(i)
  errPATBATS_sim2[[i]] <- tbatsHighFreq(Yt[(1:(96*366))+((i-1)*96)], 1, seasPA, 96, 365*96)
}
save(errPATBATS_sim2, file="simTBATS_PA_2")

# Errors  
#load("simTBATS_AP_2.RData")
#load("simTBATS_PA_2")
errModPA <- errPATBATS_sim2[[1]]$model[1:(8760-24)]
errModAP <- errAPTBATS_sim2[[1]]$model[1:(8760-24)]
errHoldPA <- matrix(NA,31,24)
errHoldAP <- matrix(NA,31,24)
for (i in 1:31){
  errModPA <- c(errModPA,errPATBATS_sim2[[i]]$model[(8760-24+1):8760])
  errModAP <- c(errModAP,errAPTBATS_sim2[[i]]$model[(8760-24+1):8760])
  errHoldPA[i,] <- errPATBATS_sim2[[i]]$hold[1:24]
  errHoldAP[i,] <- errAPTBATS_sim2[[i]]$hold[1:24]
}

fin_tbatsPAm <- Acc(errModPA, Yactual[(1:(395*24))]) #365 + 31 days
fin_tbatsAPm <- Acc(errModAP, Yactual[(1:(395*24))])
fin_tbatsPAh <- Acc(as.vector(t(errHoldPA)), Yactual[(365*24)+(1:(31*24))])
fin_tbatsAPh <- Acc(as.vector(t(errHoldAP)), Yactual[(365*24)+(1:(31*24))])

write.csv(rbind(c(fin_wkPAm, fin_wkAPm, fin_wkPAh, fin_wkAPh),
                c(fin_tbatsPAm, fin_tbatsAPm, fin_tbatsPAh, fin_tbatsAPh)), 
          file = "sim2.csv")
