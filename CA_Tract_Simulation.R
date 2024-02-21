library(tidycensus)
library(tigris)
library(dplyr)
library(tidyr)
library(spdep)
library(RSpectra)
library(Matrix)

## Prepare data
acs_20_vars = load_variables(
  year = 2021, 
  "acs5",
  cache = TRUE
)

keep_vars <- c(
  "B01001_001", #pop size
  "B17001_002", #poverty count
  "B19013_001", # median household income
  "B19058_001"
)

acs_dat <- get_acs(
  geography="tract",
  state="CA",
  year=2021,
  survey="acs5",
  variables=keep_vars
)

acs_dat <- acs_dat %>% pivot_wider(names_from="variable", values_from=c("estimate", "moe")) %>%
  mutate(SNAPRate=estimate_B19058_001/estimate_B01001_001 ,PovRate=estimate_B17001_002/estimate_B01001_001, MedInc=estimate_B19013_001, MedIncSE=(moe_B19013_001/1.645))

acs_shape <- tracts(cb=F, state="CA", year=2021)
acs_shape <- acs_shape %>% filter(GEOID %in% acs_dat$GEOID)
acs_dat <- acs_dat %>% filter(GEOID %in% acs_shape$GEOID)

racevars <- c(White = "P2_005N", 
              Black = "P2_006N", 
              Asian = "P2_008N", 
              Hispanic = "P2_002N")

demographics <- get_decennial(
  geography = "tract",
  state="CA",
  variables = racevars,
  summary_var = "P2_001N",
  year = 2020
) 

demographics <- demographics %>% mutate(pct=value/summary_value) %>%
  dplyr::select(-c(value, summary_value))%>%
  pivot_wider(names_from="variable", values_from="pct")

acs_dat <- acs_dat %>% left_join(demographics, by="GEOID")

acs_dat <- acs_dat %>% filter(!is.na(MedIncSE)) %>% filter(GEOID %in% acs_shape$GEOID)
acs_dat <- acs_dat[order(match(acs_dat$GEOID,acs_shape$GEOID)),]
acs_shape <- acs_shape %>% filter(GEOID %in% acs_dat$GEOID)




## Run Simulation
####################
set.seed(1)
nrep <- 100
iter <- 2000
burn <- 500
source('MCMC.R')


lowFH <- highFH <- lowNN30 <- highNN30 <-lowNN50 <- highNN50 <-lowNN100 <- highNN100 <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)
lowDir <- highDir <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)
predFH <- predDir <- predNN30 <- predNN50 <- predNN100 <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)
xbFH <-  xbNN30 <- xbNN50 <- xbNN100 <- matrix(NA, nrow=nrow(acs_dat), ncol=nrep)

X <- cbind(acs_dat$PovRate, acs_dat$White, acs_dat$Black, acs_dat$Hispanic, acs_dat$Asian)

for(i in 1:nrep){
  fac <- 1
  Yobs <- rnorm(nrow(X), mean=acs_dat$MedInc, sd=fac*acs_dat$MedIncSE)
  Yobs[Yobs < 0] <- 100
  
  modFH <- FH_Fit(log(Yobs), X, (fac*acs_dat$MedIncSE)^2/Yobs^2, iter=iter, burn=burn)
  modNN30 <- FH_RNN_Fit(Y=log(Yobs), X=cbind(1,scale(X)), S2=(fac*acs_dat$MedIncSE)^2/Yobs^2, nh=30, iter=iter, burn=burn)
  modNN50 <- FH_RNN_Fit(Y=log(Yobs), X=cbind(1,scale(X)), S2=(fac*acs_dat$MedIncSE)^2/Yobs^2, nh=50, iter=iter, burn=burn)
  modNN100 <- FH_RNN_Fit(Y=log(Yobs), X=cbind(1,scale(X)), S2=(fac*acs_dat$MedIncSE)^2/Yobs^2, nh=100, iter=iter, burn=burn)
  
  lowFH[,i] <- apply(exp(modFH$Preds), 1, quantile, probs=0.025)
  highFH[,i] <- apply(exp(modFH$Preds), 1, quantile, probs=0.975)
  predFH[,i] <- rowMeans(exp(modFH$Preds))
  xbFH[,i] <- rowMeans(exp(modFH$XB))
  
  
  lowNN30[,i] <- apply(exp(modNN30$Preds), 1, quantile, probs=0.025)
  highNN30[,i] <- apply(exp(modNN30$Preds), 1, quantile, probs=0.975)
  predNN30[,i] <- rowMeans(exp(modNN30$Preds))
  xbNN30[,i] <- rowMeans(exp(modNN30$XB))
  
  lowNN50[,i] <- apply(exp(modNN50$Preds), 1, quantile, probs=0.025)
  highNN50[,i] <- apply(exp(modNN50$Preds), 1, quantile, probs=0.975)
  predNN50[,i] <- rowMeans(exp(modNN50$Preds))
  xbNN50[,i] <- rowMeans(exp(modNN50$XB))
  
  lowNN100[,i] <- apply(exp(modNN100$Preds), 1, quantile, probs=0.025)
  highNN100[,i] <- apply(exp(modNN100$Preds), 1, quantile, probs=0.975)
  predNN100[,i] <- rowMeans(exp(modNN100$Preds))
  xbNN100[,i] <- rowMeans(exp(modNN100$XB))
  
  
  lowDir[,i] <- Yobs - 1.96*fac*acs_dat$MedIncSE
  highDir[,i] <- Yobs + 1.96*fac*acs_dat$MedIncSE
  predDir[,i] <- Yobs
  
  print(i)
}



######### Results


## Coverage rate
cr <- c(mean(acs_dat$MedInc < highDir & acs_dat$MedInc > lowDir),
        mean(acs_dat$MedInc < highFH & acs_dat$MedInc > lowFH),
        mean(acs_dat$MedInc < highNN30 & acs_dat$MedInc > lowNN30),
        mean(acs_dat$MedInc < highNN50 & acs_dat$MedInc > lowNN50),
        mean(acs_dat$MedInc < highNN100 & acs_dat$MedInc > lowNN100))

## Interval score (Gneiting and Raftery 2007)
IS <- c(mean((highDir -lowDir) + 2/.05*(lowDir - acs_dat$MedInc)*(acs_dat$MedInc<lowDir) + 
               2/.05*(acs_dat$MedInc-highDir)*(acs_dat$MedInc>highDir)),
        mean((highFH -lowFH) + 2/.05*(lowFH - acs_dat$MedInc)*(acs_dat$MedInc<lowFH) + 
               2/.05*(acs_dat$MedInc-highFH)*(acs_dat$MedInc>highFH)),
        mean((highNN30 -lowNN30) + 2/.05*(lowNN30 - acs_dat$MedInc)*(acs_dat$MedInc<lowNN30) + 
               2/.05*(acs_dat$MedInc-highNN30)*(acs_dat$MedInc>highNN30)),
        mean((highNN50 -lowNN50) + 2/.05*(lowNN50 - acs_dat$MedInc)*(acs_dat$MedInc<lowNN50) + 
               2/.05*(acs_dat$MedInc-highNN50)*(acs_dat$MedInc>highNN50)),
        mean((highNN100 -lowNN100) + 2/.05*(lowNN100 - acs_dat$MedInc)*(acs_dat$MedInc<lowNN100) + 
               2/.05*(acs_dat$MedInc-highNN100)*(acs_dat$MedInc>highNN100)))

## MSE
mse <- c((mean((acs_dat$MedInc - predDir)^2)),
         (mean((acs_dat$MedInc - predFH)^2)),
         (mean((acs_dat$MedInc - predNN30)^2)),
         (mean((acs_dat$MedInc - predNN50)^2)),
         (mean((acs_dat$MedInc - predNN100)^2)))

mse <- mse/mse[1]

## Abs Bias
ab <- c(mean(abs(acs_dat$MedInc - rowMeans(predDir))),
        mean(abs(acs_dat$MedInc - rowMeans(predFH))),
        mean(abs(acs_dat$MedInc - rowMeans(predNN30))),
        mean(abs(acs_dat$MedInc - rowMeans(predNN50))),
        mean(abs(acs_dat$MedInc - rowMeans(predNN100))))


## Table
resultsTab <- data.frame(Estimator=c("Direct",  "FH", "NN30", "NN50", "NN100"),
                         MSE=mse, Abs_Bias=ab/1000, Cov_Rate=cr, Int_Scor=IS/10000)

library(xtable)


print(xtable(resultsTab, digits=3, 
             caption=paste("Caption")), 
      include.rownames=F)




######### Results (Quantiles)


## Coverage rate
cr <- c(quantile(rowMeans(acs_dat$MedInc < highDir & acs_dat$MedInc > lowDir), probs=c(0.25, 0.75)),
        quantile(rowMeans(acs_dat$MedInc < highFH & acs_dat$MedInc > lowFH), probs=c(0.25, 0.75)),
        quantile(rowMeans(acs_dat$MedInc < highNN30 & acs_dat$MedInc > lowNN30), probs=c(0.25, 0.75)),
        quantile(rowMeans(acs_dat$MedInc < highNN50 & acs_dat$MedInc > lowNN50), probs=c(0.25, 0.75)),
        quantile(rowMeans(acs_dat$MedInc < highNN100 & acs_dat$MedInc > lowNN100), probs=c(0.25, 0.75)))

## Interval score (Gneiting and Raftery 2007)
IS <- c(quantile(rowMeans((highDir -lowDir) + 2/.05*(lowDir - acs_dat$MedInc)*(acs_dat$MedInc<lowDir) + 
               2/.05*(acs_dat$MedInc-highDir)*(acs_dat$MedInc>highDir)), probs=c(0.25, 0.75)),
        quantile(rowMeans((highFH -lowFH) + 2/.05*(lowFH - acs_dat$MedInc)*(acs_dat$MedInc<lowFH) + 
               2/.05*(acs_dat$MedInc-highFH)*(acs_dat$MedInc>highFH)), probs=c(0.25, 0.75)),
        quantile(rowMeans((highNN30 -lowNN30) + 2/.05*(lowNN30 - acs_dat$MedInc)*(acs_dat$MedInc<lowNN30) + 
               2/.05*(acs_dat$MedInc-highNN30)*(acs_dat$MedInc>highNN30)), probs=c(0.25, 0.75)),
        quantile(rowMeans((highNN50 -lowNN50) + 2/.05*(lowNN50 - acs_dat$MedInc)*(acs_dat$MedInc<lowNN50) + 
               2/.05*(acs_dat$MedInc-highNN50)*(acs_dat$MedInc>highNN50)), probs=c(0.25, 0.75)),
        quantile(rowMeans((highNN100 -lowNN100) + 2/.05*(lowNN100 - acs_dat$MedInc)*(acs_dat$MedInc<lowNN100) + 
               2/.05*(acs_dat$MedInc-highNN100)*(acs_dat$MedInc>highNN100)), probs=c(0.25, 0.75)))

## MSE
mse <- c(quantile(rowMeans((acs_dat$MedInc - predDir)^2), probs=c(0.25, 0.75)),
         quantile(rowMeans((acs_dat$MedInc - predFH)^2), probs=c(0.25, 0.75)),
         quantile(rowMeans((acs_dat$MedInc - predNN30)^2), probs=c(0.25, 0.75)),
         quantile(rowMeans((acs_dat$MedInc - predNN50)^2), probs=c(0.25, 0.75)),
         quantile(rowMeans((acs_dat$MedInc - predNN100)^2), probs=c(0.25, 0.75)))

mse <- round(mse/mean((acs_dat$MedInc - predDir)^2), 2)

## Abs Bias
ab <- c(quantile((abs(acs_dat$MedInc - rowMeans(predDir))), probs=c(0.25, 0.75)),
        quantile((abs(acs_dat$MedInc - rowMeans(predFH))), probs=c(0.25, 0.75)),
        quantile((abs(acs_dat$MedInc - rowMeans(predNN30))), probs=c(0.25, 0.75)),
        quantile((abs(acs_dat$MedInc - rowMeans(predNN50))), probs=c(0.25, 0.75)),
        quantile((abs(acs_dat$MedInc - rowMeans(predNN100))), probs=c(0.25, 0.75)))


## Table
resultsTabQ <- data.frame(Estimator=c("Direct",  "FH", "NN30", "NN50", "NN100"),
                         MSE=mse, Abs_Bias=ab/1000, Cov_Rate=cr, Int_Scor=IS/10000)

library(xtable)


print(xtable(resultsTab, digits=3, 
             caption=paste("Caption")), 
      include.rownames=F)
