library(tidycensus)
library(tigris)
library(dplyr)
library(tidyr)
library(spdep)
library(RSpectra)
library(Matrix)

## Prep Data
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


## Fit Models
set.seed(1)
iter <- 2000
burn <- 500
source('MCMC.R')


X <- cbind(acs_dat$PovRate, acs_dat$White, acs_dat$Black, acs_dat$Hispanic, acs_dat$Asian)
Yobs <- acs_dat$MedInc
system.time(modFH <- FH_Fit(log(Yobs), X, (acs_dat$MedIncSE)^2/Yobs^2, iter=iter, burn=burn))
system.time(modNN100 <- FH_RNN_Fit(Y=log(Yobs), X=cbind(1,scale(X)), S2=(acs_dat$MedIncSE)^2/Yobs^2, nh=30, iter=iter, burn=burn))


predFH <- rowMeans(exp(modFH$Preds))
predNN100 <- rowMeans(exp(modNN100$Preds))

seFH <- apply(exp(modFH$Preds), 1, sd)
seNN100 <- apply(exp(modNN100$Preds), 1, sd)



##### Plots
library(ggplot2)
library(ggthemes)
PlotDF <- acs_shape 
PlotDF$FH <- predFH; PlotDF$NFH <- predNN100; PlotDF$Direct <- acs_dat$MedInc
PlotDF <- PlotDF %>% pivot_longer(16:18, names_to="Model", values_to="Median Income")

ggplot(PlotDF)+
  geom_sf(linewidth=0, aes(fill=`Median Income`))+
  facet_wrap(~Model, nrow=1)+
  theme_map()+
  scale_fill_viridis_c()
ggsave("CA_ests.jpg", dpi=600)


PlotDF <- acs_shape 
PlotDF$FH <- seFH; PlotDF$NFH <- seNN100; PlotDF$Direct <- acs_dat$MedIncSE
PlotDF <- PlotDF %>% pivot_longer(16:18, names_to="Model", values_to="Standard Error")

ggplot(PlotDF)+
  geom_sf(linewidth=0, aes(fill=`Standard Error`))+
  facet_wrap(~Model, nrow=1)+
  theme_map()+
  scale_fill_viridis_c()
ggsave("CA_se.jpg", dpi=600)



