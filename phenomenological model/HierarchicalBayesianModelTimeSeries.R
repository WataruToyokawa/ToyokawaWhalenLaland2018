# Hierarchical Bayesian Model analysis -- the moderately difficult condition
# The code generating Fig. 3a and 3b

# Loading packages
library(tidyverse)
library(rstan)

# parallel computing when you run multiple chains in a MCMC simulation
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

# Loading data
allBehaviouralData2 <- read_csv("allBehaviouralData2.csv")
allBehaviouralData2 = allBehaviouralData2 %>% dplyr::mutate(T2 = allBehaviouralData2$round)
allBehaviouralData2$indiv_group = 0
allBehaviouralData2$indiv_group[which(allBehaviouralData2$sizeCategory>1)] = 1
allBehaviouralData2$room_indiv_group = as.character(allBehaviouralData2$room)
allBehaviouralData2$room_indiv_group[which(allBehaviouralData2$sizeCategory==1)] = 'Single-condition'
allBehaviouralData_25 = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%')
allBehaviouralData_50 = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%')
allBehaviouralData_78 = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%')

allBehaviouralData_50$ID_caracter = as.character(allBehaviouralData_50$ID)

## ========= 50% (moderate difficulty condition) ==================================
## ========= (50% means that 50% area of the payoff distributions are overlapped) =
##

# preparing data fed to Stan
onlineExp_stanData_50 = list(
  All = allBehaviouralData_50 %>% nrow(),
  T = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::pull(T2),
  maxT = 70,
  Group = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%') %>% dplyr::pull(room) %>% as.factor() %>% as.numeric(),
  Indv = allBehaviouralData_50 %>% dplyr::pull(ID_caracter) %>% as.factor() %>% as.numeric(),
  indiv_group = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::pull(indiv_group),
  #GroupSize_st = (allBehaviouralData_50$sizeCategory-mean(allBehaviouralData_50$sizeCategory))/sd(allBehaviouralData_50$sizeCategory),
  N_group = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::select(room) %>% dplyr::count(room) %>% spread(room, n) %>% length(),
  N_indv = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::select(amazonID) %>% dplyr::count(amazonID) %>% spread(amazonID, n) %>% length(),
  Y = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::pull(optimalChoice)
  )

# incerting payoff 
onlineExp_stanData_50$payoff = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%') %>% dplyr::pull(result)*100

# Standardizing group size
onlineExp_stanData_50$GroupSize_st = (
  ( (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%') %>% dplyr::pull(sizeCategory)) -
    (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))

# Calculating three standardised group sizes of three different categories 
onlineExp_stanData_50$smallGroupSize = (
  ( 4 -
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))
onlineExp_stanData_50$largeGroupSize = (
  ( 13 -
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))
onlineExp_stanData_50$veryLargeGroupSize = (
  ( 20 -
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))

## Compiling the Stan model
stanmodel_onlineExp = stan_model(file='./timeSeriesDiffTest_optimRate_groupSize.stan')

## Fitting
onlineExp_fitting_groupSize_50 = sampling(
  stanmodel_onlineExp, data=onlineExp_stanData_50, seed=72,
  #control = list(adapt_delta = 0.80),
  pars=c('sizeEffect','s_i','mu0','di0','mu','di','p_indiv','p_group','p_group_small','p_group_large','p_group_veryLarge'),
  init=function() {
    list(mu0=runif(1,-3,-1))
  },
  #chains=8, iter=2000, warmup=1000, thin=5
  chains=1, iter=200, warmup=100, thin=1 # debug
)

## Fitting result
onlineExp_fitting_groupSize_50

# save MCMC result (you can reuse it for the future analysis)
saveRDS(onlineExp_fitting_groupSize_50, file='./onlineExp_fitting_groupSize_50.rds')

trace_onlineExp_fitting_groupSize_50 = traceplot(onlineExp_fitting_groupSize_50, pars=c('sizeEffect', 's_i', 'mu0', 'di0'), inc_warmup=TRUE)
ggsave(file = "./trace_onlineExp_fitting_groupSize_50.pdf", plot = trace_onlineExp_fitting_groupSize_50, dpi = 600, width = 9, height = 9)

#plot(onlineExp_fitting_groupSize_50, pars=c('p_indiv','p_group_small'))
plot(onlineExp_fitting_groupSize_50, pars=c('mu','di'))
plot(onlineExp_fitting_groupSize_50, pars=c('sizeEffect','s_i'))

# Summarising the MCMC sample
data_onlineExp_fitting_groupSize_50 <- rstan::extract(onlineExp_fitting_groupSize_50)
fittedValues_onlineExp_fitting_groupSize_50 =
  as.data.frame(t(apply(data_onlineExp_fitting_groupSize_50$mu[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_50$di[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_50$p_indiv[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_50$p_group_small[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_50$p_group_large[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_50$p_group_veryLarge[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_50$p_group[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
## rename
colnames(fittedValues_onlineExp_fitting_groupSize_50) = c("mu_p2.5", "mu_p25", "mu_p50", "mu_p75", "mu_p97.5","di_p2.5", "di_p25", "di_p50", "di_p75", "di_p97.5","p_indiv_p2.5", "p_indiv_p25", "p_indiv_p50", "p_indiv_p75", "p_indiv_p97.5","p_group_small_p2.5", "p_group_small_p25", "p_group_small_p50", "p_group_small_p75", "p_group_small_p97.5","p_group_large_p2.5", "p_group_large_p25", "p_group_large_p50", "p_group_large_p75", "p_group_large_p97.5","p_group_veryLarge_p2.5", "p_group_veryLarge_p25", "p_group_veryLarge_p50", "p_group_veryLarge_p75", "p_group_veryLarge_p97.5","p_group_p2.5", "p_group_p25", "p_group_p50", "p_group_p75", "p_group_p97.5")
fittedValues_onlineExp_fitting_groupSize_50$taskDifficulty = "50%"
fittedValues_onlineExp_fitting_groupSize_50$round = rep(1:70)

# Summarising the model prediction for each group size category (i.e. large and small group)
group_level_summary_small_50 = tapply(
	subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=6&taskDifficulty=='50%')$optimalChoice, 
	list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=6&taskDifficulty=='50%')$round,
		subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=6&taskDifficulty=='50%')$room_indiv_group
		), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_large_50 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>6&taskDifficulty=='50%')$optimalChoice, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>6&taskDifficulty=='50%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory>6&taskDifficulty=='50%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)

group_level_summary = data.frame(
  round = rep(1:70),
  taskDifficulty3 = rep('Moderate Uncertainty',70),
  small = group_level_summary_small_50,
  large = group_level_summary_large_50
  )

## plot (Fig. 3a -- Moderate condition)
(timeseries_optimrate_analisys_groupSize_plot = ggplot() +
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, ymin=p_indiv_p2.5, ymax=p_indiv_p97.5),fill='black', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, ymin=p_indiv_p25, ymax=p_indiv_p75),fill='black', alpha=3/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, ymin=p_group_small_p2.5, ymax=p_group_small_p97.5),fill='orange', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, ymin=p_group_small_p25, ymax=p_group_small_p75),fill='orange', alpha=3/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, ymin=p_group_large_p2.5, ymax=p_group_large_p97.5),fill='red', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, ymin=p_group_large_p25, ymax=p_group_large_p75),fill='red', alpha=3/6)+
  stat_summary(data=subset(allBehaviouralData2, indiv_group==0), mapping=aes(x=round,optimalChoice), fun.y = mean, geom="line",colour='black',size=1,alpha=3/4)+
  geom_line(data=group_level_summary, mapping=aes(round, small), colour='orange', size=1, alpha=3/4)+
  geom_line(data=group_level_summary, mapping=aes(round, large), colour='red', size=1, alpha=3/4)+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, y=p_indiv_p50),colour='black',linetype='dashed')+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, y=p_group_small_p50),colour='orange',linetype='dashed')+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, y=p_group_large_p50),colour='red',linetype='dashed')+
  #facet_grid(.~taskDifficulty3)+
  scale_y_continuous(limits = c(-0.05,1.05), breaks=c(0.0,0.25,0.50,0.75,1.0)) +
  labs(x='Rounds', y='Proportion \nchoosing\n the best option')+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  geom_hline(yintercept=1/3, linetype="dashed")+
  theme_classic()+
  NULL)

## plot (Fig. 3b -- Moderate condition)
(difference_btwn_indiv_vs_group_plot = ggplot() +
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, ymin=di_p2.5, ymax=di_p97.5),fill='black', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, ymin=di_p25, ymax=di_p75),fill='black', alpha=3/6)+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize_50, mapping=aes(x=round, y=di_p50),colour='black')+
  labs(x='Rounds', y=expression(paste('Differences between\nsolitaries and groups ',xi[t],sep="")))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  geom_hline(yintercept=0, linetype="dashed")+
  theme_classic()+
  NULL)

