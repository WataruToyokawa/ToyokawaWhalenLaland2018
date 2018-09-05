##########################################################
## Figures in a paper entitled:
## "Social learning strategies regulate the wisdom and madness of interactive crowds"
## by Wataru Toyokawa*, Andrew Whalen and Kevin Laland
## School of Biology, University of St Andrews, Fife, Scotland
## *corresponding author: wt25@st-andrews.ac.uk (@WataruToyokawa)
##########################################################
# README
# This Script firstly loads all data and MCMC results as .rds files, and then depicts figures.
# You can re-generate all figures in the main text by just doing copy-and-paste the following script into R console.
#
# The analysis code (e.g. running a MCMC) appear after the figure-depiction scripts, where you can re-analyse the models.

# Functions
myTheme = function() {
  theme(
    legend.position = 'none',
    legend.title = element_text(size=12, family="Helvetica" ,colour='black'),
    legend.text = element_text(size=12, family="Helvetica" ,colour='black'),
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=11, family="Helvetica", vjust=1 ), #"Helvetica"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=10, family="Helvetica" ,colour='black'),
    axis.text.y = element_text(size=10, family="Helvetica" ,colour='black'),
    axis.title.x=element_text(size=12, family="Helvetica" ),
    axis.title.y=element_text(size=12, family="Helvetica" ),
    panel.background = element_rect(fill = "white", colour = NA),
    #panel.grid = element_blank(),
    panel.grid = element_line(colour = "black", size = 0.8),
    plot.title = element_text(size=12, family="Helvetica" ),
    plot.background = element_rect(colour = "white")
  )
}
myTheme_legend = function() {
  theme(
    #legend.position = 'right',
    legend.title = element_text(size=12, family="Helvetica" ,colour='black'),
    legend.text = element_text(size=12, family="Helvetica" ,colour='black'),
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=11, family="Helvetica", vjust=1 ), #"Helvetica"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=10, family="Helvetica" ,colour='black'),
    axis.text.y = element_text(size=10, family="Helvetica" ,colour='black'),
    axis.title.x=element_text(size=12, family="Helvetica" ),
    axis.title.y=element_text(size=12, family="Helvetica" ),
    panel.background = element_rect(fill = "white", colour = NA),
    #panel.grid = element_blank(),
    panel.grid = element_line(colour = "black", size = 0.8),
    plot.title = element_text(size=12, family="Helvetica" ),
    plot.background = element_rect(colour = "white")
  )
}
myTheme_inner = function() {
  theme(
    #legend.position = 'none',
    legend.title = element_text(size=12, family="Helvetica" ,colour='black'),
    legend.text = element_text(size=12, family="Helvetica" ,colour='black'),
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 1,  linetype="solid"),
    strip.text = element_text(size=11, family="Helvetica" ), #"Helvetica"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=10, family="Helvetica" ,colour='black'),
    axis.text.y = element_text(size=10, family="Helvetica" ,colour='black'),
    axis.title.x=element_text(size=12, family="Helvetica" ),
    axis.title.y=element_text(size=12, family="Helvetica" ),
    axis.line =element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_line(colour = "black", size = 0.8),
    plot.title = element_text(size=12, family="Helvetica" ),
    plot.background = element_rect(colour = "white")
  )
}
WAIC = function(fit) {
  log_lik = rstan::extract(fit)$log_lik
  lppd = sum(log(colMeans(exp(log_lik))))
  p_waic = sum(apply(log_lik, 2, var))
  waic = 2*(-lppd + p_waic)
  return(list(waic=waic, lppd=lppd, p_waic=p_waic))
}
WAIC_indv = function(fit) {
  log_lik = rstan::extract(fit)$log_lik
  lppd = log(colMeans(exp(log_lik)))
  p_waic = apply(log_lik, 2, var)
  waic = 2*(-lppd + p_waic)
  return(list(waic=waic, lppd=lppd, p_waic=p_waic))
}

# Loading stuff
library(readr)
library(cowplot)
library(ggjoy)
library(rstan)
library(ggmcmc)
library(viridisLite)
library(viridis)
library(stringr)
library(lme4)
library(data.table)
library(rstan)

## parallel computing when you run multiple chains in a MCMC simulation
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
# Loading data (.csv)
UNC_Moderate_expParameters_eachGroup <- read_csv("./data/UNC_Moderate_expParameters_eachGroup.csv")
UNC_Moderate_expParameters_timeSeries0 <- read_csv("./data/UNC_Moderate_expParameters_timeSeries.csv")
allBehaviouralData2 <- read_csv("./data/allBehaviouralData2.csv")
performanceSummary4 <- read_csv("./data/performanceSummary4.csv")
copyingProbReduction_data_fullModelOnly <- read_csv("./data/copyingProbReduction_data_fullModelOnly.csv")
questionnaire_data_exp <- read_csv("./data/questionnaire_data_exp.csv")
temp_UNC6_sReduc_annealing <- read_csv("./data/temp_UNC6_sReduc_annealing.csv")
temp_AL6_annealing_indiv <- read_csv("./data/temp_AL6_annealing_indiv.csv")
temp_UNC6_sReduc_annealing_simulation <- read_csv("./data/model6_simulation/temp_UNC6_sReduc_annealing_simulation.csv")
soc_table_UNC6_sReduc_annealing <- read_csv("./data/soc_table_UNC6_sReduc_annealing.csv")
netBeta_table_UNC6_sReduc_annealing <- read_csv("./data/netBeta_table_UNC6_sReduc_annealing.csv")
netBeta_table_AL6_annealing_indiv <- read_csv("./data/netBeta_table_AL6_annealing_indiv.csv")
UNC_Moderate_otherParam1_eachGroup_sup1 <- read_csv("./data/UNC_Moderate_otherParam1_eachGroup_sup1.csv")
UNC_Moderate_otherParam1_timeSeries_sup1 <- read_csv("./data/UNC_Moderate_otherParam1_timeSeries_sup1.csv")
UNC_Moderate_otherParam2_eachGroup_sup1 <- read_csv("./data/UNC_Moderate_otherParam2_eachGroup_sup1.csv")
UNC_Moderate_otherParam2_timeSeries_sup1 <- read_csv("./data/UNC_Moderate_otherParam2_timeSeries_sup1.csv")
posthoc_sim_summary_indivSmallLarge <- read_csv("./data/stan_mdelFitting/posthoc_sim_summary_indivSmallLarge.csv")
posthoc_sim_eachGroup_aveOptim <- read_csv("./data/stan_mdelFitting/posthoc_sim_eachGroup_aveOptim.csv")
posthoc_sim_summary_all <- read_csv("./data/stan_mdelFitting/posthoc_sim_summary_all.csv")
# Loading analysis results (.rds)
onlineExp_fitting_groupSize_25 <- readRDS("./rds/mcmc_result/onlineExp_fitting_groupSize_25.rds")
onlineExp_fitting_groupSize_50 <- readRDS("./rds/mcmc_result/onlineExp_fitting_groupSize_50.rds")
onlineExp_fitting_groupSize_78 <- readRDS("./rds/mcmc_result/onlineExp_fitting_groupSize_78.rds")
onlineExp_fitting_groupSize_supp1_25 <- readRDS("./rds/mcmc_result/onlineExp_fitting_groupSize_supp1_25.rds")
onlineExp_fitting_groupSize_supp1_50 <- readRDS("./rds/mcmc_result/onlineExp_fitting_groupSize_supp1_50.rds")
onlineExp_fitting_groupSize_supp1_78 <- readRDS("./rds/mcmc_result/onlineExp_fitting_groupSize_supp1_78.rds")
isPositiveCopier_fitting <- readRDS("./rds/mcmc_result/isPositiveCopier_fitting.rds")
alpha_fitting <- readRDS("./rds/mcmc_result/alpha_fitting.rds")
average_copy_rate_fitting <- readRDS("./rds/mcmc_result/average_copy_rate_fitting.rds")
average_copy_rate_fitting_posOnly <- readRDS("./rds/mcmc_result/average_copy_rate_fitting_posOnly.rds")
average_invTemp_fitting <- readRDS("./rds/mcmc_result/average_invTemp_fitting.rds")
conformity_exponent_fitting <- readRDS("./rds/mcmc_result/conformity_exponent_fitting.rds")
conformity_exponent_fitting_posOnly <- readRDS("./rds/mcmc_result/conformity_exponent_fitting_posOnly.rds")


## taskDifficulty3 means the three different uncertainty conditions
allBehaviouralData2$taskDifficulty3 = allBehaviouralData2$taskDifficulty
allBehaviouralData2$taskDifficulty3[which(allBehaviouralData2$taskDifficulty=='25%')] <- 'Low Uncertainty'
allBehaviouralData2$taskDifficulty3[which(allBehaviouralData2$taskDifficulty=='50%')] <- 'Moderate Uncertainty'
allBehaviouralData2$taskDifficulty3[which(allBehaviouralData2$taskDifficulty=='78%')] <- 'High Uncertainty'
allBehaviouralData2$taskDifficulty3 = factor(allBehaviouralData2$taskDifficulty3, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))

performanceSummary4$taskDifficulty3 = NA
performanceSummary4$taskDifficulty3[which(performanceSummary4$taskDifficulty=='25%')] <- 'Low Uncertainty'
performanceSummary4$taskDifficulty3[which(performanceSummary4$taskDifficulty=='50%')] <- 'Moderate Uncertainty'
performanceSummary4$taskDifficulty3[which(performanceSummary4$taskDifficulty=='78%')] <- 'High Uncertainty'
performanceSummary4$taskDifficulty3 = factor(performanceSummary4$taskDifficulty3, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))

posthoc_sim_eachGroup_aveOptim$taskDifficulty3 = factor(posthoc_sim_eachGroup_aveOptim$taskDifficulty3, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))
posthoc_sim_summary_indivSmallLarge$taskDifficulty3 = factor(posthoc_sim_summary_indivSmallLarge$taskDifficulty3, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))
posthoc_sim_eachGroup_aveOptim$indiv_small_large = factor(posthoc_sim_eachGroup_aveOptim$indiv_small_large, levels=c('Large','Small','Single'))

## where is median of the group sizes?
median_low = 9
median_moderate = 6
median_high = 11
## sizeCategoryBinary
allBehaviouralData2$sizeCategoryBinary = NA
allBehaviouralData2$sizeCategoryBinary[which(allBehaviouralData2$taskDifficulty=='25%'&allBehaviouralData2$sizeCategory<=median_low)] = 'small'
allBehaviouralData2$sizeCategoryBinary[which(allBehaviouralData2$taskDifficulty=='25%'&allBehaviouralData2$sizeCategory >median_low)] = 'large'
allBehaviouralData2$sizeCategoryBinary[which(allBehaviouralData2$taskDifficulty=='50%'&allBehaviouralData2$sizeCategory<=median_moderate)] = 'small'
allBehaviouralData2$sizeCategoryBinary[which(allBehaviouralData2$taskDifficulty=='50%'&allBehaviouralData2$sizeCategory >median_moderate)] = 'large'
allBehaviouralData2$sizeCategoryBinary[which(allBehaviouralData2$taskDifficulty=='78%'&allBehaviouralData2$sizeCategory<=median_high)] = 'small'
allBehaviouralData2$sizeCategoryBinary[which(allBehaviouralData2$taskDifficulty=='78%'&allBehaviouralData2$sizeCategory >median_high)] = 'large'
allBehaviouralData2$sizeCategoryBinary[which(allBehaviouralData2$sizeCategory==1)] = 'single'
allBehaviouralData2$sizeCategoryBinary = factor(allBehaviouralData2$sizeCategoryBinary, levels=c('single', 'small', 'large'))

allBehaviouralData2$room_indiv_group = allBehaviouralData2$room
allBehaviouralData2$room_indiv_group[which(allBehaviouralData2$sizeCategory==1)] = 'Single-condition'
allBehaviouralData2$sizeCategory_string = apply(rbind(rep('n = ', nrow(allBehaviouralData2)),allBehaviouralData2$sizeCategory), 2, str_c, collapse="")
allBehaviouralData2$sizeCategory_string = factor(allBehaviouralData2$sizeCategory_string, levels=c("n = 1","n = 2","n = 3","n = 4","n = 5","n = 6","n = 7","n = 8","n = 9","n = 10","n = 11","n = 12","n = 13","n = 14","n = 15","n = 16","n = 17","n = 18","n = 19","n = 20","n = 21","n = 22","n = 23","n = 24","n = 25","n = 26","n = 27"))

#######################################################################
##
## Figure 1 and Figure 2
## Simulation Results
##
#######################################################################
low_sigma = 0.3
high_sigma = 0.9
low_theta = 1
high_theta = 3
lifetime_simulation = 90
lifetime_experiment = 70

## To get the tidy-data frame for the time series data
UNC_Moderate_expParameters_timeSeries = subset(UNC_Moderate_expParameters_timeSeries0, theta==low_theta|theta==high_theta)
UNC_Moderate_expParameters_timeSeries2 = data.frame(
  round = rep(1:lifetime_simulation, 2*2*3),
  groupSize = c(rep(3, lifetime_simulation*2*2), rep(10, lifetime_simulation*2*2), rep(30, lifetime_simulation*2*2)),
  groupSize_factor = c(rep('n = 3', lifetime_simulation*2*2), rep('n = 10', lifetime_simulation*2*2), rep('n = 30', lifetime_simulation*2*2)),
  theta = rep(c(rep(low_theta, lifetime_simulation*2), rep(high_theta,lifetime_simulation*2)), 3),
  #theta = rep(c(rep(1, lifetime_simulation*2), rep(3,lifetime_simulation*2)), 3),
  lambda = rep(rep(c(rep(low_sigma,lifetime_simulation), rep(high_sigma,lifetime_simulation)), 2), 3)
)
UNC_Moderate_expParameters_timeSeries2$groupSize_factor = factor(UNC_Moderate_expParameters_timeSeries2$groupSize_factor, levels=c('n = 3','n = 10','n = 30'))
UNC_Moderate_expParameters_timeSeries2$optimChoiceRate = NA
for(i in 1:nrow(UNC_Moderate_expParameters_timeSeries)) {
  g = as.numeric(as.character(UNC_Moderate_expParameters_timeSeries$groupSize[i]))
  the = as.numeric(as.character(UNC_Moderate_expParameters_timeSeries$theta[i]))
  la = as.numeric(as.character(UNC_Moderate_expParameters_timeSeries$lambda[i]))
  location = which(UNC_Moderate_expParameters_timeSeries2$groupSize==g&UNC_Moderate_expParameters_timeSeries2$theta==the&UNC_Moderate_expParameters_timeSeries2$lambda==la)
  if(length(location)>0) {
  UNC_Moderate_expParameters_timeSeries2$optimChoiceRate[location] = UNC_Moderate_expParameters_timeSeries[i,4:(lifetime_simulation+3)]
  }
}
UNC_Moderate_expParameters_timeSeries2$optimChoiceRate = as.numeric(as.character(UNC_Moderate_expParameters_timeSeries2$optimChoiceRate))

for(t in 1:lifetime_simulation) {
  UNC_Moderate_expParameters_timeSeries2$optimChoiceRate[which(UNC_Moderate_expParameters_timeSeries2$groupSize==1&UNC_Moderate_expParameters_timeSeries2$round==t)] = mean(UNC_Moderate_expParameters_timeSeries2$optimChoiceRate[which(UNC_Moderate_expParameters_timeSeries2$groupSize==1&UNC_Moderate_expParameters_timeSeries2$round==t)])
}

Figure1_a =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, lambda==low_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize_factor, colour=groupSize_factor, linetype=groupSize_factor),size=1.5)+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.3\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_manual(name = '', values = c("black", "orange", "red"))+
  scale_linetype_manual(name = '', values = c("dotted", "dashed", "solid"))+
  #scale_color_distiller(name = '', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme()+
  theme(legend.position = c(0.00, 0.23),
    legend.key.width = unit(15, "mm"))
Figure1_b =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==3 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==10 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==30 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.3\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme()
Figure1_e =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==3 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==10 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==30 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme()
Figure1_f =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==3 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==10 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==30 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme()

Figure1 = plot_grid(Figure1_a, Figure1_b, Figure1_e, Figure1_f, labels = c('A','B','C','D'), ncol = 2, align = 'v', label_size=13, label_fontfamily='Helvetica',label_fontface='bold',label_x=0.25,label_y=0.85,hjust=-0.5,vjust=-0.5)
ggsave(file = "Figure1_rev.pdf", plot = Figure1, dpi = 600, width = 7.5, height = 5)
#ggsave(file = "Figure1.tiff", plot = Figure1, dpi = 600, width = 9, height = 6)

## GGjoy plot
Figure2 = ggplot()+
  geom_joy(data= UNC_Moderate_expParameters_eachGroup, aes(x = averageFirst, y = as.factor(lambda), linetype=as.factor(groupSize), colour=as.factor(groupSize),fill=as.factor(groupSize),alpha=as.factor(groupSize)), stat = "binline", bins = 40, scale = 0.98, draw_baseline = FALSE)+
  scale_color_manual(name = 'Group size:', values = c("black", "orange", "red"), labels=c("n = 3;", "n = 10;", "n = 30"))+
  scale_linetype_manual(name = 'Group size:', values = c("dotted", "dashed", "solid"), labels=c("n = 3;", "n = 10;", "n = 30"))+
  scale_fill_manual(name = 'Group size:', values = c(NA, "orange", "red"), labels=c("n = 3;", "n = 10;", "n = 30"))+
  scale_alpha_manual(name = 'Group size:', values = c(3/3, 2/3, 1/3), labels=c("n = 3;", "n = 10;", "n = 30"))+
  facet_grid(.~theta, labeller = label_bquote(cols = bar(italic(theta)) == .(theta)))+
  labs(x='Groups\' average performance', y=expression(paste('Social learning weight ',bar(sigma),sep="")),title='')+
  xlim(c(0,1))+
  theme(
    legend.position = 'top',
    axis.text.x = element_text(angle=90),
    axis.title.y = element_text(angle=90))+
  #scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  myTheme_legend()

ggsave(file = "Figure2.pdf", plot = Figure2, dpi = 600, width = 6.5, height = 5)
#ggsave(file = "Figure2.tiff", plot = Figure2, dpi = 600, width = 9, height = 6)




#######################################################################
##
## Figure 3
## Experimental Results --
## Time evolutions and distributions of decision performance for each condition.
##
#######################################################################
## preparation
allBehaviouralData2 = allBehaviouralData2 %>% dplyr::mutate(T2 = allBehaviouralData2$round)
allBehaviouralData2$indiv_group = 0
allBehaviouralData2$indiv_group[which(allBehaviouralData2$sizeCategory>1)] = 1
allBehaviouralData_25 = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%')
allBehaviouralData_50 = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%')
allBehaviouralData_78 = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%')

data_onlineExp_fitting_groupSize_25 <- rstan::extract(onlineExp_fitting_groupSize_25)
fittedValues_onlineExp_fitting_groupSize_25 =
  as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$mu[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$di[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_indiv[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_group_small[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_group_large[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_group_veryLarge[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_group[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
## rename
colnames(fittedValues_onlineExp_fitting_groupSize_25) = c("mu_p2.5", "mu_p25", "mu_p50", "mu_p75", "mu_p97.5","di_p2.5", "di_p25", "di_p50", "di_p75", "di_p97.5","p_indiv_p2.5", "p_indiv_p25", "p_indiv_p50", "p_indiv_p75", "p_indiv_p97.5","p_group_small_p2.5", "p_group_small_p25", "p_group_small_p50", "p_group_small_p75", "p_group_small_p97.5","p_group_large_p2.5", "p_group_large_p25", "p_group_large_p50", "p_group_large_p75", "p_group_large_p97.5","p_group_veryLarge_p2.5", "p_group_veryLarge_p25", "p_group_veryLarge_p50", "p_group_veryLarge_p75", "p_group_veryLarge_p97.5","p_group_p2.5", "p_group_p25", "p_group_p50", "p_group_p75", "p_group_p97.5")
fittedValues_onlineExp_fitting_groupSize_25$taskDifficulty = "25%"

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

data_onlineExp_fitting_groupSize_78 <- rstan::extract(onlineExp_fitting_groupSize_78)
fittedValues_onlineExp_fitting_groupSize_78 =
  as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$mu[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$di[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_indiv[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_group_small[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_group_large[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_group_veryLarge[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_group[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
## rename
colnames(fittedValues_onlineExp_fitting_groupSize_78) = c("mu_p2.5", "mu_p25", "mu_p50", "mu_p75", "mu_p97.5","di_p2.5", "di_p25", "di_p50", "di_p75", "di_p97.5","p_indiv_p2.5", "p_indiv_p25", "p_indiv_p50", "p_indiv_p75", "p_indiv_p97.5","p_group_small_p2.5", "p_group_small_p25", "p_group_small_p50", "p_group_small_p75", "p_group_small_p97.5","p_group_large_p2.5", "p_group_large_p25", "p_group_large_p50", "p_group_large_p75", "p_group_large_p97.5","p_group_veryLarge_p2.5", "p_group_veryLarge_p25", "p_group_veryLarge_p50", "p_group_veryLarge_p75", "p_group_veryLarge_p97.5","p_group_p2.5", "p_group_p25", "p_group_p50", "p_group_p75", "p_group_p97.5")
fittedValues_onlineExp_fitting_groupSize_78$taskDifficulty = "78%"

####
## Plot
####

fittedValues_onlineExp_fitting_groupSize_25$taskDifficulty3 = "Low Uncertainty"
fittedValues_onlineExp_fitting_groupSize_50$taskDifficulty3 = "Moderate Uncertainty"
fittedValues_onlineExp_fitting_groupSize_78$taskDifficulty3 = "High Uncertainty"

fittedValues_onlineExp_fitting_groupSize = rbind(fittedValues_onlineExp_fitting_groupSize_25, fittedValues_onlineExp_fitting_groupSize_50) %>% rbind(fittedValues_onlineExp_fitting_groupSize_78)
fittedValues_onlineExp_fitting_groupSize$round = rep(1:70, 3)
fittedValues_onlineExp_fitting_groupSize$taskDifficulty3 = factor(fittedValues_onlineExp_fitting_groupSize$taskDifficulty3, levels=c("Low Uncertainty","Moderate Uncertainty","High Uncertainty"))

fittedValues_onlineExp_fitting_groupSize$significant_indiv_vs_group = 0
fittedValues_onlineExp_fitting_groupSize$significant_indiv_vs_group[which(fittedValues_onlineExp_fitting_groupSize$di_p2.5>0)] = 1

group_level_summary_small_25 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=9&taskDifficulty=='25%')$optimalChoice, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=9&taskDifficulty=='25%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=9&taskDifficulty=='25%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_large_25 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>9&taskDifficulty=='25%')$optimalChoice, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>9&taskDifficulty=='25%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory>9&taskDifficulty=='25%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_small_50 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=6&taskDifficulty=='50%')$optimalChoice, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=6&taskDifficulty=='50%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=6&taskDifficulty=='50%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_large_50 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>6&taskDifficulty=='50%')$optimalChoice, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>6&taskDifficulty=='50%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory>6&taskDifficulty=='50%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_small_78 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=11&taskDifficulty=='78%')$optimalChoice, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=11&taskDifficulty=='78%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=11&taskDifficulty=='78%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_large_78 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>11&taskDifficulty=='78%')$optimalChoice, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>11&taskDifficulty=='78%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory>11&taskDifficulty=='78%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)

group_level_summary = data.frame(
  round = rep(1:70, 3),
  taskDifficulty3 = factor(c(rep('Low Uncertainty', 70), rep('Moderate Uncertainty',70), rep('High Uncertainty',70)),levels=c("Low Uncertainty","Moderate Uncertainty","High Uncertainty")),
  small = c(group_level_summary_small_25, group_level_summary_small_50, group_level_summary_small_78),
  large = c(group_level_summary_large_25, group_level_summary_large_50, group_level_summary_large_78)
  )

## plot
timeseries_optimrate_analisys_groupSize_plot = ggplot() +
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_indiv_p2.5, ymax=p_indiv_p97.5),fill='black', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_indiv_p25, ymax=p_indiv_p75),fill='black', alpha=3/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_group_small_p2.5, ymax=p_group_small_p97.5),fill='orange', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_group_small_p25, ymax=p_group_small_p75),fill='orange', alpha=3/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_group_large_p2.5, ymax=p_group_large_p97.5),fill='red', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_group_large_p25, ymax=p_group_large_p75),fill='red', alpha=3/6)+
  stat_summary(data=subset(allBehaviouralData2, indiv_group==0), mapping=aes(x=round,optimalChoice), fun.y = mean, geom="line",colour='black',size=1,alpha=3/4)+
  geom_line(data=group_level_summary, mapping=aes(round, small), colour='orange', size=1, alpha=3/4)+
  geom_line(data=group_level_summary, mapping=aes(round, large), colour='red', size=1, alpha=3/4)+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, y=p_indiv_p50),colour='black',linetype='dashed')+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, y=p_group_small_p50),colour='orange',linetype='dashed')+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, y=p_group_large_p50),colour='red',linetype='dashed')+
  facet_grid(.~taskDifficulty3)+
  scale_y_continuous(limits = c(-0.05,1.05), breaks=c(0.0,0.25,0.50,0.75,1.0)) +
  labs(x='Rounds', y='Proportion \nchoosing\n the best option')+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  geom_hline(yintercept=1/3, linetype="dashed")+
  myTheme_legend()+
  NULL
#ggsave(file = "timeseries_optimrate_analisys_groupSize_plot.pdf", plot = timeseries_optimrate_analisys_groupSize_plot, dpi = 600, width = 9, height = 3)

timeseries_optimrate_analisys_indiv_vs_group_plot = ggplot() +
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_indiv_p2.5, ymax=p_indiv_p97.5),fill='black', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_indiv_p25, ymax=p_indiv_p75),fill='black', alpha=3/6)+
  #geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_group_small_p2.5, ymax=p_group_small_p97.5),fill='orange', alpha=1/6)+
  #geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_group_small_p25, ymax=p_group_small_p75),fill='orange', alpha=3/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_group_p2.5, ymax=p_group_p97.5),fill='red', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=p_group_p25, ymax=p_group_p75),fill='red', alpha=3/6)+
  stat_summary(data=subset(allBehaviouralData2, indiv_group==0), mapping=aes(x=round,optimalChoice), fun.y = mean, geom="line",colour='black',size=1,alpha=3/4)+
  stat_summary(data=subset(allBehaviouralData2, indiv_group==1), mapping=aes(x=round,optimalChoice), fun.y = mean, geom="line",colour='red',size=1,alpha=3/4)+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, y=p_indiv_p50),colour='black',linetype='dashed')+
  #geom_line(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, y=p_group_small_p50),colour='orange',linetype='dashed')+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, y=p_group_p50),colour='red',linetype='dashed')+
  facet_grid(.~taskDifficulty3)+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  geom_hline(yintercept=1/3, linetype="dashed")+
  myTheme_legend()+
  NULL
#ggsave(file = "timeseries_optimrate_analisys_indiv_vs_group_plot.pdf", plot = timeseries_optimrate_analisys_indiv_vs_group_plot, dpi = 600, width = 9, height = 3)

difference_btwn_indiv_vs_group_plot = ggplot() +
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=di_p2.5, ymax=di_p97.5),fill='black', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, ymin=di_p25, ymax=di_p75),fill='black', alpha=3/6)+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize, mapping=aes(x=round, y=di_p50),colour='black')+
  facet_grid(.~taskDifficulty3)+
  labs(x='Rounds', y=expression(paste('Differences between\nsolitaries and groups ',xi[t],sep="")))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  geom_hline(yintercept=0, linetype="dashed")+
  myTheme_legend()+
  NULL

## Posthoc simulation result
posthoc_simulation_plot = ggplot() +
  geom_ribbon(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Individual'), mapping=aes(round, ymin=X25., ymax=X75.), fill='black', alpha = 1/6) +
  geom_line(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Individual'), mapping=aes(round, y=mean), colour='black', size=1, linetype='dashed') +
  geom_ribbon(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Small'), mapping=aes(round, ymin=X25., ymax=X75.), fill='orange', alpha = 1/6) +
  geom_line(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Small'), mapping=aes(round, y=mean), colour='orange', size=1, linetype='twodash') +
  geom_ribbon(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Large'), mapping=aes(round, ymin=X25., ymax=X75.), fill='red', alpha = 1/6) +
  geom_line(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Large'), mapping=aes(round, y=mean), colour='red', size=1, linetype='solid') +
  facet_grid(. ~ taskDifficulty3) +
  scale_y_continuous(limits = c(-0.05,1.05), breaks=c(0.0,0.25,0.50,0.75,1.0)) +
  scale_x_continuous(breaks=c(0,20,40,70,90)) +
  geom_vline(xintercept=70)+
  labs(x='Rounds', y='Proportion \nchoosing\n the best option')+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  geom_hline(yintercept=1/3, linetype="dashed")+
  myTheme_legend()+
  #theme(legend.position='none',
  #axis.text.y = element_text(size=11, family="Helvetica" ,colour='black'),
  #strip.text.y = element_text(angle = 0))+
  NULL

Figure_rev3_top = plot_grid(timeseries_optimrate_analisys_groupSize_plot, difference_btwn_indiv_vs_group_plot, labels = c('',''), ncol = 1, align = 'v')
#ggsave(file = "mean_optim_prob.pdf", plot = Figure_rev3_top, dpi = 600, width = 9, height = 5)

posthoc_sim_eachGroup_aveOptim_plot = ggplot() +
  geom_density(data = posthoc_sim_eachGroup_aveOptim, mapping=aes(x=meanOptimProb, y=..scaled.., colour=indiv_small_large,fill=indiv_small_large,linetype=indiv_small_large),alpha=1/3) +
  #geom_density(data = posthoc_sim_eachGroup_aveOptim, mapping=aes(x=meanOptimProb, colour=indiv_small_large,fill=indiv_small_large,linetype=indiv_small_large), stat='density',alpha=1/3) +
  labs(y='Density\n(scaled)',x='Groups\' average proportion\nchoosing the best option')+
  scale_colour_manual(name = 'Group size', values=c("Single"='black','Small'='orange', 'Large'='red'))+
  scale_fill_manual(name = 'Group size', values=c("Single"='black','Small'='orange', 'Large'='red'))+
  scale_linetype_manual(name = 'Group size', values=c("Single"='dashed','Small'='twodash', 'Large'='solid'))+
  scale_y_continuous(limits = c(-0.05,1.05), breaks=c(0.0,0.5,1.0)) +
  scale_x_continuous(limits = c(-0.05,1.05), breaks=c(0.0,0.5,1.0)) +
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  myTheme_legend()+
  theme(strip.text = element_text(size=9, family="Helvetica", lineheight=1)) +
  facet_grid(environment~taskDifficulty3, labeller = labeller(environment = label_wrap_gen(9)))+
  NULL

Figure_posthocSim = plot_grid(posthoc_simulation_plot, posthoc_sim_eachGroup_aveOptim_plot, labels = c('',''), ncol = 1, align = 'h')
#ggsave(file = "post_hoc_simulation.pdf", plot = Figure_posthocSim, dpi = 600, width = 9, height = 5)

## Figure 3 in the revised manuscript
Figure3 = plot_grid(Figure_rev3_top, Figure_posthocSim, labels=c(''),ncol=1,align='h')
ggsave(file = "Figure3.pdf", plot = Figure3, dpi = 600, width = 7.5, height = 7.5)







#######################################################################
##
## Figure 5
## Experimental Results --
## Change in fitted values (i.e. median of the Bayesian posterior distribution) of
## the social learning weight 𝜎_{𝑖,𝑡} with time for each individual, for each level of task uncertainty.
##
#######################################################################
# conservative categories (rand, neg, pos) - 25-75
copyingProbReduction_data_fullModelOnly$frequency_dependence_4 = 'Random-copying'
copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$theta_75 < 0)] = 'neg-freq-dep'
copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$theta_25 > 0)] = 'pos-freq-dep'
#copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$best_model == 'AL' | copyingProbReduction_data_fullModelOnly$best_model == 'AL_ann')] = 'AL'
#copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$best_model == 'Random')] = 'Random'
copyingProbReduction_data_fullModelOnly$frequency_dependence_4 = factor(copyingProbReduction_data_fullModelOnly$frequency_dependence_4, levels=c("Random choice","AL","Random-copying",'neg-freq-dep', "pos-freq-dep"))

## Task difficulty 2
copyingProbReduction_data_fullModelOnly$taskDifficulty2 = NA
copyingProbReduction_data_fullModelOnly$taskDifficulty2[which(copyingProbReduction_data_fullModelOnly$taskDifficulty=='25%')] <- 'Low Uncertainty'
copyingProbReduction_data_fullModelOnly$taskDifficulty2[which(copyingProbReduction_data_fullModelOnly$taskDifficulty=='50%')] <- 'Moderate Uncertainty'
copyingProbReduction_data_fullModelOnly$taskDifficulty2[which(copyingProbReduction_data_fullModelOnly$taskDifficulty=='78%')] <- 'High Uncertainty'
copyingProbReduction_data_fullModelOnly$taskDifficulty2 = factor(copyingProbReduction_data_fullModelOnly$taskDifficulty2, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))


# Using 50% CI for categorisation
# Three categories: Random-copying, neg-freq-dep, pos-freq-dep
# Presenting only the positive frequency-dependent copiers

Figure5 = ggplot() +
                    geom_line(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'&frequency_dependence_4=='pos-freq-dep'), mapping=aes(round,copyingProb,group=amazonID,colour=theta),alpha=3/4)+
                    stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'&frequency_dependence_4=='pos-freq-dep'),aes(round,copyingProb),fun.y = "median",geom="line",linetype = "dashed", colour = "black", size = 1)+
                    facet_grid(.~taskDifficulty2)+
                    labs(x='Round', y=expression(paste('Social learning weight ',sigma[i][t],sep="")), title='') +
                    ylim(c(0,1))+
                    xlim(c(2,70))+
                    #scale_colour_distiller(name=expression(paste('Conformity\nexponent',theta[i],sep="")), palette = "RdYlBu", direction = -1)+
                    scale_color_viridis(name=expression(paste('Conformity\nexponent',theta[i],sep="")), option="magma", direction = -1)+
                    myTheme_legend()+
                    #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
                    NULL
ggsave(file = "Figure5.pdf", plot = Figure5, dpi = 600, width = 7.5, height = 2.5)



#######################################################################
##
## Figure 4
## Experimental Results --
## Model fitting for the three different task’s uncertain conditions
## (the Low-, Moderate- and High-uncertainty) and the different group size.
##
#######################################################################
## STATISTICAL ANALYSIS - ALL STRATEGIES MERGED (Analyses only on the Pos-freq-dep follows below)
performanceSummary4_stanData = performanceSummary4
performanceSummary4_stanData$groupSize_standerized2 = (performanceSummary4$groupSize-mean(subset(performanceSummary4,environment=='environment 1')$groupSize))/sd(subset(performanceSummary4,environment=='environment 1')$groupSize)
performanceSummary4_stanData$groupSize_centralized2 = (performanceSummary4$groupSize-mean(subset(performanceSummary4,environment=='environment 1')$groupSize))
performanceSummary4_stanData$age_standardized = (performanceSummary4$age-mean(subset(performanceSummary4,environment=='environment 1')$age,na.rm=TRUE))/sd(subset(performanceSummary4,environment=='environment 1')$age,na.rm=TRUE)
performanceSummary4_stanData$female = performanceSummary4$male_female_other - 1
performanceSummary4_stanData$female[which(performanceSummary4$female==2)] <- NA

performanceSummary4_stanData = subset(performanceSummary4_stanData, (age_standardized > -100) & (female==0 | female==1))

  ### Setting a group level data
performanceSummary4_stanData_group = data.frame(
  groupID = names(table(performanceSummary4_stanData$groupID)))
performanceSummary4_stanData_group$environment = 'e1+e2'
performanceSummary4_stanData_group$taskDifficulty = NA
performanceSummary4_stanData_group$average_entropy = NA
performanceSummary4_stanData_group$groupSize = NA
#performanceSummary4_stanData_group$totalPositiveCopier = NA
for(g in names(table(performanceSummary4_stanData$groupID))){
  performanceSummary4_stanData_group$taskDifficulty[which(performanceSummary4_stanData_group$groupID==g)] = performanceSummary4_stanData$taskDifficulty[which(performanceSummary4_stanData$groupID==g)][1]
  performanceSummary4_stanData_group$average_entropy[which(performanceSummary4_stanData_group$groupID==g)] = performanceSummary4_stanData$average_entropy[which(performanceSummary4_stanData$groupID==g)][1]
  performanceSummary4_stanData_group$groupSize[which(performanceSummary4_stanData_group$groupID==g)] = performanceSummary4_stanData$groupSize[which(performanceSummary4_stanData$groupID==g)][1]
  #performanceSummary4_stanData_group$totalPositiveCopier[which(performanceSummary4_stanData_group$groupID==g)] = performanceSummary4_stanData$totalPositiveCopier[which(performanceSummary4_stanData$groupID==g)][1]
}
  ## Corlelations between social learning weight and group sizes
performanceSummary4_stanData_groupOnly = subset(performanceSummary4_stanData_group, groupSize>1)
performanceSummary4_stanData_groupOnly = subset(performanceSummary4_stanData_groupOnly, groupID != 'subRoom1_1_6') # no data
performanceSummary4_stanData_groupOnly$groupID_groupOnly = 1:nrow(performanceSummary4_stanData_groupOnly)
performanceSummary4_stanData_groupOnly$groupID = as.character(performanceSummary4_stanData_groupOnly$groupID)
performanceSummary4_stanData$groupID_groupOnly = NA
for(g in names(table(performanceSummary4_stanData_groupOnly$groupID))) {
  performanceSummary4_stanData$groupID_groupOnly[which(performanceSummary4_stanData$groupID==g)] <- performanceSummary4_stanData_groupOnly$groupID_groupOnly[which(performanceSummary4_stanData_groupOnly$groupID==g)]
}
performanceSummary4_stanData$taskDifficulty_num = 0
performanceSummary4_stanData$taskDifficulty_num[which(performanceSummary4_stanData$taskDifficulty=='50%')] = 0.5
performanceSummary4_stanData$taskDifficulty_num[which(performanceSummary4_stanData$taskDifficulty=='78%')] = 1.0

performanceSummary4_stanData$isPositiveCopier = 0
performanceSummary4_stanData$isPositiveCopier[which(performanceSummary4_stanData$frequency_dependence_4=='Pos-freq-dep')] = 1



  ## average_copy_rate vs group size vs task difficulty
parameters_stan_data = list(
  All = nrow(subset(subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition'))),
  N_group = length(table((subset(subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition'))$groupID))),
  groupID = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$groupID_groupOnly,
  taskDifficulty = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$taskDifficulty_num,
  groupSize_standerized = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$groupSize_standerized2,
  age_standardized = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$age_standardized,
  female = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$female,
  average_copy_rate = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$average_copy_rate,
  test_group_size = seq(-1.3,2.3,len=200),
  test_size = 200
  )
# social learning strategies using four categories
performanceSummary4$frequency_dependence_four = performanceSummary4$frequency_dependence_4
performanceSummary4$frequency_dependence_four[which(performanceSummary4$frequency_dependence_4=='Pos-freq-dep')] = 'Weak-pos-freq-dep'
performanceSummary4$frequency_dependence_four[which(performanceSummary4$frequency_dependence_4=='Pos-freq-dep'&performanceSummary4$theta>1)] = 'Strong-pos-freq-dep'

# 25% small uncertainty
strategies_proportion_25 = data.frame(groupSize = as.numeric(names(table(subset(performanceSummary4, taskDifficulty=="25%")$groupSize)))[-1])
strategies_proportion_25$proportion_random_copy = NA
strategies_proportion_25$proportion_neg_freq_dep = NA
strategies_proportion_25$proportion_pos_freq_dep = NA
strategies_proportion_25$proportion_weak_pos_freq_dep = NA
strategies_proportion_25$proportion_strong_pos_freq_dep = NA

for(i in strategies_proportion_25$groupSize){
    # proportion
    strategies_proportion_25$proportion_random_copy[which(strategies_proportion_25$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_4=="Random-copying"&taskDifficulty=='25%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='25%'&groupSize==i))
    strategies_proportion_25$proportion_neg_freq_dep[which(strategies_proportion_25$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_4=="Neg-freq-dep"&taskDifficulty=='25%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='25%'&groupSize==i))
    strategies_proportion_25$proportion_pos_freq_dep[which(strategies_proportion_25$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_4=="Pos-freq-dep"&taskDifficulty=='25%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='25%'&groupSize==i))
    strategies_proportion_25$proportion_weak_pos_freq_dep[which(strategies_proportion_25$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_four=="Weak-pos-freq-dep"&taskDifficulty=='25%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='25%'&groupSize==i))
    strategies_proportion_25$proportion_strong_pos_freq_dep[which(strategies_proportion_25$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_four=="Strong-pos-freq-dep"&taskDifficulty=='25%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='25%'&groupSize==i))
}

# 50% moderate uncertainty
strategies_proportion_50 = data.frame(groupSize = as.numeric(names(table(subset(performanceSummary4, taskDifficulty=="50%")$groupSize)))[-1])
strategies_proportion_50$proportion_random_copy = NA
strategies_proportion_50$proportion_neg_freq_dep = NA
strategies_proportion_50$proportion_pos_freq_dep = NA
strategies_proportion_50$proportion_weak_pos_freq_dep = NA
strategies_proportion_50$proportion_strong_pos_freq_dep = NA

for(i in strategies_proportion_50$groupSize){
    # proportion
    strategies_proportion_50$proportion_random_copy[which(strategies_proportion_50$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_4=="Random-copying"&taskDifficulty=='50%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='50%'&groupSize==i))
    strategies_proportion_50$proportion_neg_freq_dep[which(strategies_proportion_50$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_4=="Neg-freq-dep"&taskDifficulty=='50%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='50%'&groupSize==i))
    strategies_proportion_50$proportion_pos_freq_dep[which(strategies_proportion_50$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_4=="Pos-freq-dep"&taskDifficulty=='50%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='50%'&groupSize==i))
    strategies_proportion_50$proportion_weak_pos_freq_dep[which(strategies_proportion_50$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_four=="Weak-pos-freq-dep"&taskDifficulty=='50%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='50%'&groupSize==i))
    strategies_proportion_50$proportion_strong_pos_freq_dep[which(strategies_proportion_50$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_four=="Strong-pos-freq-dep"&taskDifficulty=='50%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='50%'&groupSize==i))
}

# 78% high uncertainty
strategies_proportion_78 = data.frame(groupSize = as.numeric(names(table(subset(performanceSummary4, taskDifficulty=="78%")$groupSize)))[-1])
strategies_proportion_78$proportion_random_copy = NA
strategies_proportion_78$proportion_neg_freq_dep = NA
strategies_proportion_78$proportion_pos_freq_dep = NA
strategies_proportion_78$proportion_weak_pos_freq_dep = NA
strategies_proportion_78$proportion_strong_pos_freq_dep = NA

for(i in strategies_proportion_78$groupSize){
    # proportion
    strategies_proportion_78$proportion_random_copy[which(strategies_proportion_78$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_4=="Random-copying"&taskDifficulty=='78%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='78%'&groupSize==i))
    strategies_proportion_78$proportion_neg_freq_dep[which(strategies_proportion_78$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_4=="Neg-freq-dep"&taskDifficulty=='78%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='78%'&groupSize==i))
    strategies_proportion_78$proportion_pos_freq_dep[which(strategies_proportion_78$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_4=="Pos-freq-dep"&taskDifficulty=='78%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='78%'&groupSize==i))
    strategies_proportion_78$proportion_weak_pos_freq_dep[which(strategies_proportion_78$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_four=="Weak-pos-freq-dep"&taskDifficulty=='78%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='78%'&groupSize==i))
    strategies_proportion_78$proportion_strong_pos_freq_dep[which(strategies_proportion_78$groupSize==i)] <- nrow(subset(performanceSummary4,environment=='environment 1'&frequency_dependence_four=="Strong-pos-freq-dep"&taskDifficulty=='78%'&groupSize==i))/nrow(subset(performanceSummary4,environment=='environment 1'&taskDifficulty=='78%'&groupSize==i))
}

isPositiveCopier_mcmc_result <- rstan::extract(isPositiveCopier_fitting)
isPositiveCopier_pred_data = cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(isPositiveCopier_mcmc_result$isPositiveCopier_pred_25[,1:ncol(isPositiveCopier_mcmc_result$isPositiveCopier_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
isPositiveCopier_pred_data = rbind(isPositiveCopier_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(isPositiveCopier_mcmc_result$isPositiveCopier_pred_50[,1:ncol(isPositiveCopier_mcmc_result$isPositiveCopier_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
isPositiveCopier_pred_data = rbind(isPositiveCopier_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(isPositiveCopier_mcmc_result$isPositiveCopier_pred_78[,1:ncol(isPositiveCopier_mcmc_result$isPositiveCopier_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
isPositiveCopier_pred_data$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(isPositiveCopier_pred_data) <- c("groupSize_standerized2",
              "isPositiveCopier_p2.5", "isPositiveCopier_p25", "isPositiveCopier_p50", "isPositiveCopier_p75", "isPositiveCopier_p97.5","taskDifficulty")

lsfit(performanceSummary4$groupSize_standerized2, performanceSummary4$groupSize)$coefficients
#Intercept         X
#10.421192  7.808047

isPositiveCopier_pred_data$groupSize = isPositiveCopier_pred_data$groupSize_standerized2*7.808047 + 10.421192

## Figure 4 top row
GS_proportionSL_25 = ggplot() +
  geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_25, aes(groupSize, proportion_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_25,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_25,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  labs(x='Group size', y='Proportion of \nstrategies', title='Low Uncertainty') +
  myTheme_legend()+ylim(c(NA,1))+xlim(c(2, 28))+
  theme(legend.position = c(0.05, 0.5))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL
GS_proportionSL_50 = ggplot() +
  geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_50, aes(groupSize, proportion_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_50,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_50,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  #labs(x='Group size', y='Proportion of \nstrategies', title='Moderate Uncertainty') +
  labs(x='Group size', y='', title='Moderate Uncertainty') +
  myTheme_legend()+ylim(c(NA,1))+xlim(c(2, 28))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL
GS_proportionSL_78 = ggplot() +
  geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=2.2), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=2.2), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_78, aes(groupSize, proportion_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_78,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_78,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  #labs(x='Group size', y='Proportion of \nstrategies', title='High Uncertainty') +
  labs(x='Group size', y='', title='High Uncertainty') +
  myTheme_legend()+ylim(c(NA,1))+xlim(c(2, 28))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL

# Summary plot
#Figure4_top = plot_grid(GS_proportionSL_25, GS_proportionSL_50, GS_proportionSL_78, labels = c("","",""), ncol = 3, align = 'v')

average_copy_rate_mcmc_result <- rstan::extract(average_copy_rate_fitting)

average_copy_rate_pred_data = cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_copy_rate_mcmc_result$average_copy_rate_pred_25[,1:ncol(average_copy_rate_mcmc_result$average_copy_rate_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
average_copy_rate_pred_data = rbind(average_copy_rate_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_copy_rate_mcmc_result$average_copy_rate_pred_50[,1:ncol(average_copy_rate_mcmc_result$average_copy_rate_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
average_copy_rate_pred_data = rbind(average_copy_rate_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_copy_rate_mcmc_result$average_copy_rate_pred_78[,1:ncol(average_copy_rate_mcmc_result$average_copy_rate_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
average_copy_rate_pred_data$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(average_copy_rate_pred_data) <- c("groupSize_standerized2",
              "average_copy_rate_p2.5", "average_copy_rate_p25", "average_copy_rate_p50", "average_copy_rate_p75", "average_copy_rate_p97.5","taskDifficulty")

Figure4_b_25 = ggplot() +
  geom_ribbon(data=subset(average_copy_rate_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=average_copy_rate_p25,ymax=average_copy_rate_p75), fill = "grey80")+ #fill = "lightpink"
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&frequency_dependence_4!='Single-condition'&frequency_dependence_4!='--'&taskDifficulty=='25%'),mapping=aes(groupSize_standerized2, average_copy_rate,colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_path(data=subset(average_copy_rate_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=average_copy_rate_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  theme(legend.position = c(-0.05, 0.9), legend.text = element_text(size=9, family="Helvetica" ))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL
Figure4_b_50 = ggplot() +
  geom_ribbon(data=subset(average_copy_rate_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=average_copy_rate_p25,ymax=average_copy_rate_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='50%'),mapping=aes(groupSize_standerized2, average_copy_rate,colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_path(data=subset(average_copy_rate_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=average_copy_rate_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies") +
  labs(x='Group size\n (standardized)', y='', title='') +
  #labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL
Figure4_b_78 = ggplot() +
  geom_ribbon(data=subset(average_copy_rate_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, ymin=average_copy_rate_p25,ymax=average_copy_rate_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='78%'),mapping=aes(groupSize_standerized2, average_copy_rate,colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_path(data=subset(average_copy_rate_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, y=average_copy_rate_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies") +
  labs(x='Group size\n (standardized)', y='', title='') +
  #labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL

#Figure4_b = plot_grid(Figure4_b_25, Figure4_b_50, Figure4_b_78, labels = c('','',''), ncol = 3, align = 'v')

## conformity exponent
conformity_exponent_mcmc_result <- rstan::extract(conformity_exponent_fitting)
conformity_exponent_pred_data = cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(conformity_exponent_mcmc_result$conformity_exponent_pred_25[,1:ncol(conformity_exponent_mcmc_result$conformity_exponent_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
conformity_exponent_pred_data = rbind(conformity_exponent_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(conformity_exponent_mcmc_result$conformity_exponent_pred_50[,1:ncol(conformity_exponent_mcmc_result$conformity_exponent_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
conformity_exponent_pred_data = rbind(conformity_exponent_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(conformity_exponent_mcmc_result$conformity_exponent_pred_78[,1:ncol(conformity_exponent_mcmc_result$conformity_exponent_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
conformity_exponent_pred_data$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(conformity_exponent_pred_data) <- c("groupSize_standerized2",
              "conformity_exponent_p2.5", "conformity_exponent_p25", "conformity_exponent_p50", "conformity_exponent_p75", "conformity_exponent_p97.5","taskDifficulty")

Figure4_c_25 = ggplot() +
	geom_ribbon(data=subset(conformity_exponent_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=conformity_exponent_p25,ymax=conformity_exponent_p75), fill = "grey80")+#fill = "lightpink"
	geom_hline(yintercept=1, linetype='solid', size=0.7, alpha=0.7)+
	geom_hline(yintercept=-1, linetype='solid', size=0.7, alpha=0.7)+
	geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&frequency_dependence_4!='Single-condition'&frequency_dependence_4!='--'&taskDifficulty=='25%'),mapping=aes(groupSize_standerized2, theta, colour=frequency_dependence_4, shape=frequency_dependence_4))+
	geom_line(data=subset(conformity_exponent_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=conformity_exponent_p50), linetype='dashed', size=0.6,color='grey20')+
	scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
	scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
	labs(x='Group size\n (standardized)', y=bquote(atop('Conformity', 'exponent' ~ theta[i])), title='') +
	myTheme()+ylim(c(-7,12))+xlim(c(-1.3,2.3))+
	theme(legend.position = c(-0.05, 0.9), legend.text = element_text(size=9, family="Helvetica" ))+
	#panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL
Figure4_c_50 = ggplot() +
	geom_ribbon(data=subset(conformity_exponent_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=conformity_exponent_p25,ymax=conformity_exponent_p75), fill = "grey80")+
	geom_hline(yintercept=1, linetype='solid', size=0.7, alpha=0.7)+
	geom_hline(yintercept=-1, linetype='solid', size=0.7, alpha=0.7)+
	geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='50%'),mapping=aes(groupSize_standerized2, theta, colour=frequency_dependence_4, shape=frequency_dependence_4))+
	geom_line(data=subset(conformity_exponent_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=conformity_exponent_p50), linetype='dashed', size=0.6,color='grey20')+
	scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies") +
	scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies") +
	#labs(x='Group size\n (standardized)', y=bquote(atop('Conformity', 'exponent' ~ theta[i])), title='') +
  labs(x='Group size\n (standardized)', y='', title='') +
	myTheme()+ylim(c(-7,12))+xlim(c(-1.3,2.3))+
	#panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL
Figure4_c_78 = ggplot() +
	geom_ribbon(data=subset(conformity_exponent_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, ymin=conformity_exponent_p25,ymax=conformity_exponent_p75), fill = "grey80")+
	geom_hline(yintercept=1, linetype='solid', size=0.7, alpha=0.7)+
	geom_hline(yintercept=-1, linetype='solid', size=0.7, alpha=0.7)+
	geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='78%'),mapping=aes(groupSize_standerized2, theta, colour=frequency_dependence_4, shape=frequency_dependence_4))+
	geom_line(data=subset(conformity_exponent_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, y=conformity_exponent_p50), linetype='dashed', size=0.6,color='grey20')+
	scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies") +
	scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies") +
	#labs(x='Group size\n (standardized)', y=bquote(atop('Conformity', 'exponent' ~ theta[i])), title='') +
  labs(x='Group size\n (standardized)', y='', title='') +
	myTheme()+ylim(c(-7,12))+xlim(c(-1.3,2.3))+
	#panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL


Figure4_rev_right = plot_grid(GS_proportionSL_25,Figure4_b_25,Figure4_c_25, labels=c("","",""), ncol = 1, align = 'v')
Figure4_rev_left = plot_grid(GS_proportionSL_50,GS_proportionSL_78, Figure4_b_50,Figure4_b_78, Figure4_c_50,Figure4_c_78, labels=c("","",""), ncol = 2, align = 'v')
Figure4_rev = plot_grid(Figure4_rev_right,Figure4_rev_left, labels=c("","",""), ncol = 2, rel_widths = c(1,1.8))

ggsave(file = "Figure4_rev.pdf", plot = Figure4_rev, dpi = 600, width = 9, height = 7.5)



































#######################################################################
##
## Analysis for Figure 3
## Experimental Results --
## Time evolutions and distributions of decision performance for each condition.
##
#######################################################################
## _25%
onlineExp_stanData_25 = list(
  All = allBehaviouralData_25 %>% nrow(),
  T = allBehaviouralData2 %>% subset(taskDifficulty=='25%') %>% dplyr::pull(T2),
  maxT = 70,
  Group = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%') %>% dplyr::pull(room) %>% as.factor() %>% as.numeric(),
  Indv = allBehaviouralData_25 %>% dplyr::pull(amazonID) %>% as.factor() %>% as.numeric(),
  indiv_group = allBehaviouralData2 %>% subset(taskDifficulty=='25%') %>% dplyr::pull(indiv_group),
  #GroupSize_st = (allBehaviouralData_25$sizeCategory-mean(allBehaviouralData_25$sizeCategory))/sd(allBehaviouralData_25$sizeCategory),
  N_group = allBehaviouralData2 %>% subset(taskDifficulty=='25%') %>% dplyr::select(room) %>% dplyr::count(room) %>% spread(room, n) %>% length(),
  N_indv = allBehaviouralData2 %>% subset(taskDifficulty=='25%') %>% dplyr::select(amazonID) %>% dplyr::count(amazonID) %>% spread(amazonID, n) %>% length(),
  Y = allBehaviouralData2 %>% subset(taskDifficulty=='25%') %>% dplyr::pull(optimalChoice)
  )
## _50%
onlineExp_stanData_50 = list(
  All = allBehaviouralData_50 %>% nrow(),
  T = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::pull(T2),
  maxT = 70,
  Group = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%') %>% dplyr::pull(room) %>% as.factor() %>% as.numeric(),
  Indv = allBehaviouralData_50 %>% dplyr::pull(amazonID) %>% as.factor() %>% as.numeric(),
  indiv_group = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::pull(indiv_group),
  #GroupSize_st = (allBehaviouralData_50$sizeCategory-mean(allBehaviouralData_50$sizeCategory))/sd(allBehaviouralData_50$sizeCategory),
  N_group = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::select(room) %>% dplyr::count(room) %>% spread(room, n) %>% length(),
  N_indv = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::select(amazonID) %>% dplyr::count(amazonID) %>% spread(amazonID, n) %>% length(),
  Y = allBehaviouralData2 %>% subset(taskDifficulty=='50%') %>% dplyr::pull(optimalChoice)
  )
## _78%
onlineExp_stanData_78 = list(
  All = allBehaviouralData_78 %>% nrow(),
  T = allBehaviouralData2 %>% subset(taskDifficulty=='78%') %>% dplyr::pull(T2),
  maxT = 70,
  Group = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%') %>% dplyr::pull(room) %>% as.factor() %>% as.numeric(),
  Indv = allBehaviouralData_78 %>% dplyr::pull(amazonID) %>% as.factor() %>% as.numeric(),
  indiv_group = allBehaviouralData2 %>% subset(taskDifficulty=='78%') %>% dplyr::pull(indiv_group),
  #GroupSize_st = (allBehaviouralData_78$sizeCategory-mean(allBehaviouralData_78$sizeCategory))/sd(allBehaviouralData_78$sizeCategory),
  N_group = allBehaviouralData2 %>% subset(taskDifficulty=='78%') %>% dplyr::select(room) %>% dplyr::count(room) %>% spread(room, n) %>% length(),
  N_indv = allBehaviouralData2 %>% subset(taskDifficulty=='78%') %>% dplyr::select(amazonID) %>% dplyr::count(amazonID) %>% spread(amazonID, n) %>% length(),
  Y = allBehaviouralData2 %>% subset(taskDifficulty=='78%') %>% dplyr::pull(optimalChoice)
  )


#################################################################################################
###
### Statistical analysis on the average performance between conditons and different group sizes
###
### Summarised performance analysis 2 --
### -- With effect of group size and social learning
###
### Figure 3 A-B (revision)
#################################################################################################

## Compiling the model
stanmodel_onlineExp = stan_model(file='stan_model6/timeSeriesDiffTest_optimRate_groupSize.stan')


## _25%
## Counting social learner
#allBehaviouralData_25$pos_freq_dep = 0
#socialLearnersName_25 = performanceSummary4 %>% dplyr::filter(frequency_dependence_4 != '--' & frequency_dependence_4 != 'Single-condition' & taskDifficulty=='25%') %>% dplyr::count(amazonID) %>% spread(amazonID, n) %>% names()
#allBehaviouralData_25$pos_freq_dep[which(allBehaviouralData_25$amazonID %in% socialLearnersName_25)] = 1

onlineExp_stanData_25$GroupSize_st = (
  ( (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%') %>% dplyr::pull(sizeCategory)) -
    (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))

onlineExp_stanData_25$smallGroupSize = (
  ( 6 -
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))
onlineExp_stanData_25$largeGroupSize = (
  ( 14 -
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))
onlineExp_stanData_25$veryLargeGroupSize = (
  ( 20 -
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))

## Fitting
## You can skip the following model fitting if you already have done this before. In that case, you can load RDS file:
## onlineExp_fitting_groupSize_25 <- readRDS("stan_model6/mcmc_result/onlineExp_fitting_groupSize_25.rds")

onlineExp_fitting_groupSize_25 = sampling(
  stanmodel_onlineExp, data=onlineExp_stanData_25, seed=72,
  #control = list(adapt_delta = 0.80),
  pars=c('sizeEffect','s_i','mu0','di0','mu','di','p_indiv','p_group','p_group_small','p_group_large','p_group_veryLarge'),
  init=function() {
    list(mu0=runif(1,-3,-1))
  },
  chains=8, iter=2000, warmup=1000, thin=5
)
onlineExp_fitting_groupSize_25

saveRDS(onlineExp_fitting_groupSize_25, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/onlineExp_fitting_groupSize_25.rds')

trace_onlineExp_fitting_groupSize_25 = traceplot(onlineExp_fitting_groupSize_25, pars=c('sizeEffect', 's_i', 'mu0', 'di0'), inc_warmup=TRUE)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/trace_onlineExp_fitting_groupSize_25.pdf", plot = trace_onlineExp_fitting_groupSize_25, dpi = 600, width = 9, height = 9)

#plot(onlineExp_fitting_groupSize_25, pars=c('p_indiv','p_group_small'))
plot(onlineExp_fitting_groupSize_25, pars=c('mu','di'))
plot(onlineExp_fitting_groupSize_25, pars=c('sizeEffect','s_i'))

data_onlineExp_fitting_groupSize_25 <- rstan::extract(onlineExp_fitting_groupSize_25)
fittedValues_onlineExp_fitting_groupSize_25 =
  as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$mu[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$di[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_indiv[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_group_small[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_group_large[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_group_veryLarge[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_25$p_group[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
## rename
colnames(fittedValues_onlineExp_fitting_groupSize_25) = c("mu_p2.5", "mu_p25", "mu_p50", "mu_p75", "mu_p97.5","di_p2.5", "di_p25", "di_p50", "di_p75", "di_p97.5","p_indiv_p2.5", "p_indiv_p25", "p_indiv_p50", "p_indiv_p75", "p_indiv_p97.5","p_group_small_p2.5", "p_group_small_p25", "p_group_small_p50", "p_group_small_p75", "p_group_small_p97.5","p_group_large_p2.5", "p_group_large_p25", "p_group_large_p50", "p_group_large_p75", "p_group_large_p97.5","p_group_veryLarge_p2.5", "p_group_veryLarge_p25", "p_group_veryLarge_p50", "p_group_veryLarge_p75", "p_group_veryLarge_p97.5","p_group_p2.5", "p_group_p25", "p_group_p50", "p_group_p75", "p_group_p97.5")
fittedValues_onlineExp_fitting_groupSize_25$taskDifficulty = "25%"




## _50%
## Counting social learner
#allBehaviouralData_50$pos_freq_dep = 0
#socialLearnersName_50 = performanceSummary4 %>% dplyr::filter(frequency_dependence_4 != '--' & frequency_dependence_4 != 'Single-condition' & taskDifficulty=='50%') %>% dplyr::count(amazonID) %>% spread(amazonID, n) %>% names()
#allBehaviouralData_50$pos_freq_dep[which(allBehaviouralData_50$amazonID %in% socialLearnersName_50)] = 1

onlineExp_stanData_50$GroupSize_st = (
  ( (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%') %>% dplyr::pull(sizeCategory)) -
    (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))

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

## Fitting
## You can skip the following model fitting if you already have done this before. In that case, you can load RDS file:
## onlineExp_fitting_groupSize_50 <- readRDS("stan_model6/mcmc_result/onlineExp_fitting_groupSize_50.rds")
onlineExp_fitting_groupSize_50 = sampling(
  stanmodel_onlineExp, data=onlineExp_stanData_50, seed=72,
  #control = list(adapt_delta = 0.80),
  pars=c('sizeEffect','s_i','mu0','di0','mu','di','p_indiv','p_group','p_group_small','p_group_large','p_group_veryLarge'),
  init=function() {
    list(mu0=runif(1,-3,-1))
  },
  chains=8, iter=2000, warmup=1000, thin=5
)
onlineExp_fitting_groupSize_50

saveRDS(onlineExp_fitting_groupSize_50, file='mcmc_result/onlineExp_fitting_groupSize_50.rds')

trace_onlineExp_fitting_groupSize_50 = traceplot(onlineExp_fitting_groupSize_50, pars=c('sizeEffect', 's_i', 'mu0', 'di0'), inc_warmup=TRUE)
ggsave(file = "mcmc_result/trace_onlineExp_fitting_groupSize_50.pdf", plot = trace_onlineExp_fitting_groupSize_50, dpi = 600, width = 9, height = 9)

#plot(onlineExp_fitting_groupSize_50, pars=c('p_indiv','p_group_small'))
plot(onlineExp_fitting_groupSize_50, pars=c('mu','di'))
plot(onlineExp_fitting_groupSize_50, pars=c('sizeEffect','s_i'))

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






## _78%
## Counting social learner
#allBehaviouralData_78$pos_freq_dep = 0
#socialLearnersName_78 = performanceSummary4 %>% dplyr::filter(frequency_dependence_4 != '--' & frequency_dependence_4 != 'Single-condition' & taskDifficulty=='78%') %>% dplyr::count(amazonID) %>% spread(amazonID, n) %>% names()
#allBehaviouralData_78$pos_freq_dep[which(allBehaviouralData_78$amazonID %in% socialLearnersName_78)] = 1

onlineExp_stanData_78$GroupSize_st = (
  ( (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%') %>% dplyr::pull(sizeCategory)) -
    (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))

onlineExp_stanData_78$smallGroupSize = (
  ( 4 -
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))
onlineExp_stanData_78$largeGroupSize = (
  ( 15 -
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))
onlineExp_stanData_78$veryLargeGroupSize = (
  ( 20 -
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% mean()))/
  (allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%'&sizeCategory>1) %>% dplyr::pull(sizeCategory) %>% sd()))

## Fitting
## You can skip the following model fitting if you already have done this before. In that case, you can load RDS file:
## onlineExp_fitting_groupSize_78 <- readRDS("stan_model6/mcmc_result/onlineExp_fitting_groupSize_78.rds")
onlineExp_fitting_groupSize_78 = sampling(
  stanmodel_onlineExp, data=onlineExp_stanData_78, seed=72,
  #control = list(adapt_delta = 0.80),
  pars=c('sizeEffect','s_i','mu0','di0','mu','di','p_indiv','p_group','p_group_small','p_group_large','p_group_veryLarge'),
  init=function() {
    list(mu0=runif(1,-3,-1))
  },
  chains=8, iter=2000, warmup=1000, thin=5
)
onlineExp_fitting_groupSize_78

saveRDS(onlineExp_fitting_groupSize_78, file='mcmc_result/onlineExp_fitting_groupSize_78.rds')

trace_onlineExp_fitting_groupSize_78 = traceplot(onlineExp_fitting_groupSize_78, pars=c('sizeEffect', 's_i', 'mu0', 'di0'), inc_warmup=TRUE)
ggsave(file = "mcmc_result/trace_onlineExp_fitting_groupSize_78.pdf", plot = trace_onlineExp_fitting_groupSize_78, dpi = 600, width = 9, height = 9)

#plot(onlineExp_fitting_groupSize_78, pars=c('p_indiv','p_group_small'))
plot(onlineExp_fitting_groupSize_78, pars=c('mu','di'))
plot(onlineExp_fitting_groupSize_78, pars=c('sizeEffect','s_i'))

data_onlineExp_fitting_groupSize_78 <- rstan::extract(onlineExp_fitting_groupSize_78)
fittedValues_onlineExp_fitting_groupSize_78 =
  as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$mu[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$di[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_indiv[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_group_small[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_group_large[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_group_veryLarge[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_78$p_group[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
## rename
colnames(fittedValues_onlineExp_fitting_groupSize_78) = c("mu_p2.5", "mu_p25", "mu_p50", "mu_p75", "mu_p97.5","di_p2.5", "di_p25", "di_p50", "di_p75", "di_p97.5","p_indiv_p2.5", "p_indiv_p25", "p_indiv_p50", "p_indiv_p75", "p_indiv_p97.5","p_group_small_p2.5", "p_group_small_p25", "p_group_small_p50", "p_group_small_p75", "p_group_small_p97.5","p_group_large_p2.5", "p_group_large_p25", "p_group_large_p50", "p_group_large_p75", "p_group_large_p97.5","p_group_veryLarge_p2.5", "p_group_veryLarge_p25", "p_group_veryLarge_p50", "p_group_veryLarge_p75", "p_group_veryLarge_p97.5","p_group_p2.5", "p_group_p25", "p_group_p50", "p_group_p75", "p_group_p97.5")
fittedValues_onlineExp_fitting_groupSize_78$taskDifficulty = "78%"






## Investigating the post-hoc simulation
## How likely it is that larger groups can perform better than smaller or solitaries?
sampleSize = 10000
### 25% 1st - 40th 'Single' vs 'Small'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='25%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='25%' & indiv_small_large=='Small') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 25% 41st - 70th 'Single' vs 'Small'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='25%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='25%' & indiv_small_large=='Small') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 25% 1st - 40th 'Single' vs 'Large'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='25%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='25%' & indiv_small_large=='Large') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 25% 41st - 70th 'Single' vs 'Large'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='25%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='25%' & indiv_small_large=='Large') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 25% 41st - 70th 'Small' vs 'Large'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='25%' & indiv_small_large=='Small') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='25%' & indiv_small_large=='Large') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 50% 1st - 40th 'Single' vs 'Small'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='50%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='50%' & indiv_small_large=='Small') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 50% 41st - 70th 'Single' vs 'Small'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='50%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='50%' & indiv_small_large=='Small') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 50% 1st - 40th 'Single' vs 'Large'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='50%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='50%' & indiv_small_large=='Large') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 50% 41st - 70th 'Single' vs 'Large'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='50%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='50%' & indiv_small_large=='Large') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 50% 41st - 70th 'Small' vs 'Large'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='50%' & indiv_small_large=='Small') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='50%' & indiv_small_large=='Large') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 78% 1st - 40th 'Single' vs 'Small'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='78%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='78%' & indiv_small_large=='Small') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 78% 41st - 70th 'Single' vs 'Small'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='78%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='78%' & indiv_small_large=='Small') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 78% 1st - 40th 'Single' vs 'Large'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='78%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='1st - 40th' & taskDifficulty=='78%' & indiv_small_large=='Large') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 78% 41st - 70th 'Single' vs 'Large'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='78%' & indiv_small_large=='Single') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='78%' & indiv_small_large=='Large') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize
### 78% 41st - 70th 'Small' vs 'Large'
(sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='78%' & indiv_small_large=='Small') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE) < sample(posthoc_sim_eachGroup_aveOptim %>% dplyr::filter(environment=='41st - 70th' & taskDifficulty=='78%' & indiv_small_large=='Large') %>% dplyr::pull(meanOptimProb),sampleSize,replace=TRUE)) %>% table()/sampleSize

# Is the shape likely to be a normal distribution?
#shapiro.test()








#######################################################################
##
## Figure 5 - social learning weight changing
## Model fitting results
##
#######################################################################


## Time evolution of the Social Learning Weight
copyingProbReduction_data_fullModelOnly = data.frame( amazonID = rep(temp_UNC6_sReduc_annealing$amazonID,70)
                                    ,   taskDifficulty = rep(temp_UNC6_sReduc_annealing$taskDifficulty,70)
                                    ,   correctProp = rep(temp_UNC6_sReduc_annealing$correctProp,70)
                                    ,   groupSize = NA
                                    ,   round = NA
                                    ,   best_model = NA
                                    ,   alpha = NA
                                    ,   beta = NA
                                    ,   beta_low = NA
                                    ,   beta_high = NA
                                    ,   annealing = NA
                                    ,   theta = NA
                                    ,   low_theta = NA
                                    ,   high_theta = NA
                                    ,   copyingProb_low = NA
                                    ,   copyingProb = NA
                                    ,   copyingProb_high = NA
                                    ,   theta_25 = NA
                                    ,   theta_75 = NA
                                    )
for(t in 1:70) {
    copyingProbReduction_data_fullModelOnly$round[(1+(t-1)*nrow(temp_UNC6_sReduc_annealing)):(t*nrow(temp_UNC6_sReduc_annealing))] <- t
}

for(i in names(table(copyingProbReduction_data_fullModelOnly$amazonID))) {
    #thisSubjectModel = as.character(performanceSummary4$smallestWaicModel[performanceSummary4$subjectID==i])[1]
    copyingProbReduction_data_fullModelOnly$best_model[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] = "UNC6_sReduc_ann"
    #if(thisSubjectModel == 'UNC6_sReduc_ann') {
        # copying rate
        for(t in 1:70) {
            copyingProbReduction_data_fullModelOnly$copyingProb[which(copyingProbReduction_data_fullModelOnly$amazonID==i&copyingProbReduction_data_fullModelOnly$round==t)] <- soc_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_p50[which(soc_table_UNC6_sReduc_annealing$round==t&soc_table_UNC6_sReduc_annealing$amazonID==i)]
            copyingProbReduction_data_fullModelOnly$copyingProb_low[which(copyingProbReduction_data_fullModelOnly$amazonID==i&copyingProbReduction_data_fullModelOnly$round==t)] <- soc_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_p2.5[which(soc_table_UNC6_sReduc_annealing$round==t&soc_table_UNC6_sReduc_annealing$amazonID==i)]
            copyingProbReduction_data_fullModelOnly$copyingProb_high[which(copyingProbReduction_data_fullModelOnly$amazonID==i&copyingProbReduction_data_fullModelOnly$round==t)] <- soc_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_p97.5[which(soc_table_UNC6_sReduc_annealing$round==t&soc_table_UNC6_sReduc_annealing$amazonID==i)]
        }
        # netBeta
        for(t in 1:70) {
            copyingProbReduction_data_fullModelOnly$beta[which(copyingProbReduction_data_fullModelOnly$amazonID==i&copyingProbReduction_data_fullModelOnly$round==t)] <- netBeta_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_netBeta_p50[which(netBeta_table_UNC6_sReduc_annealing$round==t&netBeta_table_UNC6_sReduc_annealing$amazonID==i)]
            copyingProbReduction_data_fullModelOnly$beta_low[which(copyingProbReduction_data_fullModelOnly$amazonID==i&copyingProbReduction_data_fullModelOnly$round==t)] <- netBeta_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_netBeta_p2.5[which(netBeta_table_UNC6_sReduc_annealing$round==t&netBeta_table_UNC6_sReduc_annealing$amazonID==i)]
            copyingProbReduction_data_fullModelOnly$beta_high[which(copyingProbReduction_data_fullModelOnly$amazonID==i&copyingProbReduction_data_fullModelOnly$round==t)] <- netBeta_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_netBeta_p97.5[which(netBeta_table_UNC6_sReduc_annealing$round==t&netBeta_table_UNC6_sReduc_annealing$amazonID==i)]
        }
        # other parameters
        copyingProbReduction_data_fullModelOnly$alpha[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] = performanceSummary4$alpha[which(performanceSummary4$subjectID==i)][1]
        #copyingProbReduction_data_fullModelOnly$beta[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] = performanceSummary4$beta[which(performanceSummary4$subjectID==i)][1]
        copyingProbReduction_data_fullModelOnly$annealing[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] = performanceSummary4$annealing[which(performanceSummary4$subjectID==i)][1]
        copyingProbReduction_data_fullModelOnly$theta[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] = performanceSummary4$theta[which(performanceSummary4$subjectID==i)][1]
        copyingProbReduction_data_fullModelOnly$low_theta[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] = performanceSummary4$low_theta[which(performanceSummary4$subjectID==i)][1]
        copyingProbReduction_data_fullModelOnly$high_theta[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] = performanceSummary4$high_theta[which(performanceSummary4$subjectID==i)][1]
        copyingProbReduction_data_fullModelOnly$theta_25[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] = performanceSummary4$theta_25[which(performanceSummary4$subjectID==i)][1]
        copyingProbReduction_data_fullModelOnly$theta_75[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] = performanceSummary4$theta_75[which(performanceSummary4$subjectID==i)][1]
    #}
}

    # group size
for(i in names(table(copyingProbReduction_data_fullModelOnly$amazonID))){
    copyingProbReduction_data_fullModelOnly$groupSize[which(copyingProbReduction_data_fullModelOnly$amazonID==i)] <- performanceSummary4$groupSize[which(performanceSummary4$subjectID==i)][1]
}

# conservative categories (rand, neg, pos) - 25-75
copyingProbReduction_data_fullModelOnly$frequency_dependence_4 = 'Random-copying'
copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$theta_75 < 0)] = 'neg-freq-dep'
copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$theta_25 > 0)] = 'pos-freq-dep'
#copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$best_model == 'AL' | copyingProbReduction_data_fullModelOnly$best_model == 'AL_ann')] = 'AL'
#copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$best_model == 'Random')] = 'Random'
copyingProbReduction_data_fullModelOnly$frequency_dependence_4 = factor(copyingProbReduction_data_fullModelOnly$frequency_dependence_4, levels=c("Random choice","AL","Random-copying",'neg-freq-dep', "pos-freq-dep"))

## Task difficulty 2
copyingProbReduction_data_fullModelOnly$taskDifficulty2 = NA
copyingProbReduction_data_fullModelOnly$taskDifficulty2[which(copyingProbReduction_data_fullModelOnly$taskDifficulty=='25%')] <- 'Low Uncertainty'
copyingProbReduction_data_fullModelOnly$taskDifficulty2[which(copyingProbReduction_data_fullModelOnly$taskDifficulty=='50%')] <- 'Moderate Uncertainty'
copyingProbReduction_data_fullModelOnly$taskDifficulty2[which(copyingProbReduction_data_fullModelOnly$taskDifficulty=='78%')] <- 'High Uncertainty'
copyingProbReduction_data_fullModelOnly$taskDifficulty2 = factor(copyingProbReduction_data_fullModelOnly$taskDifficulty2, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))

write.csv(copyingProbReduction_data_fullModelOnly,
      "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/copyingProbReduction_data_fullModelOnly.csv",
      row.names=FALSE)

# Using 50% CI for categorisation
# Three categories: Random-copying, neg-freq-dep, pos-freq-dep
# Presenting only the positive frequency-dependent copiers

Figure4_toprow = ggplot() +
                    geom_line(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'&frequency_dependence_4=='pos-freq-dep'), mapping=aes(round,copyingProb,group=amazonID,colour=theta),alpha=3/4)+
                    stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'&frequency_dependence_4=='pos-freq-dep'),aes(round,copyingProb),fun.y = "median",geom="line",linetype = "dashed", colour = "black", size = 1)+
                    facet_grid(.~taskDifficulty2)+
                    labs(x='Round', y=expression(paste('Social learning weight ',sigma[i][t],sep="")), title='') +
                    ylim(c(0,1))+
                    xlim(c(2,70))+
                    #scale_colour_distiller(name=expression(paste('Conformity\nexponent',theta[i],sep="")), palette = "RdYlBu", direction = -1)+
                    scale_color_viridis(name=expression(paste('Conformity\nexponent',theta[i],sep="")), option="magma", direction = -1)+
                    myTheme_legend()+
                    #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
                    NULL
ggsave(file = "Figure5_rev.pdf", plot = Figure4_toprow, dpi = 600, width = 7.5, height = 2.5)

round_copyingRate_theta_fullModel = ggplot() +
                    geom_line(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'), mapping=aes(round,copyingProb,group=amazonID,colour=frequency_dependence_4),alpha=2/4)+
                    stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,copyingProb),fun.y = "median",geom="line", colour = "black", size = 1)+
                    #stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,copyingProb_low),fun.y = "median",geom="line", colour = "grey50", size = 1)+
                    #stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,copyingProb_high),fun.y = "median",geom="line", colour = "grey50", size = 1)+
                    facet_grid(taskDifficulty2~frequency_dependence_4)+
                    labs(x='Round', y='Copying probability', title='') +
                    ylim(c(0,1))+
                    xlim(c(2,70))+
                    scale_colour_manual(values=c("random-choice"="#999999","AL"="#999999","Random-copying"="#999999",'neg-freq-dep'="blue", "pos-freq-dep"="red"), name="Strength of \nfrequency dependence") +
                    #scale_colour_distiller(name="Frequency \ndependence", palette = "Spectral", direction = -1)+
                    myTheme()+
                    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

ggsave(file = "round_copyingRate_theta_fullModel.pdf", plot = round_copyingRate_theta_fullModel, dpi = 600, width = 9*1.2, height = 9*1.2)


#ggsave(file = "testp.pdf", plot = testp, dpi = 600, width = 9*1.2, height = 3*1.2)
theta_soc_change_posOnly = ggplot() +
                    geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4=='Pos-freq-dep'), mapping=aes(theta, soc_change, colour=theta)) +
                    facet_grid(.~taskDifficulty3) +
                    labs(x=expression(paste('Conformity exponent ', theta[i],sep="")), y=expression(paste('Social learning weight\'s slope', delta[i], sep="")), title='') +
                    geom_hline(yintercept = 0, linetype='dashed')+
                    scale_color_viridis(name=expression(paste('Conformity\nexponent',theta[i],sep="")), option="magma", direction = -1)+
                    myTheme_legend()+
                    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

ggsave(file = "theta_soc_change_posOnly.pdf", plot = theta_soc_change_posOnly, dpi = 600, width = 9*1.2, height = 3*1.2)


















#######################################################################
##
## Figure 4
## Model Fitting Results
##
#######################################################################






###################################################################
## STATISTICAL ANALYSIS about Fig 3_top
## Does being positive copier depend on the factors?
###################################################################
parameters_stan_data$isPositiveCopier = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$isPositiveCopier

  ## FITTING
isPositiveCopier_stan_model = stan_model(file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/isPositiveCopier.stan') # debug
isPositiveCopier_fitting = sampling(
  isPositiveCopier_stan_model, data=parameters_stan_data, seed=77,
  #control = list(adapt_delta = 0.90, max_treedepth=13),
  pars=c('beta','sigma_r','r_group','isPositiveCopier_pred_25','isPositiveCopier_pred_50','isPositiveCopier_pred_78'),
  init=function() {
    list(
      sigma_r=runif(1,0,3)
      )
  },
  chains=8, iter=5000, warmup=1500, thin=2
  #chains=8, iter=500, warmup=150, thin=2
)
isPositiveCopier_fitting
saveRDS(isPositiveCopier_fitting, file='mcmc_result/isPositiveCopier_fitting.rds')

plot(isPositiveCopier_fitting, pars=c('beta','sigma_r','sigma_r'))
plot(isPositiveCopier_fitting, pars=c('isPositiveCopier_pred_25'))
  ## 事後診断
traceplot(isPositiveCopier_fitting, pars=c('beta','sigma_e'), inc_warmup=FALSE)

isPositiveCopier_mcmc_result <- rstan::extract(isPositiveCopier_fitting)
isPositiveCopier_pred_data = cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(isPositiveCopier_mcmc_result$isPositiveCopier_pred_25[,1:ncol(isPositiveCopier_mcmc_result$isPositiveCopier_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
isPositiveCopier_pred_data = rbind(isPositiveCopier_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(isPositiveCopier_mcmc_result$isPositiveCopier_pred_50[,1:ncol(isPositiveCopier_mcmc_result$isPositiveCopier_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
isPositiveCopier_pred_data = rbind(isPositiveCopier_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(isPositiveCopier_mcmc_result$isPositiveCopier_pred_78[,1:ncol(isPositiveCopier_mcmc_result$isPositiveCopier_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
isPositiveCopier_pred_data$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(isPositiveCopier_pred_data) <- c("groupSize_standerized2",
              "isPositiveCopier_p2.5", "isPositiveCopier_p25", "isPositiveCopier_p50", "isPositiveCopier_p75", "isPositiveCopier_p97.5","taskDifficulty")

lsfit(performanceSummary4$groupSize_standerized2, performanceSummary4$groupSize)$coefficients
#Intercept         X
#10.421192  7.808047

isPositiveCopier_pred_data$groupSize = isPositiveCopier_pred_data$groupSize_standerized2*7.808047 + 10.421192

## Figure 3 top row
GS_proportionSL_25 = ggplot() +
  geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_25, aes(groupSize, proportion_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_25,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_25,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  labs(x='Group size', y='Proportion of \nstrategies', title='Low Uncertainty') +
  myTheme_legend()+ylim(c(NA,1))+xlim(c(2, 28))+
  theme(legend.position = c(0.05, 0.5))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL
GS_proportionSL_50 = ggplot() +
  geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_50, aes(groupSize, proportion_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_50,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_50,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  #labs(x='Group size', y='Proportion of \nstrategies', title='Moderate Uncertainty') +
  labs(x='Group size', y='', title='Moderate Uncertainty') +
  myTheme_legend()+ylim(c(NA,1))+xlim(c(2, 28))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL
GS_proportionSL_78 = ggplot() +
  geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=2.2), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=2.2), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_78, aes(groupSize, proportion_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_78,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_78,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  #labs(x='Group size', y='Proportion of \nstrategies', title='High Uncertainty') +
  labs(x='Group size', y='', title='High Uncertainty') +
  myTheme_legend()+ylim(c(NA,1))+xlim(c(2, 28))+
  #panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
  NULL

# Summary plot
Figure3_top = plot_grid(GS_proportionSL_25, GS_proportionSL_50, GS_proportionSL_78, labels = c("","",""), ncol = 3, align = 'v')


## Figure S??
## Strategies' proportion using four categories (i.e. negative, random, weak-positive and strong-positive)
GS_proportionSL_supp_25 = ggplot() +
  #geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  #geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_25, aes(groupSize, proportion_strong_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_strong_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_25, aes(groupSize, proportion_weak_pos_freq_dep),colour='orange') + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_weak_pos_freq_dep),colour='orange',shape=2) +
  geom_line(data=strategies_proportion_25,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_25,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_25, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  labs(x='Group size', y='Proportion of \nstrategies', title='Low Uncertainty') +
  myTheme_inner()+ylim(c(NA,1))+xlim(c(2, 28))+
  theme(legend.position = c(0.05, 0.5))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
GS_proportionSL_supp_50 = ggplot() +
  #geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  #geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_50, aes(groupSize, proportion_strong_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_strong_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_50, aes(groupSize, proportion_weak_pos_freq_dep),colour='orange') + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_weak_pos_freq_dep),colour='orange',shape=2) +
  geom_line(data=strategies_proportion_50,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_50,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  #labs(x='Group size', y='Proportion of \nstrategies', title='Moderate Uncertainty') +
  labs(x='Group size', y='', title='Moderate Uncertainty') +
  myTheme()+ylim(c(NA,1))+xlim(c(2, 28))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
GS_proportionSL_supp_78 = ggplot() +
  #geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=2.2), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  #geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=2.2), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_78, aes(groupSize, proportion_strong_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_strong_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_78, aes(groupSize, proportion_weak_pos_freq_dep),colour='orange') + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_weak_pos_freq_dep),colour='orange',shape=2) +
  geom_line(data=strategies_proportion_78,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_78,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  #labs(x='Group size', y='Proportion of \nstrategies', title='High Uncertainty') +
  labs(x='Group size', y='', title='High Uncertainty') +
  myTheme()+ylim(c(NA,1))+xlim(c(2, 28))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

# Summary plot
FigureS_socialLearnersProp = plot_grid(GS_proportionSL_supp_25, GS_proportionSL_supp_50, GS_proportionSL_supp_78, labels = c("","",""), ncol = 3, align = 'v')
ggsave(file = "FigureS_socialLearnersProp.pdf", plot = FigureS_socialLearnersProp, dpi = 700, width = 10, height = 3.33)
#ggsave(file = "FigureS_socialLearnersProp.tiff", plot = FigureS_socialLearnersProp, dpi = 700, width = 10, height = 10)
###################################################################
## STATISTICAL ANALYSIS about Fig 3
## Does being positive copier depend on the factors?  -- END
###################################################################


###################################################################
## STATISTICAL ANALYSIS about Fig 3
## Average copying Rate
###################################################################





  ## average_copy_rate vs group size vs task difficulty
parameters_stan_data = list(
  All = nrow(subset(subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition'))),
  N_group = length(table((subset(subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition'))$groupID))),
  groupID = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$groupID_groupOnly,
  taskDifficulty = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$taskDifficulty_num,
  groupSize_standerized = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$groupSize_standerized2,
  age_standardized = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$age_standardized,
  female = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$female,
  average_copy_rate = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$average_copy_rate,
  test_group_size = seq(-1.3,2.3,len=200),
  test_size = 200
  )

  ## FITTING
average_copy_rate_stan_model = stan_model(file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/average_copy_rate.stan') # debug
average_copy_rate_fitting = sampling(
  average_copy_rate_stan_model, data=parameters_stan_data, seed=77,
  #control = list(adapt_delta = 0.90, max_treedepth=13),
  pars=c('beta','gamma','sigma_e','sigma_r','r_group','average_copy_rate_pred_25','average_copy_rate_pred_50','average_copy_rate_pred_78'),
  init=function() {
    list(
      sigma_e=runif(1,0,3),
      sigma_r=runif(1,0,3)
      )
  },
  chains=8, iter=5000, warmup=1500, thin=2
  #chains=8, iter=500, warmup=150, thin=2
)
average_copy_rate_fitting
saveRDS(average_copy_rate_fitting, file='mcmc_result/average_copy_rate_fitting.rds')

###################################################################
## STATISTICAL ANALYSIS about Fig 3
## Average copying Rate               -- END
###################################################################


################################
## Supporting figure:
## Analysis only on Pos-freq-dep
################################
performanceSummary4_stanData_posOnly = subset(performanceSummary4_stanData, (age_standardized > -100) & (female==0 | female==1) & frequency_dependence_4=='Pos-freq-dep')

  ### Setting a group level data
performanceSummary4_stanData_posOnly_group = data.frame(
  groupID = names(table(performanceSummary4_stanData_posOnly$groupID)))
performanceSummary4_stanData_posOnly_group$environment = 'e1+e2'
performanceSummary4_stanData_posOnly_group$taskDifficulty = NA
performanceSummary4_stanData_posOnly_group$average_entropy = NA
performanceSummary4_stanData_posOnly_group$groupSize = NA
#performanceSummary4_stanData_posOnly_group$totalPositiveCopier = NA
for(g in names(table(performanceSummary4_stanData_posOnly$groupID))){
  performanceSummary4_stanData_posOnly_group$taskDifficulty[which(performanceSummary4_stanData_posOnly_group$groupID==g)] = performanceSummary4_stanData_posOnly$taskDifficulty[which(performanceSummary4_stanData_posOnly$groupID==g)][1]
  performanceSummary4_stanData_posOnly_group$average_entropy[which(performanceSummary4_stanData_posOnly_group$groupID==g)] = performanceSummary4_stanData_posOnly$average_entropy[which(performanceSummary4_stanData_posOnly$groupID==g)][1]
  performanceSummary4_stanData_posOnly_group$groupSize[which(performanceSummary4_stanData_posOnly_group$groupID==g)] = performanceSummary4_stanData_posOnly$groupSize[which(performanceSummary4_stanData_posOnly$groupID==g)][1]
  #performanceSummary4_stanData_posOnly_group$totalPositiveCopier[which(performanceSummary4_stanData_posOnly_group$groupID==g)] = performanceSummary4_stanData_posOnly$totalPositiveCopier[which(performanceSummary4_stanData_posOnly$groupID==g)][1]
}
  ## Corlelations between social learning weight and group sizes
performanceSummary4_stanData_posOnly_groupOnly = subset(performanceSummary4_stanData_posOnly_group, groupSize>1)
performanceSummary4_stanData_posOnly_groupOnly = subset(performanceSummary4_stanData_posOnly_groupOnly, groupID != 'subRoom1_1_6') # no data
performanceSummary4_stanData_posOnly_groupOnly$groupID_groupOnly = 1:nrow(performanceSummary4_stanData_posOnly_groupOnly)
performanceSummary4_stanData_posOnly_groupOnly$groupID = as.character(performanceSummary4_stanData_posOnly_groupOnly$groupID)
performanceSummary4_stanData_posOnly$groupID_groupOnly = NA
for(g in names(table(performanceSummary4_stanData_posOnly_groupOnly$groupID))) {
  performanceSummary4_stanData_posOnly$groupID_groupOnly[which(performanceSummary4_stanData_posOnly$groupID==g)] <- performanceSummary4_stanData_posOnly_groupOnly$groupID_groupOnly[which(performanceSummary4_stanData_posOnly_groupOnly$groupID==g)]
}
performanceSummary4_stanData_posOnly$taskDifficulty_num = 0
performanceSummary4_stanData_posOnly$taskDifficulty_num[which(performanceSummary4_stanData_posOnly$taskDifficulty=='50%')] = 0.5
performanceSummary4_stanData_posOnly$taskDifficulty_num[which(performanceSummary4_stanData_posOnly$taskDifficulty=='78%')] = 1.0



  ## average_copy_rate vs group size vs task difficulty
parameters_stan_data_posOnly = list(
  All = nrow(subset(subset(performanceSummary4_stanData_posOnly, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition'))),
  N_group = length(table((subset(subset(performanceSummary4_stanData_posOnly, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition'))$groupID))),
  groupID = subset(performanceSummary4_stanData_posOnly, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$groupID_groupOnly,
  taskDifficulty = subset(performanceSummary4_stanData_posOnly, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$taskDifficulty_num,
  groupSize_standerized = subset(performanceSummary4_stanData_posOnly, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$groupSize_standerized2,
  age_standardized = subset(performanceSummary4_stanData_posOnly, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$age_standardized,
  female = subset(performanceSummary4_stanData_posOnly, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$female,
  average_copy_rate = subset(performanceSummary4_stanData_posOnly, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$average_copy_rate,
  test_group_size = seq(-1.3,2.3,len=200),
  test_size = 200
  )

  ## FITTING
average_copy_rate_stan_model = stan_model(file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/average_copy_rate.stan') # debug
average_copy_rate_fitting_posOnly = sampling(
  average_copy_rate_stan_model, data=parameters_stan_data_posOnly, seed=77,
  #control = list(adapt_delta = 0.90, max_treedepth=13),
  pars=c('beta','gamma','sigma_e','sigma_r','r_group','average_copy_rate_pred_25','average_copy_rate_pred_50','average_copy_rate_pred_78'),
  init=function() {
    list(
      sigma_e=runif(1,0,3),
      sigma_r=runif(1,0,3)
      )
  },
  chains=8, iter=5000, warmup=1500, thin=4
  #chains=8, iter=500, warmup=150, thin=2
)
average_copy_rate_fitting_posOnly
saveRDS(average_copy_rate_fitting_posOnly, file='mcmc_result/average_copy_rate_fitting_posOnly.rds')
average_copy_rate_mcmc_result_posOnly <- rstan::extract(average_copy_rate_fitting_posOnly)
plot(average_copy_rate_fitting_posOnly, pars=c('beta','gamma','sigma_e','sigma_r'))
#plot(average_copy_rate_fitting_posOnly, pars=c('average_copy_rate_pred_25'))

  ## 事後診断
traceplot(average_copy_rate_fitting_posOnly, pars=c('beta','sigma_e','sigma_r'), inc_warmup=FALSE)

average_copy_rate_pred_data_posOnly = cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_copy_rate_mcmc_result_posOnly$average_copy_rate_pred_25[,1:ncol(average_copy_rate_mcmc_result_posOnly$average_copy_rate_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
average_copy_rate_pred_data_posOnly = rbind(average_copy_rate_pred_data_posOnly, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_copy_rate_mcmc_result_posOnly$average_copy_rate_pred_50[,1:ncol(average_copy_rate_mcmc_result_posOnly$average_copy_rate_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
average_copy_rate_pred_data_posOnly = rbind(average_copy_rate_pred_data_posOnly, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_copy_rate_mcmc_result_posOnly$average_copy_rate_pred_78[,1:ncol(average_copy_rate_mcmc_result_posOnly$average_copy_rate_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
average_copy_rate_pred_data_posOnly$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(average_copy_rate_pred_data_posOnly) <- c("groupSize_standerized2",
              "average_copy_rate_p2.5", "average_copy_rate_p25", "average_copy_rate_p50", "average_copy_rate_p75", "average_copy_rate_p97.5","taskDifficulty")

FigureS_copyingRate_posOnly_25 = ggplot() +
  geom_ribbon(data=subset(average_copy_rate_pred_data_posOnly,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=average_copy_rate_p25,ymax=average_copy_rate_p75), fill = "lightpink")+ #fill = "lightpink"
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&frequency_dependence_4!='Single-condition'&frequency_dependence_4!='--'&taskDifficulty=='25%'),mapping=aes(groupSize_standerized2, average_copy_rate,colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_path(data=subset(average_copy_rate_pred_data_posOnly,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=average_copy_rate_p50), linetype='dashed', size=0.6,color='red')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='Low Uncertainty') +
  myTheme_inner()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  theme(legend.position = c(0.02, 0.75))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
FigureS_copyingRate_posOnly_50 = ggplot() +
  geom_ribbon(data=subset(average_copy_rate_pred_data_posOnly,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=average_copy_rate_p25,ymax=average_copy_rate_p75), fill = "lightpink")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='50%'),mapping=aes(groupSize_standerized2, average_copy_rate,colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_path(data=subset(average_copy_rate_pred_data_posOnly,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=average_copy_rate_p50), linetype='dashed', size=0.6,color='red')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='Moderate Uncertainty') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
FigureS_copyingRate_posOnly_78 = ggplot() +
  geom_ribbon(data=subset(average_copy_rate_pred_data_posOnly,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, ymin=average_copy_rate_p25,ymax=average_copy_rate_p75), fill = "lightpink")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='78%'),mapping=aes(groupSize_standerized2, average_copy_rate,colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_path(data=subset(average_copy_rate_pred_data_posOnly,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, y=average_copy_rate_p50), linetype='dashed', size=0.6,color='red')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='High Uncertainty') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

FigureS_copyingRate_posOnly = plot_grid(FigureS_copyingRate_posOnly_25, FigureS_copyingRate_posOnly_50, FigureS_copyingRate_posOnly_78, labels = c('','',''), ncol = 3, align = 'v')


################################
## Supporting figure:
## Analysis only on Pos-freq-dep - END
################################



###################################################################
## STATISTICAL ANALYSIS about Fig 3
## Conformity Exponent
###################################################################
parameters_stan_data$conformity_exponent = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$theta

  ## FITTING
conformity_exponent_stan_model = stan_model(file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/conformity_exponent.stan') # debug
conformity_exponent_fitting = sampling(
  conformity_exponent_stan_model, data=parameters_stan_data, seed=77,
  #control = list(adapt_delta = 0.90, max_treedepth=13),
  pars=c('beta','gamma','sigma_e','sigma_r','r_group','conformity_exponent_pred_25','conformity_exponent_pred_50','conformity_exponent_pred_78'),
  init=function() {
    list(
      sigma_e=runif(1,0,3),
      sigma_r=runif(1,0,3)
      )
  },
  chains=8, iter=5000, warmup=1500, thin=2
)
conformity_exponent_fitting
saveRDS(conformity_exponent_fitting, file='mcmc_result/conformity_exponent_fitting.rds')

plot(conformity_exponent_fitting, pars=c('beta','gamma','sigma_e','sigma_r'))
plot(conformity_exponent_fitting, pars=c('conformity_exponent_pred_25'))

  ## 事後診断
traceplot(conformity_exponent_fitting, pars=c('beta','sigma_e'), inc_warmup=TRUE)



################################
## Supporting figure:
## Analysis only on Pos-freq-dep
################################
parameters_stan_data_posOnly$conformity_exponent = subset(performanceSummary4_stanData_posOnly, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$theta

  ## FITTING
conformity_exponent_stan_model = stan_model(file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/conformity_exponent.stan') # debug
conformity_exponent_fitting_posOnly = sampling(
  conformity_exponent_stan_model, data=parameters_stan_data_posOnly, seed=77,
  #control = list(adapt_delta = 0.90, max_treedepth=13),
  pars=c('beta','gamma','sigma_e','sigma_r','r_group','conformity_exponent_pred_25','conformity_exponent_pred_50','conformity_exponent_pred_78'),
  init=function() {
    list(
      sigma_e=runif(1,0,3),
      sigma_r=runif(1,0,3)
      )
  },
  chains=8, iter=5000, warmup=1500, thin=4
)
conformity_exponent_fitting_posOnly
saveRDS(conformity_exponent_fitting_posOnly, file='mcmc_result/conformity_exponent_fitting_posOnly.rds')
conformity_exponent_mcmc_result_posOnly <- rstan::extract(conformity_exponent_fitting_posOnly)
plot(conformity_exponent_fitting_posOnly, pars=c('beta','gamma','sigma_e','sigma_r'))
#plot(conformity_exponent_fitting_posOnly, pars=c('conformity_exponent_pred_25'))

  ## 事後診断
traceplot(conformity_exponent_fitting_posOnly, pars=c('beta','sigma_e','sigma_r'), inc_warmup=TRUE)

conformity_exponent_pred_data_posOnly = cbind(parameters_stan_data_posOnly$test_group_size, as.data.frame(t(apply(conformity_exponent_mcmc_result_posOnly$conformity_exponent_pred_25[,1:ncol(conformity_exponent_mcmc_result_posOnly$conformity_exponent_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
conformity_exponent_pred_data_posOnly = rbind(conformity_exponent_pred_data_posOnly, cbind(parameters_stan_data_posOnly$test_group_size, as.data.frame(t(apply(conformity_exponent_mcmc_result_posOnly$conformity_exponent_pred_50[,1:ncol(conformity_exponent_mcmc_result_posOnly$conformity_exponent_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
conformity_exponent_pred_data_posOnly = rbind(conformity_exponent_pred_data_posOnly, cbind(parameters_stan_data_posOnly$test_group_size, as.data.frame(t(apply(conformity_exponent_mcmc_result_posOnly$conformity_exponent_pred_78[,1:ncol(conformity_exponent_mcmc_result_posOnly$conformity_exponent_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
conformity_exponent_pred_data_posOnly$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(conformity_exponent_pred_data_posOnly) <- c("groupSize_standerized2",
              "conformity_exponent_p2.5", "conformity_exponent_p25", "conformity_exponent_p50", "conformity_exponent_p75", "conformity_exponent_p97.5","taskDifficulty")


## conformity exponent
FigureS_theta_posOnly_25 = ggplot() +
  geom_ribbon(data=subset(conformity_exponent_pred_data_posOnly,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=conformity_exponent_p25,ymax=conformity_exponent_p75), fill = "lightpink")+#fill = "lightpink"
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&frequency_dependence_4!='Single-condition'&frequency_dependence_4!='--'&taskDifficulty=='25%'),mapping=aes(groupSize_standerized2, theta, colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_line(data=subset(conformity_exponent_pred_data_posOnly,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=conformity_exponent_p50), linetype='dashed', size=0.6,color='red')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  labs(x='Group size\n (standardized)', y=bquote(atop('Conformity', 'exponent' ~ theta[i])), title='') +
  myTheme_inner()+ylim(c(-7,12))+xlim(c(-1.3,2.3))+
  theme(legend.position = c(0.02, 0.75))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
FigureS_theta_posOnly_50 = ggplot() +
  geom_ribbon(data=subset(conformity_exponent_pred_data_posOnly,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=conformity_exponent_p25,ymax=conformity_exponent_p75), fill = "lightpink")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='50%'),mapping=aes(groupSize_standerized2, theta, colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_line(data=subset(conformity_exponent_pred_data_posOnly,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=conformity_exponent_p50), linetype='dashed', size=0.6,color='red')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Conformity', 'exponent' ~ theta[i])), title='') +
  myTheme()+ylim(c(-7,12))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
FigureS_theta_posOnly_78 = ggplot() +
  geom_ribbon(data=subset(conformity_exponent_pred_data_posOnly,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, ymin=conformity_exponent_p25,ymax=conformity_exponent_p75), fill = "lightpink")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='78%'),mapping=aes(groupSize_standerized2, theta, colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_line(data=subset(conformity_exponent_pred_data_posOnly,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, y=conformity_exponent_p50), linetype='dashed', size=0.6,color='red')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Conformity', 'exponent' ~ theta[i])), title='') +
  myTheme()+ylim(c(-7,12))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

FigureS_theta_posOnly = plot_grid(FigureS_theta_posOnly_25, FigureS_theta_posOnly_50, FigureS_theta_posOnly_78, labels = c('','',''), ncol = 3, align = 'v')
#Figure3_c

## summary
FigureS_param_groupSize_posOnly = plot_grid(FigureS_copyingRate_posOnly, FigureS_theta_posOnly, labels=c("a","b"), ncol = 1, align = 'v')
#Figure3
ggsave(file = "FigureS_param_groupSize_posOnly.pdf", plot = FigureS_param_groupSize_posOnly, dpi = 700, width = 10, height = 6.67)
#ggsave(file = "FigureS_param_groupSize_posOnly.tiff", plot = FigureS_param_groupSize_posOnly, dpi = 700, width = 10, height = 6.67)
################################
## Supporting figure:
## Analysis only on Pos-freq-dep - END
################################


###################################################################
## Other parameters
## Learning Rate (alpha)
###################################################################
parameters_stan_data$alpha = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$alpha
  ## FITTING
alpha_stan_model = stan_model(file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/alpha.stan') # debug
alpha_fitting = sampling(
  alpha_stan_model, data=parameters_stan_data, seed=77,
  #control = list(adapt_delta = 0.90, max_treedepth=13),
  pars=c('beta','sigma_e','sigma_r','r_group','alpha_pred_25','alpha_pred_50','alpha_pred_78'),
  init=function() {
    list(
      sigma_e=runif(1,0,3),
      sigma_r=runif(1,0,3)
      )
  },
  chains=8, iter=5000, warmup=1500, thin=4
)
alpha_fitting
saveRDS(alpha_fitting, file='mcmc_result/alpha_fitting.rds')
alpha_mcmc_result <- rstan::extract(alpha_fitting)
plot(alpha_fitting, pars=c('beta','sigma_e','sigma_r'))

  ## 事後診断
traceplot(alpha_fitting, pars=c('beta','sigma_e','sigma_r'), inc_warmup=TRUE)

alpha_pred_data = cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(alpha_mcmc_result$alpha_pred_25[,1:ncol(alpha_mcmc_result$alpha_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
alpha_pred_data = rbind(alpha_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(alpha_mcmc_result$alpha_pred_50[,1:ncol(alpha_mcmc_result$alpha_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
alpha_pred_data = rbind(alpha_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(alpha_mcmc_result$alpha_pred_78[,1:ncol(alpha_mcmc_result$alpha_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
alpha_pred_data$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(alpha_pred_data) <- c("groupSize_standerized2",
              "alpha_p2.5", "alpha_p25", "alpha_p50", "alpha_p75", "alpha_p97.5","taskDifficulty")

## alpha
Figure_s_alpha_groupSize_25 = ggplot() +
  geom_ribbon(data=subset(alpha_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=alpha_p25,ymax=alpha_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='25%'),mapping=aes(groupSize_standerized2, alpha, colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_line(data=subset(alpha_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=alpha_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Learning', 'rate' ~ alpha[i])), title='Low Uncertainty') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
Figure_s_alpha_groupSize_50 = ggplot() +
  geom_ribbon(data=subset(alpha_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1.5), aes(x=groupSize_standerized2, ymin=alpha_p25,ymax=alpha_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='50%'),mapping=aes(groupSize_standerized2, alpha, colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_line(data=subset(alpha_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1.5), aes(x=groupSize_standerized2, y=alpha_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Learning', 'rate' ~ alpha[i])), title='Moderate Uncertainty') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
Figure_s_alpha_groupSize_78 = ggplot() +
  geom_ribbon(data=subset(alpha_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, ymin=alpha_p25,ymax=alpha_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='78%'),mapping=aes(groupSize_standerized2, alpha, colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_line(data=subset(alpha_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, y=alpha_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Learning', 'rate' ~ alpha[i])), title='High Uncertainty') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

Figure_s_alpha_groupSize = plot_grid(Figure_s_alpha_groupSize_25, Figure_s_alpha_groupSize_50, Figure_s_alpha_groupSize_78, labels = c('','',''), ncol = 3, align = 'v')

###################################################################
## Other parameters????
## average inverse temperature
###################################################################
parameters_stan_data$average_invTemp = subset(performanceSummary4_stanData, environment=='e1+e2'&frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition')$average_invTemp
  ## FITTING
average_invTemp_stan_model = stan_model(file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/average_invTemp.stan') # debug
average_invTemp_fitting = sampling(
  average_invTemp_stan_model, data=parameters_stan_data, seed=77,
  #control = list(adapt_delta = 0.90, max_treedepth=13),
  pars=c('beta','sigma_e','sigma_r','r_group','average_invTemp_pred_25','average_invTemp_pred_50','average_invTemp_pred_78'),
  init=function() {
    list(
      sigma_e=runif(1,0,3),
      sigma_r=runif(1,0,3)
      )
  },
  chains=8, iter=5000, warmup=1500, thin=4
)
average_invTemp_fitting
saveRDS(average_invTemp_fitting, file='mcmc_result/average_invTemp_fitting.rds')
average_invTemp_mcmc_result <- rstan::extract(average_invTemp_fitting)
plot(average_invTemp_fitting, pars=c('beta','sigma_e','sigma_r'))

  ## 事後診断
traceplot(average_invTemp_fitting, pars=c('beta','sigma_e','sigma_r'), inc_warmup=TRUE)

average_invTemp_pred_data = cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_invTemp_mcmc_result$average_invTemp_pred_25[,1:ncol(average_invTemp_mcmc_result$average_invTemp_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
average_invTemp_pred_data = rbind(average_invTemp_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_invTemp_mcmc_result$average_invTemp_pred_50[,1:ncol(average_invTemp_mcmc_result$average_invTemp_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
average_invTemp_pred_data = rbind(average_invTemp_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_invTemp_mcmc_result$average_invTemp_pred_78[,1:ncol(average_invTemp_mcmc_result$average_invTemp_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
average_invTemp_pred_data$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(average_invTemp_pred_data) <- c("groupSize_standerized2",
              "average_invTemp_p2.5", "average_invTemp_p25", "average_invTemp_p50", "average_invTemp_p75", "average_invTemp_p97.5","taskDifficulty")

## conformity exponent
Figure_s_average_invTemp_groupSize_25 = ggplot() +
  geom_ribbon(data=subset(average_invTemp_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=average_invTemp_p25,ymax=average_invTemp_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='25%'),mapping=aes(groupSize_standerized2, average_invTemp, colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_line(data=subset(average_invTemp_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=average_invTemp_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Average', 'inverse temperature' ~ bar(beta[i]))), title='') +
  myTheme()+ylim(c(-0.5,10))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
Figure_s_average_invTemp_groupSize_50 = ggplot() +
  geom_ribbon(data=subset(average_invTemp_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1.5), aes(x=groupSize_standerized2, ymin=average_invTemp_p25,ymax=average_invTemp_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='50%'),mapping=aes(groupSize_standerized2, average_invTemp, colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_line(data=subset(average_invTemp_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1.5), aes(x=groupSize_standerized2, y=average_invTemp_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Average', 'inverse temperature' ~ bar(beta[i]))), title='') +
  myTheme()+ylim(c(-0.5,10))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
Figure_s_average_invTemp_groupSize_78 = ggplot() +
  geom_ribbon(data=subset(average_invTemp_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, ymin=average_invTemp_p25,ymax=average_invTemp_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='78%'),mapping=aes(groupSize_standerized2, average_invTemp, colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_line(data=subset(average_invTemp_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, y=average_invTemp_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Average', 'inverse temperature' ~ bar(beta[i]))), title='') +
  myTheme()+ylim(c(-0.5,10))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

Figure_s_average_invTemp_groupSize = plot_grid(Figure_s_average_invTemp_groupSize_25, Figure_s_average_invTemp_groupSize_50, Figure_s_average_invTemp_groupSize_78, labels = c('','',''), ncol = 3, align = 'v')
## summary
Figure_s_otherParams_vs_groupSize = plot_grid(Figure_s_alpha_groupSize, Figure_s_average_invTemp_groupSize, labels=c("a","b"), ncol = 1, align = 'v')
ggsave(file = "Figure_s_otherParams_vs_groupSize.pdf", plot = Figure_s_otherParams_vs_groupSize, dpi = 700, width = 10, height = 6.67)
#ggsave(file = "Figure_s_otherParams_vs_groupSize.tiff", plot = Figure_s_otherParams_vs_groupSize, dpi = 700, width = 10, height = 6.67)





## changing beta
round_invTemp_posOnly = ggplot() +
                    geom_line(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'&frequency_dependence_4=='pos-freq-dep'), mapping=aes(round,beta,group=amazonID,colour=alpha),alpha=3/4)+
                    stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'&frequency_dependence_4=='pos-freq-dep'),aes(round,beta),fun.y = "median",geom="line",linetype = "dashed", colour = "black", size = 1)+
                    facet_grid(.~taskDifficulty2)+
                    labs(x='Round', y=expression(paste('Inverse temperature ', beta[i][t],sep="")), title='') +
                    xlim(c(1,70))+
                    #scale_colour_distiller(name=expression(paste('Learning\nrate',alpha[i],sep="")), palette = "RdYlBu", direction = -1)+
                    scale_color_viridis(name=expression(paste('Learning\nrate',alpha[i],sep="")), option="plasma")+
                    myTheme_legend()+
                    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

round_beta_fullModel = ggplot() +
                    geom_line(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'), mapping=aes(round,beta,group=amazonID,colour=frequency_dependence_4), alpha=2/4)+
                    stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,beta),fun.y = "median",geom="line", colour = "black", size = 1)+
                    facet_grid(taskDifficulty2~frequency_dependence_4)+
                    labs(x='Round', y='Inverse temperature', title='') +
                    #ylim(c(0,1))+
                    xlim(c(1,70))+
                    scale_colour_manual(values=c("random-choice"="#999999","AL"="#999999","Random-copying"="#999999",'neg-freq-dep'="blue", "pos-freq-dep"="red"), name="Strength of \nfrequency dependence") +
                    #scale_colour_distiller(name="Frequency \ndependence", palette = "Spectral", direction = -1)+
                    myTheme()+
                    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

ggsave(file = "round_beta_fullModel.pdf", plot = round_beta_fullModel, dpi = 600, width = 9*1.2, height = 9*1.2)




#######################################################################
##
## Supplementary figure -- simulation using other parameters 1
## Simulation Results
##
#######################################################################
#low_sigma = 0.3
#high_sigma = 0.9
#low_theta = 1
#high_theta = 3

## To get the tidy-data frame for the time series data
UNC_Moderate_otherParam1_timeSeries_sup1 = subset(UNC_Moderate_otherParam1_timeSeries_sup1, theta==low_theta|theta==high_theta)
UNC_Moderate_otherParam1_timeSeries_sup12 = data.frame(
  round = rep(1:lifetime_simulation, 2*2*3),
  groupSize = c(rep(3, lifetime_simulation*2*2), rep(10, lifetime_simulation*2*2), rep(30, lifetime_simulation*2*2)),
  groupSize_factor = c(rep('n = 3', lifetime_simulation*2*2), rep('n = 10', lifetime_simulation*2*2), rep('n = 30', lifetime_simulation*2*2)),
  theta = rep(c(rep(low_theta, lifetime_simulation*2), rep(high_theta,lifetime_simulation*2)), 3),
  #theta = rep(c(rep(1, lifetime_simulation*2), rep(3,lifetime_simulation*2)), 3),
  lambda = rep(rep(c(rep(low_sigma,lifetime_simulation), rep(high_sigma,lifetime_simulation)), 2), 3)
)
UNC_Moderate_otherParam1_timeSeries_sup12$groupSize_factor = factor(UNC_Moderate_otherParam1_timeSeries_sup12$groupSize_factor, levels=c('n = 3','n = 10','n = 30'))
UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate = NA
for(i in 1:nrow(UNC_Moderate_otherParam1_timeSeries_sup1)) {
  g = as.numeric(as.character(UNC_Moderate_otherParam1_timeSeries_sup1$groupSize[i]))
  the = as.numeric(as.character(UNC_Moderate_otherParam1_timeSeries_sup1$theta[i]))
  la = as.numeric(as.character(UNC_Moderate_otherParam1_timeSeries_sup1$lambda[i]))
  location = which(UNC_Moderate_otherParam1_timeSeries_sup12$groupSize==g&UNC_Moderate_otherParam1_timeSeries_sup12$theta==the&UNC_Moderate_otherParam1_timeSeries_sup12$lambda==la)
  if(length(location)>0) {
  UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate[location] = UNC_Moderate_otherParam1_timeSeries_sup1[i,4:(lifetime_simulation+3)]
  }
}
UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate = as.numeric(as.character(UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate))

for(t in 1:lifetime_simulation) {
  UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate[which(UNC_Moderate_otherParam1_timeSeries_sup12$groupSize==1&UNC_Moderate_otherParam1_timeSeries_sup12$round==t)] = mean(UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate[which(UNC_Moderate_otherParam1_timeSeries_sup12$groupSize==1&UNC_Moderate_otherParam1_timeSeries_sup12$round==t)])
}

Figure1_otherParam1_a =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, lambda==low_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize_factor, colour=groupSize_factor, linetype=groupSize_factor),size=1.5)+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.3\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_manual(name = '', values = c("black", "orange", "red"))+
  scale_linetype_manual(name = '', values = c("dotted", "dashed", "solid"))+
  #scale_color_distiller(name = '', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme_inner()+
  theme(legend.position = c(0.00, 0.23),
    legend.key.width = unit(15, "mm"))
Figure1_otherParam1_b =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==3 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==10 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==30 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.3\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme_inner()
Figure1_otherParam1_e =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==3 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==10 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==30 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme_inner()
Figure1_otherParam1_f =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==3 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==10 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==30 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme_inner()

Figure1_otherParam1 = plot_grid(Figure1_otherParam1_a, Figure1_otherParam1_b, Figure1_otherParam1_e, Figure1_otherParam1_f, labels = c('A','','',''), ncol = 2, align = 'v', label_size=13, label_fontfamily='Helvetica',label_fontface='bold',label_x=0.25,label_y=0.85,hjust=-0.5,vjust=-0.5)
ggsave(file = "Figure1_otherParam1_rev.pdf", plot = Figure1_otherParam1, dpi = 600, width = 7.5, height = 5)
#ggsave(file = "Figure1.tiff", plot = Figure1, dpi = 600, width = 9, height = 6)




## GGjoy plot
Figure2_otherParam1 = ggplot()+
  geom_joy(data= UNC_Moderate_otherParam1_eachGroup_sup1, aes(x = averageFirst, y = as.factor(lambda), linetype=as.factor(groupSize), colour=as.factor(groupSize),fill=as.factor(groupSize),alpha=as.factor(groupSize)), stat = "binline", bins = 40, scale = 0.98, draw_baseline = FALSE)+
  scale_color_manual(name = 'Group size:', values = c("black", "orange", "red"), labels=c("n = 3;", "n = 10;", "n = 30"))+
  scale_linetype_manual(name = 'Group size:', values = c("dotted", "dashed", "solid"), labels=c("n = 3;", "n = 10;", "n = 30"))+
  scale_fill_manual(name = 'Group size:', values = c(NA, "orange", "red"), labels=c("n = 3;", "n = 10;", "n = 30"))+
  scale_alpha_manual(name = 'Group size:', values = c(3/3, 2/3, 1/3), labels=c("n = 3;", "n = 10;", "n = 30"))+
  facet_grid(.~theta, labeller = label_bquote(cols = bar(italic(theta)) == .(theta)))+
  labs(x='Groups\' average performance', y=expression(paste('Social learning weight ',bar(sigma),sep="")),title='')+
  xlim(c(0,1))+
  theme(
    legend.position = 'top',
    axis.text.x = element_text(angle=90),
    axis.title.y = element_text(angle=90))+
  #scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  myTheme_legend()

ggsave(file = "Figure2_otherParam1.pdf", plot = Figure2_otherParam1, dpi = 600, width = 9, height = 6)
#ggsave(file = "Figure2_otherParam1.tiff", plot = Figure2_otherParam1, dpi = 600, width = 9, height = 6)




#######################################################################
##
## Supplementary figure -- simulation using other parameters 2
## Simulation Results
##
#######################################################################
## To get the tidy-data frame for the time series data
UNC_Moderate_otherParam2_timeSeries_sup1 = subset(UNC_Moderate_otherParam2_timeSeries_sup1, theta==low_theta|theta==high_theta)
UNC_Moderate_otherParam2_timeSeries_sup12 = data.frame(
  round = rep(1:lifetime_simulation, 2*2*3),
  groupSize = c(rep(3, lifetime_simulation*2*2), rep(10, lifetime_simulation*2*2), rep(30, lifetime_simulation*2*2)),
  groupSize_factor = c(rep('n = 3', lifetime_simulation*2*2), rep('n = 10', lifetime_simulation*2*2), rep('n = 30', lifetime_simulation*2*2)),
  theta = rep(c(rep(low_theta, lifetime_simulation*2), rep(high_theta,lifetime_simulation*2)), 3),
  #theta = rep(c(rep(1, lifetime_simulation*2), rep(3,lifetime_simulation*2)), 3),
  lambda = rep(rep(c(rep(low_sigma,lifetime_simulation), rep(high_sigma,lifetime_simulation)), 2), 3)
)
UNC_Moderate_otherParam2_timeSeries_sup12$groupSize_factor = factor(UNC_Moderate_otherParam2_timeSeries_sup12$groupSize_factor, levels=c('n = 3','n = 10','n = 30'))
UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate = NA
for(i in 1:nrow(UNC_Moderate_otherParam2_timeSeries_sup1)) {
  g = as.numeric(as.character(UNC_Moderate_otherParam2_timeSeries_sup1$groupSize[i]))
  the = as.numeric(as.character(UNC_Moderate_otherParam2_timeSeries_sup1$theta[i]))
  la = as.numeric(as.character(UNC_Moderate_otherParam2_timeSeries_sup1$lambda[i]))
  location = which(UNC_Moderate_otherParam2_timeSeries_sup12$groupSize==g&UNC_Moderate_otherParam2_timeSeries_sup12$theta==the&UNC_Moderate_otherParam2_timeSeries_sup12$lambda==la)
  if(length(location)>0) {
  UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate[location] = UNC_Moderate_otherParam2_timeSeries_sup1[i,4:(lifetime_simulation+3)]
  }
}
UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate = as.numeric(as.character(UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate))

for(t in 1:lifetime_simulation) {
  UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate[which(UNC_Moderate_otherParam2_timeSeries_sup12$groupSize==1&UNC_Moderate_otherParam2_timeSeries_sup12$round==t)] = mean(UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate[which(UNC_Moderate_otherParam2_timeSeries_sup12$groupSize==1&UNC_Moderate_otherParam2_timeSeries_sup12$round==t)])
}

Figure1_otherParam2_a =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, lambda==low_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize_factor, colour=groupSize_factor, linetype=groupSize_factor),size=1.5)+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.3\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_manual(name = '', values = c("black", "orange", "red"))+
  scale_linetype_manual(name = '', values = c("dotted", "dashed", "solid"))+
  #scale_color_distiller(name = '', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme_inner()+
  theme(legend.position = c(0.00, 0.23),
    legend.key.width = unit(15, "mm"))
Figure1_otherParam2_b =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==3 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==10 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==30 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.3\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme_inner()
Figure1_otherParam2_e =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==3 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==10 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==30 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme_inner()
Figure1_otherParam2_f =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==3 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==10 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==30 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 70, y = 0.93, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=5, family="Helvetica")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed", size=0.5)+
  scale_x_continuous(limits = c(0, lifetime_simulation), breaks=c(0,20,40,70,90)) +
  #geom_vline(xintercept = 70, linetype="solid", size=0.5)+
  myTheme_inner()

Figure1_otherParam2 = plot_grid(Figure1_otherParam2_a, Figure1_otherParam2_b, Figure1_otherParam2_e, Figure1_otherParam2_f, labels = c('B','','',''), ncol = 2, align = 'v', label_size=13, label_fontfamily='Helvetica',label_fontface='bold',label_x=0.25,label_y=0.85,hjust=-0.5,vjust=-0.5)
ggsave(file = "Figure1_otherParam2_rev.pdf", plot = Figure1_otherParam2, dpi = 600, width = 7.5, height = 5)
#ggsave(file = "Figure1.tiff", plot = Figure1, dpi = 600, width = 9, height = 6)

## merging
Figure1_otherParam_merged = plot_grid(Figure1_otherParam1,Figure1_otherParam2, labels = c('',''), ncol = 2, align = 'v')
ggsave(file = "Figure1_otherParam_merged.pdf", plot = Figure1_otherParam_merged, dpi = 600, width = 15, height = 5)

## GGjoy plot
Figure2_otherParam2 = ggplot()+
  geom_joy(data= UNC_Moderate_otherParam2_eachGroup_sup1, aes(x = averageFirst, y = as.factor(lambda), linetype=as.factor(groupSize), colour=as.factor(groupSize),fill=as.factor(groupSize),alpha=as.factor(groupSize)), stat = "binline", bins = 40, scale = 0.98, draw_baseline = FALSE)+
  scale_color_manual(name = 'Group size:', values = c("black", "orange", "red"), labels=c("n = 3;", "n = 10;", "n = 30"))+
  scale_linetype_manual(name = 'Group size:', values = c("dotted", "dashed", "solid"), labels=c("n = 3;", "n = 10;", "n = 30"))+
  scale_fill_manual(name = 'Group size:', values = c(NA, "orange", "red"), labels=c("n = 3;", "n = 10;", "n = 30"))+
  scale_alpha_manual(name = 'Group size:', values = c(3/3, 2/3, 1/3), labels=c("n = 3;", "n = 10;", "n = 30"))+
  facet_grid(.~theta, labeller = label_bquote(cols = bar(italic(theta)) == .(theta)))+
  labs(x='Groups\' average performance', y=expression(paste('Social learning weight ',bar(sigma),sep="")),title='')+
  xlim(c(0,1))+
  theme(
    legend.position = 'top',
    axis.text.x = element_text(angle=90),
    axis.title.y = element_text(angle=90))+
  #scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  myTheme_legend()

ggsave(file = "Figure2_otherParam2.pdf", plot = Figure2_otherParam2, dpi = 600, width = 9, height = 6)
#ggsave(file = "Figure2_otherParam2.tiff", plot = Figure2_otherParam2, dpi = 600, width = 9, height = 6)

### merge all
Figure2_otherParam_merge = plot_grid(Figure2_otherParam1,Figure2_otherParam2, labels = c('',''), ncol = 2, align = 'v')
FigureS_otherParam = plot_grid(Figure1_otherParam_merged,Figure2_otherParam_merge, labels = c('',''), ncol = 1, align = 'v')

ggsave(file = "FigureS_otherParam.pdf", plot = FigureS_otherParam, dpi = 600, width = 15, height = 11)


## Figure S2
FigureS2a = ggplot(data=subset(performanceSummary4, environment=='e1+e2'), aes(age)) +
  geom_histogram() +
  geom_vline(xintercept=32, linetype='dashed', colour='blue') +
  labs(x='Age', y='Count',title='')+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  myTheme()
performanceSummary4$sex_factor = performanceSummary4$male_female_other
performanceSummary4$sex_factor[which(performanceSummary4$sex_factor==1)] = 'Male'
performanceSummary4$sex_factor[which(performanceSummary4$sex_factor==2)] = 'Female'
performanceSummary4$sex_factor[which(performanceSummary4$sex_factor==3)] = 'Other'

FigureS2b = ggplot(data=subset(performanceSummary4, environment=='e1+e2'), aes(as.factor(sex_factor))) +
  geom_bar(stat='count') +
  labs(x='Gender', y='Count',title='')+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  myTheme()
FigureS2 = plot_grid(FigureS2a, FigureS2b, labels = c('a','b'), ncol = 2, align = 'v')
ggsave(file = "age_gender.pdf", plot = FigureS2, dpi = 600, width = 9, height = 3)








## Figure: Time evolution of the proportion of choosing the best option for each group
## Predicted Choice probability for each individual
## Using the best fit parameters and actual pay-off history for each subject,
## we can fit a most likely optimal-choice probability
softmax = function(Q, beta) {
  return( exp(beta*Q)/sum(exp(beta*Q)) )
}
frequencyDependentCopy = function(F, theta) {
  f = F + 1e-01
  return( f^theta/sum(f^theta) )
}
allBehaviouralData2$pred_p0 = NA
allBehaviouralData2$pred_p1 = NA
allBehaviouralData2$pred_p2 = NA
allBehaviouralData2$pred_p_optimal = NA
for(i in names(table(performanceSummary4$amazonID))) {
  thisSubject = subset(performanceSummary4, amazonID==i&environment=="environment 1")
  if(thisSubject$frequency_dependence_4 == "Single-condition") {
    ## Setting parameter value
    thisSubjectAlpha = thisSubject$alpha
    thisSubjectBeta = subset(netBeta_table_AL6_annealing_indiv, amazonID==i)$AL6_annealing_indiv_netBeta_p50
    thisSubjectPayoff = subset(allBehaviouralData2, amazonID==i)$result * 100
    thisSubjectChoice = subset(allBehaviouralData2, amazonID==i)$machine + 1
    thisSubjectStartT = min(subset(allBehaviouralData2, amazonID==i)$round) # Some participant missed the 1st round
    thisSubjectEndT = max(subset(allBehaviouralData2, amazonID==i)$round)
    thisSubjectRound = subset(allBehaviouralData2, amazonID==i)$round
    ## Setting other initial values
    thisSubjectQ = matrix(nrow=3, ncol=thisSubjectEndT)
    thisSubjectP = matrix(nrow=3, ncol=thisSubjectEndT)
    thisSubjectQ[,1] = 1e-10
    thisSubjectP[,1] = 1/3
    ## Calculating choice probability
    for(t in 2:thisSubjectEndT) {
      if((t-1) %in% thisSubjectRound) { # if this subject participated in the round t-1
        thisSubjectQ[,t] = thisSubjectQ[,t-1]
        ## Rescorla-Wagner rule
        thisSubjectQ[thisSubjectChoice[t-1],t] = (1-thisSubjectAlpha) * thisSubjectQ[thisSubjectChoice[t-1],t-1] + thisSubjectAlpha * thisSubjectPayoff[which(thisSubjectRound==t-1)]
      } else {
        thisSubjectQ[,t] = thisSubjectQ[,t-1]
      }
      if((t) %in% thisSubjectRound) { # if this subject participated in this round t
        thisSubjectP[,t] = softmax(thisSubjectQ[,t], thisSubjectBeta[t])
      }
    }
  } else if (thisSubject$frequency_dependence_4 != "--" & thisSubject$frequency_dependence_4 != "Single-condition") {
    ## Setting parameter value
    thisSubjectAlpha = thisSubject$alpha
    thisSubjectBeta = subset(netBeta_table_UNC6_sReduc_annealing, amazonID==i)$UNC6_sReduc_annealing_netBeta_p50
    thisSubjectTheta = thisSubject$theta
    thisSubjectSigma = subset(soc_table_UNC6_sReduc_annealing, amazonID==i)$UNC6_sReduc_annealing_soc_p50
    thisSubjectPayoff = subset(allBehaviouralData2, amazonID==i)$result * 100
    thisSubjectChoice = subset(allBehaviouralData2, amazonID==i)$machine + 1
    thisSubjectF = rbind(rbind(subset(allBehaviouralData2, amazonID==i)$socialFreq0, subset(allBehaviouralData2, amazonID==i)$socialFreq1),subset(allBehaviouralData2, amazonID==i)$socialFreq2)
    thisSubjectStartT = min(subset(allBehaviouralData2, amazonID==i)$round) # Some participant missed the 1st round
    thisSubjectEndT = max(subset(allBehaviouralData2, amazonID==i)$round)
    thisSubjectRound = subset(allBehaviouralData2, amazonID==i)$round
    ## Setting other initial values
    thisSubjectQ = matrix(nrow=3, ncol=thisSubjectEndT)
    thisSubjectP = matrix(nrow=3, ncol=thisSubjectEndT)
    thisSubjectQ[,1] = 1e-10
    thisSubjectP[,1] = 1/3
    ## Calculating choice probability
    for(t in 2:thisSubjectEndT) {
      if((t-1) %in% thisSubjectRound) { # if this subject participated in the round t-1
        thisSubjectQ[,t] = thisSubjectQ[,t-1]
        ## Rescorla-Wagner rule
        thisSubjectQ[thisSubjectChoice[t-1],t] = (1-thisSubjectAlpha) * thisSubjectQ[thisSubjectChoice[t-1],t-1] + thisSubjectAlpha * thisSubjectPayoff[which(thisSubjectRound==t-1)]
      } else {
        thisSubjectQ[,t] = thisSubjectQ[,t-1]
      }
      # subtracting this subject's own choice from social frequency information
      thisSubjectThisRoundF = thisSubjectF[,which(thisSubjectRound==t)]
      if((t-1) %in% thisSubjectRound) { thisSubjectThisRoundF[thisSubjectChoice[t-1]] = thisSubjectThisRoundF[thisSubjectChoice[t-1]] - 1}
      ## Calculating choice probability only if this subject chose at this round t
      if((t) %in% thisSubjectRound) {
        thisSubjectP[,t] = (1-thisSubjectSigma[t]) * softmax(thisSubjectQ[,t], thisSubjectBeta[t]) + thisSubjectSigma[t] * frequencyDependentCopy(thisSubjectThisRoundF, thisSubjectTheta)
      }
    }
  } else {
    thisSubjectStartT = min(subset(allBehaviouralData2, amazonID==i)$round) # Some participant missed the 1st round
    thisSubjectEndT = max(subset(allBehaviouralData2, amazonID==i)$round)
    thisSubjectRound = subset(allBehaviouralData2, amazonID==i)$round
    thisSubjectQ = matrix(nrow=3, ncol=thisSubjectEndT)
    thisSubjectP = matrix(nrow=3, ncol=thisSubjectEndT)
  }
  ## Insert the result
  allBehaviouralData2$pred_p0[which(allBehaviouralData2$amazonID==i)] = thisSubjectP[1,thisSubjectRound]
  allBehaviouralData2$pred_p1[which(allBehaviouralData2$amazonID==i)] = thisSubjectP[2,thisSubjectRound]
  allBehaviouralData2$pred_p2[which(allBehaviouralData2$amazonID==i)] = thisSubjectP[3,thisSubjectRound]
}

allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==1&allBehaviouralData2$round<=40)] = allBehaviouralData2$pred_p2[which(allBehaviouralData2$condition==1&allBehaviouralData2$round<=40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==2&allBehaviouralData2$round<=40)] = allBehaviouralData2$pred_p2[which(allBehaviouralData2$condition==2&allBehaviouralData2$round<=40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==3&allBehaviouralData2$round<=40)] = allBehaviouralData2$pred_p1[which(allBehaviouralData2$condition==3&allBehaviouralData2$round<=40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==4&allBehaviouralData2$round<=40)] = allBehaviouralData2$pred_p1[which(allBehaviouralData2$condition==4&allBehaviouralData2$round<=40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==5&allBehaviouralData2$round<=40)] = allBehaviouralData2$pred_p0[which(allBehaviouralData2$condition==5&allBehaviouralData2$round<=40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==6&allBehaviouralData2$round<=40)] = allBehaviouralData2$pred_p0[which(allBehaviouralData2$condition==6&allBehaviouralData2$round<=40)]

allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==1&allBehaviouralData2$round> 40)] = allBehaviouralData2$pred_p0[which(allBehaviouralData2$condition==1&allBehaviouralData2$round> 40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==2&allBehaviouralData2$round> 40)] = allBehaviouralData2$pred_p1[which(allBehaviouralData2$condition==2&allBehaviouralData2$round> 40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==3&allBehaviouralData2$round> 40)] = allBehaviouralData2$pred_p0[which(allBehaviouralData2$condition==3&allBehaviouralData2$round> 40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==4&allBehaviouralData2$round> 40)] = allBehaviouralData2$pred_p2[which(allBehaviouralData2$condition==4&allBehaviouralData2$round> 40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==5&allBehaviouralData2$round> 40)] = allBehaviouralData2$pred_p1[which(allBehaviouralData2$condition==5&allBehaviouralData2$round> 40)]
allBehaviouralData2$pred_p_optimal[which(allBehaviouralData2$condition==6&allBehaviouralData2$round> 40)] = allBehaviouralData2$pred_p2[which(allBehaviouralData2$condition==6&allBehaviouralData2$round> 40)]

## Black lines indicate means for each group of group condition (i.e. n >= 2) or for all session of solitary condition (n = 1)
## Red lines indicate the fitted model's prediction
## Light and dark shades indicate the 95% and 50% quantile of posthoc simulation, respectively. Note that the posthoc simulation was run for 90 rounds.
FigureS_eachGroupsPerformance = ggplot() +
  geom_ribbon(data=posthoc_sim_summary_all %>% dplyr::filter(value=="optimalChoiceProb"), mapping=aes(round, ymin=X2.5., ymax=X97.5., group=groupSize), alpha = 1/6) +
  geom_ribbon(data=posthoc_sim_summary_all %>% dplyr::filter(value=="optimalChoiceProb"), mapping=aes(round, ymin=X25., ymax=X75., group=groupSize), alpha = 3/6) +
  stat_summary(data=allBehaviouralData2, mapping=aes(round,optimalChoice, group=room_indiv_group), fun.y = mean, geom="line")+
  stat_summary(data=allBehaviouralData2, mapping=aes(round,pred_p_optimal, group=room_indiv_group), fun.y = mean, geom="line", colour="red")+
  facet_grid(sizeCategory_string ~ taskDifficulty3) +
  scale_y_continuous(limits = c(-0.05,1.05), breaks=c(0.0,1.0)) +
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  geom_hline(yintercept=1/3, linetype="dashed")+
  myTheme()+
  theme(legend.position='none',
  axis.text.y = element_text(size=11, family="Helvetica" ,colour='black'),
  strip.text.y = element_text(angle = 0))+
  NULL

ggsave(file = "data_and_analysis/FigureS_eachGroupsPerformance.pdf", plot = FigureS_eachGroupsPerformance, dpi = 600, width = 9, height = 9)



#################################################################################################
###
### Supplementary statistical analysis
###
### Performance comparison using *EARNINGS* as a measure of choice quality
### -- With effect of group size and social learning
###
### Figure S?? (revision)
#################################################################################################
## earnings (payoff)
onlineExp_stanData_25$payoff = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='25%') %>% dplyr::pull(result)*100
onlineExp_stanData_50$payoff = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='50%') %>% dplyr::pull(result)*100
onlineExp_stanData_78$payoff = allBehaviouralData2 %>% dplyr::filter(taskDifficulty=='78%') %>% dplyr::pull(result)*100


## Compiling the model
stanmodel_onlineExp_supp1 = stan_model(file='stan_model6/timeSeriesDiffTest_optimRate_groupSize_supp1.stan')



## _25%
## Fitting
onlineExp_fitting_groupSize_supp1_25 = sampling(
  stanmodel_onlineExp_supp1, data=onlineExp_stanData_25, seed=72,
  #control = list(adapt_delta = 0.80),
  pars=c('sizeEffect','s_i','s_p','mu0','di0','mu','di','p_indiv','p_group','p_group_small','p_group_large','p_group_veryLarge'),
  init=function() {
    list(mu0=runif(1,0,4))
  },
  chains=8, iter=2000, warmup=1000, thin=5
)
onlineExp_fitting_groupSize_supp1_25

saveRDS(onlineExp_fitting_groupSize_supp1_25, file='mcmc_result/onlineExp_fitting_groupSize_supp1_25.rds')

trace_onlineExp_fitting_groupSize_supp1_25 = traceplot(onlineExp_fitting_groupSize_supp1_25, pars=c('sizeEffect', 's_i', 'mu0', 'di0'), inc_warmup=TRUE)
ggsave(file = "mcmc_result/trace_onlineExp_fitting_groupSize_supp1_25.pdf", plot = trace_onlineExp_fitting_groupSize_supp1_25, dpi = 600, width = 9, height = 9)

#plot(onlineExp_fitting_groupSize_supp1_25, pars=c('p_indiv','p_group_small'))
plot(onlineExp_fitting_groupSize_supp1_25, pars=c('mu','di'))
plot(onlineExp_fitting_groupSize_supp1_25, pars=c('sizeEffect','s_i','s_p'))

data_onlineExp_fitting_groupSize_supp1_25 <- rstan::extract(onlineExp_fitting_groupSize_supp1_25)
fittedValues_onlineExp_fitting_groupSize_supp1_25 =
  as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_25$mu[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_25$di[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_25$p_indiv[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_25$p_group_small[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_25$p_group_large[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_25$p_group_veryLarge[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_25$p_group[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
## rename
colnames(fittedValues_onlineExp_fitting_groupSize_supp1_25) = c("mu_p2.5", "mu_p25", "mu_p50", "mu_p75", "mu_p97.5","di_p2.5", "di_p25", "di_p50", "di_p75", "di_p97.5","p_indiv_p2.5", "p_indiv_p25", "p_indiv_p50", "p_indiv_p75", "p_indiv_p97.5","p_group_small_p2.5", "p_group_small_p25", "p_group_small_p50", "p_group_small_p75", "p_group_small_p97.5","p_group_large_p2.5", "p_group_large_p25", "p_group_large_p50", "p_group_large_p75", "p_group_large_p97.5","p_group_veryLarge_p2.5", "p_group_veryLarge_p25", "p_group_veryLarge_p50", "p_group_veryLarge_p75", "p_group_veryLarge_p97.5","p_group_p2.5", "p_group_p25", "p_group_p50", "p_group_p75", "p_group_p97.5")
fittedValues_onlineExp_fitting_groupSize_supp1_25$taskDifficulty = "25%"

## _50%
## Fitting
onlineExp_fitting_groupSize_supp1_50 = sampling(
  stanmodel_onlineExp_supp1, data=onlineExp_stanData_50, seed=72,
  #control = list(adapt_delta = 0.80),
  pars=c('sizeEffect','s_i','s_p','mu0','di0','mu','di','p_indiv','p_group','p_group_small','p_group_large','p_group_veryLarge'),
  init=function() {
    list(mu0=runif(1,0,4))
  },
  chains=8, iter=2000, warmup=1000, thin=5
)
onlineExp_fitting_groupSize_supp1_50

saveRDS(onlineExp_fitting_groupSize_supp1_50, file='mcmc_result/onlineExp_fitting_groupSize_supp1_50.rds')

trace_onlineExp_fitting_groupSize_supp1_50 = traceplot(onlineExp_fitting_groupSize_supp1_50, pars=c('sizeEffect', 's_i', 'mu0', 'di0'), inc_warmup=TRUE)
ggsave(file = "mcmc_result/trace_onlineExp_fitting_groupSize_supp1_50.pdf", plot = trace_onlineExp_fitting_groupSize_supp1_50, dpi = 600, width = 9, height = 9)

#plot(onlineExp_fitting_groupSize_supp1_50, pars=c('p_indiv','p_group_small'))
plot(onlineExp_fitting_groupSize_supp1_50, pars=c('mu','di'))
plot(onlineExp_fitting_groupSize_supp1_50, pars=c('sizeEffect','s_i','s_p'))

data_onlineExp_fitting_groupSize_supp1_50 <- rstan::extract(onlineExp_fitting_groupSize_supp1_50)
fittedValues_onlineExp_fitting_groupSize_supp1_50 =
  as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_50$mu[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_50$di[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_50$p_indiv[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_50$p_group_small[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_50$p_group_large[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_50$p_group_veryLarge[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_50$p_group[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
## rename
colnames(fittedValues_onlineExp_fitting_groupSize_supp1_50) = c("mu_p2.5", "mu_p25", "mu_p50", "mu_p75", "mu_p97.5","di_p2.5", "di_p25", "di_p50", "di_p75", "di_p97.5","p_indiv_p2.5", "p_indiv_p25", "p_indiv_p50", "p_indiv_p75", "p_indiv_p97.5","p_group_small_p2.5", "p_group_small_p25", "p_group_small_p50", "p_group_small_p75", "p_group_small_p97.5","p_group_large_p2.5", "p_group_large_p25", "p_group_large_p50", "p_group_large_p75", "p_group_large_p97.5","p_group_veryLarge_p2.5", "p_group_veryLarge_p25", "p_group_veryLarge_p50", "p_group_veryLarge_p75", "p_group_veryLarge_p97.5","p_group_p2.5", "p_group_p25", "p_group_p50", "p_group_p75", "p_group_p97.5")
fittedValues_onlineExp_fitting_groupSize_supp1_50$taskDifficulty = "25%"


## _78%
## Fitting
onlineExp_fitting_groupSize_supp1_78 = sampling(
  stanmodel_onlineExp_supp1, data=onlineExp_stanData_78, seed=72,
  #control = list(adapt_delta = 0.80),
  pars=c('sizeEffect','s_i','s_p','mu0','di0','mu','di','p_indiv','p_group','p_group_small','p_group_large','p_group_veryLarge'),
  init=function() {
    list(mu0=runif(1,0,4))
  },
  chains=8, iter=2000, warmup=1000, thin=5
)
onlineExp_fitting_groupSize_supp1_78

saveRDS(onlineExp_fitting_groupSize_supp1_78, file='mcmc_result/onlineExp_fitting_groupSize_supp1_78.rds')

trace_onlineExp_fitting_groupSize_supp1_78 = traceplot(onlineExp_fitting_groupSize_supp1_78, pars=c('sizeEffect', 's_i', 'mu0', 'di0'), inc_warmup=TRUE)
ggsave(file = "mcmc_result/trace_onlineExp_fitting_groupSize_supp1_78.pdf", plot = trace_onlineExp_fitting_groupSize_supp1_78, dpi = 600, width = 9, height = 9)

#plot(onlineExp_fitting_groupSize_supp1_78, pars=c('p_indiv','p_group_small'))
plot(onlineExp_fitting_groupSize_supp1_78, pars=c('mu','di'))
plot(onlineExp_fitting_groupSize_supp1_78, pars=c('sizeEffect','s_i','s_p'))

data_onlineExp_fitting_groupSize_supp1_78 <- rstan::extract(onlineExp_fitting_groupSize_supp1_78)
fittedValues_onlineExp_fitting_groupSize_supp1_78 =
  as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_78$mu[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_78$di[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_78$p_indiv[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_78$p_group_small[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_78$p_group_large[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_78$p_group_veryLarge[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) %>%
  cbind(as.data.frame(t(apply(data_onlineExp_fitting_groupSize_supp1_78$p_group[,1:70], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
## rename
colnames(fittedValues_onlineExp_fitting_groupSize_supp1_78) = c("mu_p2.5", "mu_p25", "mu_p50", "mu_p75", "mu_p97.5","di_p2.5", "di_p25", "di_p50", "di_p75", "di_p97.5","p_indiv_p2.5", "p_indiv_p25", "p_indiv_p50", "p_indiv_p75", "p_indiv_p97.5","p_group_small_p2.5", "p_group_small_p25", "p_group_small_p50", "p_group_small_p75", "p_group_small_p97.5","p_group_large_p2.5", "p_group_large_p25", "p_group_large_p50", "p_group_large_p75", "p_group_large_p97.5","p_group_veryLarge_p2.5", "p_group_veryLarge_p25", "p_group_veryLarge_p50", "p_group_veryLarge_p75", "p_group_veryLarge_p97.5","p_group_p2.5", "p_group_p25", "p_group_p50", "p_group_p75", "p_group_p97.5")
fittedValues_onlineExp_fitting_groupSize_supp1_78$taskDifficulty = "78%"



####
## Plot
####

fittedValues_onlineExp_fitting_groupSize_supp1_25$taskDifficulty3 = "Low Uncertainty"
fittedValues_onlineExp_fitting_groupSize_supp1_50$taskDifficulty3 = "Moderate Uncertainty"
fittedValues_onlineExp_fitting_groupSize_supp1_78$taskDifficulty3 = "High Uncertainty"

fittedValues_onlineExp_fitting_groupSize_supp1 = rbind(fittedValues_onlineExp_fitting_groupSize_supp1_25, fittedValues_onlineExp_fitting_groupSize_supp1_50) %>% rbind(fittedValues_onlineExp_fitting_groupSize_supp1_78)
fittedValues_onlineExp_fitting_groupSize_supp1$round = rep(1:70, 3)
fittedValues_onlineExp_fitting_groupSize_supp1$taskDifficulty3 = factor(fittedValues_onlineExp_fitting_groupSize_supp1$taskDifficulty3, levels=c("Low Uncertainty","Moderate Uncertainty","High Uncertainty"))

fittedValues_onlineExp_fitting_groupSize_supp1$significant_indiv_vs_group = 0
fittedValues_onlineExp_fitting_groupSize_supp1$significant_indiv_vs_group[which(fittedValues_onlineExp_fitting_groupSize_supp1$di_p2.5>0)] = 1

group_level_summary_supp1_small_25 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=9&taskDifficulty=='25%')$result*100, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=9&taskDifficulty=='25%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=9&taskDifficulty=='25%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_supp1_large_25 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>9&taskDifficulty=='25%')$result*100, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>9&taskDifficulty=='25%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory>9&taskDifficulty=='25%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_supp1_small_50 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=6&taskDifficulty=='50%')$result*100, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=6&taskDifficulty=='50%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=6&taskDifficulty=='50%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_supp1_large_50 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>6&taskDifficulty=='50%')$result*100, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>6&taskDifficulty=='50%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory>6&taskDifficulty=='50%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_supp1_small_78 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=11&taskDifficulty=='78%')$result*100, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=11&taskDifficulty=='78%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory<=11&taskDifficulty=='78%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)
group_level_summary_supp1_large_78 = tapply(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>11&taskDifficulty=='78%')$result*100, list(subset(allBehaviouralData2, sizeCategory>1&sizeCategory>11&taskDifficulty=='78%')$round,subset(allBehaviouralData2, sizeCategory>1&sizeCategory>11&taskDifficulty=='78%')$room_indiv_group), mean, na.rm=TRUE) %>% apply(1,mean, na.rm=TRUE)

group_level_summary_supp1 = data.frame(
  round = rep(1:70, 3),
  taskDifficulty3 = factor(c(rep('Low Uncertainty', 70), rep('Moderate Uncertainty',70), rep('High Uncertainty',70)),levels=c("Low Uncertainty","Moderate Uncertainty","High Uncertainty")),
  small = c(group_level_summary_supp1_small_25, group_level_summary_supp1_small_50, group_level_summary_supp1_small_78),
  large = c(group_level_summary_supp1_large_25, group_level_summary_supp1_large_50, group_level_summary_supp1_large_78)
  )

## plot indiv_group
timeseries_payoff_analisys_groupSize_plot = ggplot() +
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, ymin=p_indiv_p2.5, ymax=p_indiv_p97.5),fill='black', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, ymin=p_indiv_p25, ymax=p_indiv_p75),fill='black', alpha=3/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, ymin=p_group_small_p2.5, ymax=p_group_small_p97.5),fill='orange', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, ymin=p_group_small_p25, ymax=p_group_small_p75),fill='orange', alpha=3/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, ymin=p_group_large_p2.5, ymax=p_group_large_p97.5),fill='red', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, ymin=p_group_large_p25, ymax=p_group_large_p75),fill='red', alpha=3/6)+
  stat_summary(data=subset(allBehaviouralData2, indiv_group==0), mapping=aes(x=round,result*100), fun.y = mean, geom="line",colour='black',size=1,alpha=3/4)+
  geom_line(data=group_level_summary_supp1, mapping=aes(round, small), colour='orange', size=1, alpha=3/4)+
  geom_line(data=group_level_summary_supp1, mapping=aes(round, large), colour='red', size=1, alpha=3/4)+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, y=p_indiv_p50),colour='black',linetype='dashed')+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, y=p_group_small_p50),colour='orange',linetype='dashed')+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, y=p_group_large_p50),colour='red',linetype='dashed')+
  facet_grid(.~taskDifficulty3)+
  #scale_y_continuous(limits = c(-0.05,1.05), breaks=c(0.0,0.25,0.50,0.75,1.0)) +
  labs(x='Rounds', y='Average earned payoff\n(US cents)')+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  #geom_hline(yintercept=1/3, linetype="dashed")+
  myTheme_legend()+
  NULL

ggsave(file = "data_and_analysis/timeseries_payoff_analisys_groupSize_plot.pdf", plot = timeseries_payoff_analisys_groupSize_plot, dpi = 600, width = 9, height = 3)


difference_btwn_indiv_vs_group_plot_supp1 = ggplot() +
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, ymin=di_p2.5, ymax=di_p97.5),fill='black', alpha=1/6)+
  geom_ribbon(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, ymin=di_p25, ymax=di_p75),fill='black', alpha=3/6)+
  geom_line(data=fittedValues_onlineExp_fitting_groupSize_supp1, mapping=aes(x=round, y=di_p50),colour='black')+
  facet_grid(.~taskDifficulty3)+
  labs(x='Rounds', y=expression(paste('Differences between\nsolitaries and groups ',xi[t],sep="")))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  geom_hline(yintercept=0, linetype="dashed")+
  myTheme_legend()+
  NULL


Figure_payoff_comparison = plot_grid(timeseries_payoff_analisys_groupSize_plot, difference_btwn_indiv_vs_group_plot_supp1, labels = c('',''), ncol = 1, align = 'v')
ggsave(file = "data_and_analysis/mean_payoff_comparison.pdf", plot = Figure_payoff_comparison, dpi = 600, width = 9, height = 5)



## Reviewer 2's comment 9
## what is the normalized distribution of social learning rates in the three group size and 3 task uncertainties after 40 rounds?
meanCopyingRate_at_2 =
  ggplot(performanceSummary4 %>% dplyr::filter(frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition'&environment=='e1+e2'))+
  geom_point(aes(groupSize, soc2,colour=frequency_dependence_4, shape=frequency_dependence_4),alpha=3/4)+
  facet_grid(. ~ taskDifficulty3)+
  labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='Snapshots of social learning weight\n 1st environment') +
  myTheme()+ylim(c(0,1))+#xlim(c(-1.3,2.3))+
  geom_hline(yintercept=0,size=0.5,linetype='dotted')+
  geom_hline(yintercept=1,size=0.5,linetype='dotted')+
  theme(legend.position = c(-0.00, 0.80), legend.text = element_text(size=9, family="Helvetica" ))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep."))+
  NULL
meanCopyingRate_at_70 =
  ggplot(performanceSummary4 %>% dplyr::filter(frequency_dependence_4!='--'&frequency_dependence_4!='Single-condition'&environment=='e1+e2'))+
  geom_point(aes(groupSize, soc70,colour=frequency_dependence_4, shape=frequency_dependence_4),alpha=3/4)+
  facet_grid(. ~ taskDifficulty3)+
  labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='Snapshots of social learning weight\n 2nd environment') +
  myTheme()+ylim(c(0,1))+#xlim(c(-1.3,2.3))+
  geom_hline(yintercept=0,size=0.5,linetype='dotted')+
  geom_hline(yintercept=1,size=0.5,linetype='dotted')+
  theme(legend.position = c(-0.00, 0.80), legend.text = element_text(size=9, family="Helvetica" ))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep."))+
  NULL


meanCopyingRate_snapshot = plot_grid(meanCopyingRate_at_2, meanCopyingRate_at_70, labels = c('',''), ncol = 2, align = 'v')
ggsave(file = "data_and_analysis/meanCopyingRate_snapshot.pdf", plot = meanCopyingRate_snapshot, dpi = 600, width = 18, height = 3)


## Reviewer 3 comment no 6
#I suspect that this may result from the "assumption" that players are competing with each other for a limited return.

competitive_vs_strategies = performanceSummary4 %>%
    dplyr::filter(environment=='e1+e2') %>%
    dplyr::select(amazonID,frequency_dependence_4) %>%
    dplyr::inner_join(questionnaire_data_exp, by='amazonID')









