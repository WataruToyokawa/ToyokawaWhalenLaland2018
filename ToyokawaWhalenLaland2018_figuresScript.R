##########################################################
## Figures and GLMM analyses in
## Social learning strategies regulate the wisdom and madness of interactive online crowds
## Wataru Toyokawa*, Andrew Whalen and Kevin Laland
## University of St Andrews
## *corresponding author: wt25@st-andrews.ac.uk
##########################################################
# README
# The following script will generate Figure 1, 2, 4, and 3 in this order. 
# You can re-generate both all the figures in the main text and some in the supplementary material
# by just doing copy-and-paste.
# 

# Functions
myTheme = function() {
  theme(
    legend.position = 'none',
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Times" ), #"Times"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=15, family="Times" ,colour='black'),
    axis.text.y = element_text(size=15, family="Times" ,colour='black'),
    axis.title.x=element_text(size=16, family="Times" ),
    axis.title.y=element_text(size=16, family="Times" ),
    axis.line =element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    #panel.grid = element_blank(),
    panel.grid = element_line(colour = "black", size = 0.8),
    plot.title = element_text(size=15, family="Times" ),
    plot.background = element_rect(colour = "white")
  )
}
myTheme_small = function() {
    theme(
      legend.position = 'none',
      strip.background = element_rect(fill = NA),
      panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
      strip.text = element_text(size=10, family="Times" ), #"Times"
      #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
      #axis.text.x = element_text(size=9, family="Times" ,colour='black'),
      axis.text.x = element_blank(),
      #axis.text.y = element_text(size=9, family="Times" ,colour='black'),
      axis.text.y = element_blank(),
      axis.title.x=element_text(size=10, family="Times" ),
      axis.title.y=element_text(size=10, family="Times" ),
      axis.line =element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      #panel.grid = element_blank(),
      panel.grid = element_line(colour = "black", size = 0.8),
      plot.title = element_text(size=9, family="Times" ),
      plot.background = element_rect(colour = "white")
    )
  }
myTheme_legend = function() {
  theme(
    #legend.position = 'right',
    legend.title = element_text(size=16, family="Times" ,colour='black'),
    legend.text = element_text(size=15, family="Times" ,colour='black'),
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Times" ), #"Times"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=15, family="Times" ,colour='black'),
    axis.text.y = element_text(size=15, family="Times" ,colour='black'),
    axis.title.x=element_text(size=16, family="Times" ),
    axis.title.y=element_text(size=16, family="Times" ),
    panel.background = element_rect(fill = "white", colour = NA),
    #panel.grid = element_blank(),
    panel.grid = element_line(colour = "black", size = 0.8),
    plot.title = element_text(size=15, family="Times" ),
    plot.background = element_rect(colour = "white")
  )
}
myTheme_inner = function() {
  theme(
    #legend.position = 'none',
    legend.title = element_text(size=13, family="Times" ,colour='black'),
    legend.text = element_text(size=13, family="Times" ,colour='black'),
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 1,  linetype="solid"),
    strip.text = element_text(size=16, family="Times" ), #"Times"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=15, family="Times" ,colour='black'),
    axis.text.y = element_text(size=15, family="Times" ,colour='black'),
    axis.title.x=element_text(size=16, family="Times" ),
    axis.title.y=element_text(size=16, family="Times" ),
    axis.line =element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_line(colour = "black", size = 0.8),
    plot.title = element_text(size=15, family="Times" ),
    plot.background = element_rect(colour = "white")
  )
}
# Loading stuff
library(readr)
library(cowplot)
library(ggjoy)
library(rstan)
library(ggmcmc)
library(viridisLite)
library(viridis)
## parallel computing
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
# Reading data
UNC_Moderate_expParameters_eachGroup <- read_csv("UNC_Moderate_expParameters_eachGroup.csv")
UNC_Moderate_expParameters_timeSeries0 <- read_csv("UNC_Moderate_expParameters_timeSeries.csv")
allBehaviouralData2 <- read_csv("allBehaviouralData2.csv")
simulation_fullModel_25 <- read_csv("simulation_fullModel_25.csv")
simulation_fullModel_50 <- read_csv("simulation_fullModel_50.csv")
simulation_fullModel_78 <- read_csv("simulation_fullModel_78.csv")
performanceSummary4 <- read_csv("performanceSummary4.csv")
temp_UNC6_sReduc_annealing <- read_csv("temp_UNC6_sReduc_annealing.csv")
temp_AL6_annealing_indiv <- read_csv("temp_AL6_annealing_indiv.csv")
temp_UNC6_sReduc_annealing_simulation <- read_csv("temp_UNC6_sReduc_annealing_simulation.csv")
soc_table_UNC6_sReduc_annealing <- read_csv("soc_table_UNC6_sReduc_annealing.csv")
netBeta_table_UNC6_sReduc_annealing <- read_csv("netBeta_table_UNC6_sReduc_annealing.csv")
netBeta_table_AL6_annealing_indiv <- read_csv("netBeta_table_AL6_annealing_indiv.csv")
UNC_Moderate_otherParam1_eachGroup_sup1 <- read_csv("UNC_Moderate_otherParam1_eachGroup_sup1.csv")
UNC_Moderate_otherParam1_timeSeries_sup1 <- read_csv("UNC_Moderate_otherParam1_timeSeries_sup1.csv")
UNC_Moderate_otherParam2_eachGroup_sup1 <- read_csv("UNC_Moderate_otherParam2_eachGroup_sup1.csv")
UNC_Moderate_otherParam2_timeSeries_sup1 <- read_csv("UNC_Moderate_otherParam2_timeSeries_sup1.csv")

#######################################################################
##
## Figure 1 and Figure 2
## Simulation Results
##
#######################################################################
low_sigma = 0.4
high_sigma = 0.9
low_theta = 1
high_theta = 3

## To get the tidy-data frame for the time series data
UNC_Moderate_expParameters_timeSeries = subset(UNC_Moderate_expParameters_timeSeries0, theta==low_theta|theta==high_theta)
UNC_Moderate_expParameters_timeSeries2 = data.frame(
  round = rep(1:70, 2*2*3),
  groupSize = c(rep(3, 70*2*2), rep(10, 70*2*2), rep(30, 70*2*2)),
  groupSize_factor = c(rep('n = 3', 70*2*2), rep('n = 10', 70*2*2), rep('n = 30', 70*2*2)),
  theta = rep(c(rep(low_theta, 70*2), rep(high_theta,70*2)), 3),
  #theta = rep(c(rep(1, 70*2), rep(3,70*2)), 3),
  lambda = rep(rep(c(rep(low_sigma,70), rep(high_sigma,70)), 2), 3)
)
UNC_Moderate_expParameters_timeSeries2$groupSize_factor = factor(UNC_Moderate_expParameters_timeSeries2$groupSize_factor, levels=c('n = 3','n = 10','n = 30'))
UNC_Moderate_expParameters_timeSeries2$optimChoiceRate = NA
for(i in 1:nrow(UNC_Moderate_expParameters_timeSeries)) {
  g = as.numeric(as.character(UNC_Moderate_expParameters_timeSeries$groupSize[i]))
  the = as.numeric(as.character(UNC_Moderate_expParameters_timeSeries$theta[i]))
  la = as.numeric(as.character(UNC_Moderate_expParameters_timeSeries$lambda[i]))
  location = which(UNC_Moderate_expParameters_timeSeries2$groupSize==g&UNC_Moderate_expParameters_timeSeries2$theta==the&UNC_Moderate_expParameters_timeSeries2$lambda==la)
  if(length(location)>0) {
  UNC_Moderate_expParameters_timeSeries2$optimChoiceRate[location] = UNC_Moderate_expParameters_timeSeries[i,4:73]
  }
}
UNC_Moderate_expParameters_timeSeries2$optimChoiceRate = as.numeric(as.character(UNC_Moderate_expParameters_timeSeries2$optimChoiceRate))

for(t in 1:70) {
  UNC_Moderate_expParameters_timeSeries2$optimChoiceRate[which(UNC_Moderate_expParameters_timeSeries2$groupSize==1&UNC_Moderate_expParameters_timeSeries2$round==t)] = mean(UNC_Moderate_expParameters_timeSeries2$optimChoiceRate[which(UNC_Moderate_expParameters_timeSeries2$groupSize==1&UNC_Moderate_expParameters_timeSeries2$round==t)])
}

Figure1_a =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, lambda==low_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize_factor, colour=groupSize_factor, linetype=groupSize_factor),size=1.5)+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.4\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_manual(name = 'Group size', values = c("black", "orange", "red"))+
  scale_linetype_manual(name = 'Group size', values = c("dotted", "dashed", "solid"))+
  #scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()+
  theme(legend.position = c(0.06, 0.2), legend.key.width = unit(15, "mm"))
Figure1_b =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==3 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==10 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==30 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.4\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()
Figure1_e =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==3 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==10 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==30 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()
Figure1_f =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==3 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==10 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_expParameters_timeSeries2, groupSize==30 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()

Figure1 = plot_grid(Figure1_a, Figure1_b, Figure1_e, Figure1_f, labels = c('','','',''), ncol = 2, align = 'v')
ggsave(file = "Figure1.pdf", plot = Figure1, dpi = 600, width = 9, height = 6)
ggsave(file = "Figure1.tiff", plot = Figure1, dpi = 600, width = 9, height = 6)

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
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
  myTheme_legend()

ggsave(file = "Figure2.pdf", plot = Figure2, dpi = 600, width = 9, height = 6)
ggsave(file = "Figure2.tiff", plot = Figure2, dpi = 600, width = 9, height = 6)



#######################################################################
##
## Figure 4
## Experimental Results
##
#######################################################################
## taskDifficulty3
allBehaviouralData2$taskDifficulty3 = allBehaviouralData2$taskDifficulty
allBehaviouralData2$taskDifficulty3[which(allBehaviouralData2$taskDifficulty=='25%')] <- 'Low Uncertainty'
allBehaviouralData2$taskDifficulty3[which(allBehaviouralData2$taskDifficulty=='50%')] <- 'Moderate Uncertainty'
allBehaviouralData2$taskDifficulty3[which(allBehaviouralData2$taskDifficulty=='78%')] <- 'High Uncertainty'
allBehaviouralData2$taskDifficulty3 = factor(allBehaviouralData2$taskDifficulty3, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))
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

## Post-hoc simulation
## import simulation_fulModel_all.csv
simulation_fullModel_all = rbind(rbind(simulation_fullModel_25, simulation_fullModel_50), simulation_fullModel_78)
simulation_fullModel_all$taskDifficulty = c(rep('Low Uncertainty', nrow(simulation_fullModel_25)),rep('Moderate Uncertainty', nrow(simulation_fullModel_50)),rep('High Uncertainty', nrow(simulation_fullModel_78)))
simulation_fullModel_all$taskDifficulty = factor(simulation_fullModel_all$taskDifficulty, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))

simulation_fullModel_all$sizeCategoryBinary = NA
simulation_fullModel_all$sizeCategoryBinary[which(simulation_fullModel_all$taskDifficulty=='Low Uncertainty'&simulation_fullModel_all$groupSize<=median_low)] = 'small'
simulation_fullModel_all$sizeCategoryBinary[which(simulation_fullModel_all$taskDifficulty=='Low Uncertainty'&simulation_fullModel_all$groupSize >median_low)] = 'large'
simulation_fullModel_all$sizeCategoryBinary[which(simulation_fullModel_all$taskDifficulty=='Moderate Uncertainty'&simulation_fullModel_all$groupSize<=median_moderate)] = 'small'
simulation_fullModel_all$sizeCategoryBinary[which(simulation_fullModel_all$taskDifficulty=='Moderate Uncertainty'&simulation_fullModel_all$groupSize >median_moderate)] = 'large'
simulation_fullModel_all$sizeCategoryBinary[which(simulation_fullModel_all$taskDifficulty=='High Uncertainty'&simulation_fullModel_all$groupSize<=median_high)] = 'small'
simulation_fullModel_all$sizeCategoryBinary[which(simulation_fullModel_all$taskDifficulty=='High Uncertainty'&simulation_fullModel_all$groupSize >median_high)] = 'large'
simulation_fullModel_all$sizeCategoryBinary[which(simulation_fullModel_all$groupSize==1)] = 'single'
simulation_fullModel_all$sizeCategoryBinary = factor(simulation_fullModel_all$sizeCategoryBinary, levels=c('single', 'small', 'large'))
## taskDifficulty3
simulation_fullModel_all$taskDifficulty3 = simulation_fullModel_all$taskDifficulty
#simulation_fullModel_all$taskDifficulty3[which(simulation_fullModel_all$taskDifficulty=='Easy')] <- 'Low'
#simulation_fullModel_all$taskDifficulty3[which(simulation_fullModel_all$taskDifficulty=='Moderate')] <- 'Moderate'
#simulation_fullModel_all$taskDifficulty3[which(simulation_fullModel_all$taskDifficulty=='Difficult')] <- 'High'
#simulation_fullModel_all$taskDifficulty3 = factor(simulation_fullModel_all$taskDifficulty3, levels=c('Low','Moderate','High'))
## Line plot of simulation





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

performanceSummary4$taskDifficulty3 = NA
performanceSummary4$taskDifficulty3[which(performanceSummary4$taskDifficulty=='25%')] <- 'Low Uncertainty'
performanceSummary4$taskDifficulty3[which(performanceSummary4$taskDifficulty=='50%')] <- 'Moderate Uncertainty'
performanceSummary4$taskDifficulty3[which(performanceSummary4$taskDifficulty=='78%')] <- 'High Uncertainty'
performanceSummary4$taskDifficulty3 = factor(performanceSummary4$taskDifficulty3, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))

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
                    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

Figure4_bottomrow = ggplot() +
  #stat_summary(data=allBehaviouralData2, mapping=aes(round, optimalChoice, group=sizeCategoryBinary, colour=sizeCategoryBinary),fun.y = mean,geom="line",size=0.6,alpha=0.4)+
  #stat_summary(data=allBehaviouralData2, mapping=aes(round, optimalChoice, group=sizeCategoryBinary, colour=sizeCategoryBinary, shape=sizeCategoryBinary),fun.y = mean,geom="point",size=0.9,alpha=0.7)+
  stat_summary(data=simulation_fullModel_all, mapping=aes(t, isThisBestOption, group=sizeCategoryBinary, colour=sizeCategoryBinary, linetype=sizeCategoryBinary),fun.y = mean,geom="line",size=1.2)+
  #geom_point(data=allBehaviouralData2, mapping=aes(round, optimalChoice, group=sizeCategoryBinary, colour=sizeCategoryBinary, shape=sizeCategoryBinary), size=2)+
  myTheme_legend()+
  ylim(c(0,1))+
    labs(x = "Rounds", y = "Probability of choosing\n the best option")+
    facet_grid(.~taskDifficulty3)+
    scale_color_manual(name = 'Group size', values = c("single" = "grey30", "small" = "orange", "large" = "red"))+
    scale_linetype_manual(name = 'Group size', values=c("single" = "dashed", "small" = "twodash", "large" = "solid"))+
    scale_shape_manual(name = 'Group size', values=c("single" = 15, "small" = 16, "large" = 17))+
    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)


Figure_4_pre = plot_grid(Figure4_toprow, Figure4_bottomrow, labels=c("",""), ncol = 1, align = 'v')

fig4bottom_lowUncertain_exp = ggplot() +
    stat_summary(data=subset(allBehaviouralData2, taskDifficulty3=='Low Uncertainty'), mapping=aes(round, optimalChoice, group=sizeCategoryBinary, colour=sizeCategoryBinary),fun.y = mean,geom="line",size=0.6,alpha=0.4)+
    stat_summary(data=subset(allBehaviouralData2, taskDifficulty3=='Low Uncertainty'), mapping=aes(round, optimalChoice, group=sizeCategoryBinary, colour=sizeCategoryBinary, shape=sizeCategoryBinary),fun.y = mean,geom="point",size=0.9,alpha=0.7)+
    myTheme_small()+
    ylim(c(0,1))+
      labs(x = "", y = "")+
      scale_color_manual(name = 'Group size', values = c("single" = "grey30", "small" = "orange", "large" = "red"))+
      scale_linetype_manual(name = 'Group size', values=c("single" = "dashed", "small" = "twodash", "large" = "solid"))+
      scale_shape_manual(name = 'Group size', values=c("single" = 15, "small" = 16, "large" = 17))+
      panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
  fig4bottom_moderateUncertain_exp = ggplot() +
    stat_summary(data=subset(allBehaviouralData2, taskDifficulty3=='Moderate Uncertainty'), mapping=aes(round, optimalChoice, group=sizeCategoryBinary, colour=sizeCategoryBinary),fun.y = mean,geom="line",size=0.6,alpha=0.4)+
    stat_summary(data=subset(allBehaviouralData2, taskDifficulty3=='Moderate Uncertainty'), mapping=aes(round, optimalChoice, group=sizeCategoryBinary, colour=sizeCategoryBinary, shape=sizeCategoryBinary),fun.y = mean,geom="point",size=0.9,alpha=0.7)+
    myTheme_small()+
    ylim(c(0,1))+
      labs(x = "", y = "")+
      scale_color_manual(name = 'Group size', values = c("single" = "grey30", "small" = "orange", "large" = "red"))+
      scale_linetype_manual(name = 'Group size', values=c("single" = "dashed", "small" = "twodash", "large" = "solid"))+
      scale_shape_manual(name = 'Group size', values=c("single" = 15, "small" = 16, "large" = 17))+
      panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
  fig4bottom_highUncertain_exp = ggplot() +
    stat_summary(data=subset(allBehaviouralData2, taskDifficulty3=='High Uncertainty'), mapping=aes(round, optimalChoice, group=sizeCategoryBinary, colour=sizeCategoryBinary),fun.y = mean,geom="line",size=0.6,alpha=0.4)+
    stat_summary(data=subset(allBehaviouralData2, taskDifficulty3=='High Uncertainty'), mapping=aes(round, optimalChoice, group=sizeCategoryBinary, colour=sizeCategoryBinary, shape=sizeCategoryBinary),fun.y = mean,geom="point",size=0.9,alpha=0.7)+
    myTheme_small()+
    ylim(c(0,1))+
      labs(x = "", y = "")+
      scale_color_manual(name = 'Group size', values = c("single" = "grey30", "small" = "orange", "large" = "red"))+
      scale_linetype_manual(name = 'Group size', values=c("single" = "dashed", "small" = "twodash", "large" = "solid"))+
      scale_shape_manual(name = 'Group size', values=c("single" = 15, "small" = 16, "large" = 17))+
      panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

Figure4_merged = ggdraw() +
  #draw_plot(fig4bottom_lowUncertain_sim + theme(legend.justification = "bottom"), 0, 0, 1, 1) +
  #draw_plot(Figure4_toprow + theme(legend.justification = "bottom"), 0, 0, 1, 1) +
  draw_plot(Figure_4_pre, 0, 0, 1, 1) +
  draw_plot(fig4bottom_lowUncertain_exp + 
              theme(legend.justification = "top"), 0.089, 0.0425, 0.155, 0.2)+#0.065, 0.085, 0.16, 0.4
  draw_plot(fig4bottom_moderateUncertain_exp + 
              theme(legend.justification = "top"), 0.345, 0.0425, 0.155, 0.2)+#0.33, 0.085, 0.16, 0.4
  draw_plot(fig4bottom_highUncertain_exp + 
              theme(legend.justification = "top"), 0.59, 0.219, 0.155, 0.2)#0.59, 0.085, 0.16, 0.4

ggsave(file = "Figure4_merged.pdf", plot = Figure4_merged, dpi = 600, width = 9*1.2, height = 6*1.2)
ggsave(file = "Figure4_merged.tiff", plot = Figure4_merged, dpi = 600, width = 9*1.2, height = 6*1.2)

#ggsave(file = "Figure4.pdf", plot = Figure4, dpi = 600, width = 9*1.2, height = 6*1.2)
#ggsave(file = "Figure4.tiff", plot = Figure4, dpi = 600, width = 9*1.2, height = 6*1.2)


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
## Figure 3
## Model Fitting Results
##
#######################################################################

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
  chains=4, iter=5000, warmup=1500, thin=2
  #chains=4, iter=500, warmup=150, thin=2
)
isPositiveCopier_fitting
saveRDS(isPositiveCopier_fitting, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/isPositiveCopier_fitting.rds')
isPositiveCopier_mcmc_result <- rstan::extract(isPositiveCopier_fitting)
plot(isPositiveCopier_fitting, pars=c('beta','sigma_r','sigma_r'))
plot(isPositiveCopier_fitting, pars=c('isPositiveCopier_pred_25'))
  ## 事後診断
traceplot(isPositiveCopier_fitting, pars=c('beta','sigma_e'), inc_warmup=FALSE)

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
  myTheme_inner()+ylim(c(NA,1))+xlim(c(2, 28))+
  theme(legend.position = c(0.05, 0.5))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
GS_proportionSL_50 = ggplot() +
  geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=1), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_50, aes(groupSize, proportion_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_50,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_50,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_50, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size', y='Proportion of \nstrategies', title='Moderate Uncertainty') +
  myTheme()+ylim(c(NA,1))+xlim(c(2, 28))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
GS_proportionSL_78 = ggplot() +
  geom_ribbon(data=subset(isPositiveCopier_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=2.2), aes(x=groupSize, ymin=isPositiveCopier_p2.5,ymax=isPositiveCopier_p97.5), fill = "lightpink")+ #fill = "lightpink"
  geom_path(data=subset(isPositiveCopier_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.4&groupSize_standerized2<=2.2), aes(x=groupSize, y=isPositiveCopier_p50), linetype='dashed', size=0.6,color='red')+
  geom_line(data=strategies_proportion_78, aes(groupSize, proportion_pos_freq_dep),colour='red') + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_pos_freq_dep),colour='red',shape=2) +
  geom_line(data=strategies_proportion_78,aes(groupSize, proportion_neg_freq_dep),colour='blue') + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_neg_freq_dep),colour='blue',shape=18) +
  geom_line(data=strategies_proportion_78,aes(groupSize, proportion_random_copy),colour='grey50', linetype="dashed") + geom_point(data=strategies_proportion_78, aes(groupSize, proportion_random_copy),colour='grey50') +
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Strength of \nfrequency dependence") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Strength of \nfrequency dependence") +
  labs(x='Group size', y='Proportion of \nstrategies', title='High Uncertainty') +
  myTheme()+ylim(c(NA,1))+xlim(c(2, 28))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

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
  labs(x='Group size', y='Proportion of \nstrategies', title='Moderate Uncertainty') +
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
  labs(x='Group size', y='Proportion of \nstrategies', title='High Uncertainty') +
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
  chains=4, iter=5000, warmup=1500, thin=2
  #chains=4, iter=500, warmup=150, thin=2
)
average_copy_rate_fitting
saveRDS(average_copy_rate_fitting, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/average_copy_rate_fitting.rds')
average_copy_rate_mcmc_result <- rstan::extract(average_copy_rate_fitting)
plot(average_copy_rate_fitting, pars=c('beta','gamma','sigma_e','sigma_r'))
plot(average_copy_rate_fitting, pars=c('average_copy_rate_pred_25'))

  ## 事後診断
traceplot(average_copy_rate_fitting, pars=c('beta','sigma_e'), inc_warmup=FALSE)

average_copy_rate_pred_data = cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_copy_rate_mcmc_result$average_copy_rate_pred_25[,1:ncol(average_copy_rate_mcmc_result$average_copy_rate_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
average_copy_rate_pred_data = rbind(average_copy_rate_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_copy_rate_mcmc_result$average_copy_rate_pred_50[,1:ncol(average_copy_rate_mcmc_result$average_copy_rate_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
average_copy_rate_pred_data = rbind(average_copy_rate_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(average_copy_rate_mcmc_result$average_copy_rate_pred_78[,1:ncol(average_copy_rate_mcmc_result$average_copy_rate_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
average_copy_rate_pred_data$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(average_copy_rate_pred_data) <- c("groupSize_standerized2",
              "average_copy_rate_p2.5", "average_copy_rate_p25", "average_copy_rate_p50", "average_copy_rate_p75", "average_copy_rate_p97.5","taskDifficulty")

Figure3_b_25 = ggplot() +
  geom_ribbon(data=subset(average_copy_rate_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=average_copy_rate_p25,ymax=average_copy_rate_p75), fill = "grey80")+ #fill = "lightpink"
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&frequency_dependence_4!='Single-condition'&frequency_dependence_4!='--'&taskDifficulty=='25%'),mapping=aes(groupSize_standerized2, average_copy_rate,colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_path(data=subset(average_copy_rate_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=average_copy_rate_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
  labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='') +
  myTheme_inner()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  theme(legend.position = c(0.02, 0.8))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
Figure3_b_50 = ggplot() +
  geom_ribbon(data=subset(average_copy_rate_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=average_copy_rate_p25,ymax=average_copy_rate_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='50%'),mapping=aes(groupSize_standerized2, average_copy_rate,colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_path(data=subset(average_copy_rate_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=average_copy_rate_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
Figure3_b_78 = ggplot() +
  geom_ribbon(data=subset(average_copy_rate_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, ymin=average_copy_rate_p25,ymax=average_copy_rate_p75), fill = "grey80")+
  geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='78%'),mapping=aes(groupSize_standerized2, average_copy_rate,colour=frequency_dependence_4, shape=frequency_dependence_4))+
  geom_path(data=subset(average_copy_rate_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, y=average_copy_rate_p50), linetype='dashed', size=0.6,color='grey20')+
  scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies") +
  scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies") +
  labs(x='Group size\n (standardized)', y=bquote(atop('Mean','social learning weight' ~ bar(sigma[i]))), title='') +
  myTheme()+ylim(c(0,1))+xlim(c(-1.3,2.3))+
  panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

Figure3_b = plot_grid(Figure3_b_25, Figure3_b_50, Figure3_b_78, labels = c('','',''), ncol = 3, align = 'v')
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
  #chains=4, iter=500, warmup=150, thin=2
)
average_copy_rate_fitting_posOnly
saveRDS(average_copy_rate_fitting_posOnly, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/average_copy_rate_fitting_posOnly.rds')
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
  chains=4, iter=5000, warmup=1500, thin=2
)
conformity_exponent_fitting
saveRDS(conformity_exponent_fitting, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/conformity_exponent_fitting.rds')
conformity_exponent_mcmc_result <- rstan::extract(conformity_exponent_fitting)
plot(conformity_exponent_fitting, pars=c('beta','gamma','sigma_e','sigma_r'))
plot(conformity_exponent_fitting, pars=c('conformity_exponent_pred_25'))

  ## 事後診断
traceplot(conformity_exponent_fitting, pars=c('beta','sigma_e'), inc_warmup=TRUE)

conformity_exponent_pred_data = cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(conformity_exponent_mcmc_result$conformity_exponent_pred_25[,1:ncol(conformity_exponent_mcmc_result$conformity_exponent_pred_25)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
conformity_exponent_pred_data = rbind(conformity_exponent_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(conformity_exponent_mcmc_result$conformity_exponent_pred_50[,1:ncol(conformity_exponent_mcmc_result$conformity_exponent_pred_50)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
conformity_exponent_pred_data = rbind(conformity_exponent_pred_data, cbind(parameters_stan_data$test_group_size, as.data.frame(t(apply(conformity_exponent_mcmc_result$conformity_exponent_pred_78[,1:ncol(conformity_exponent_mcmc_result$conformity_exponent_pred_78)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))) )
conformity_exponent_pred_data$taskDifficulty = c(rep('25%',200),rep('50%',200),rep('78%',200))
colnames(conformity_exponent_pred_data) <- c("groupSize_standerized2",
              "conformity_exponent_p2.5", "conformity_exponent_p25", "conformity_exponent_p50", "conformity_exponent_p75", "conformity_exponent_p97.5","taskDifficulty")


## conformity exponent
Figure3_c_25 = ggplot() +
	geom_ribbon(data=subset(conformity_exponent_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=conformity_exponent_p25,ymax=conformity_exponent_p75), fill = "grey80")+#fill = "lightpink"
	geom_hline(yintercept=1, linetype='solid', size=0.7, alpha=0.7)+
	geom_hline(yintercept=-1, linetype='solid', size=0.7, alpha=0.7)+
	geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&frequency_dependence_4!='Single-condition'&frequency_dependence_4!='--'&taskDifficulty=='25%'),mapping=aes(groupSize_standerized2, theta, colour=frequency_dependence_4, shape=frequency_dependence_4))+
	geom_line(data=subset(conformity_exponent_pred_data,taskDifficulty=='25%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=conformity_exponent_p50), linetype='dashed', size=0.6,color='grey20')+
	scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
	scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="", labels=c("Random-copying"="Random choice", "Neg-freq-dep"="Negative freq. dep.","Pos-freq-dep"="Positive freq. dep.")) +
	labs(x='Group size\n (standardized)', y=bquote(atop('Conformity', 'exponent' ~ theta[i])), title='') +
	myTheme_inner()+ylim(c(-7,12))+xlim(c(-1.3,2.3))+
	theme(legend.position = c(0.02, 0.8))+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
Figure3_c_50 = ggplot() +
	geom_ribbon(data=subset(conformity_exponent_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, ymin=conformity_exponent_p25,ymax=conformity_exponent_p75), fill = "grey80")+
	geom_hline(yintercept=1, linetype='solid', size=0.7, alpha=0.7)+
	geom_hline(yintercept=-1, linetype='solid', size=0.7, alpha=0.7)+
	geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='50%'),mapping=aes(groupSize_standerized2, theta, colour=frequency_dependence_4, shape=frequency_dependence_4))+
	geom_line(data=subset(conformity_exponent_pred_data,taskDifficulty=='50%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=1), aes(x=groupSize_standerized2, y=conformity_exponent_p50), linetype='dashed', size=0.6,color='grey20')+
	scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies") +
	scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies") +
	labs(x='Group size\n (standardized)', y=bquote(atop('Conformity', 'exponent' ~ theta[i])), title='') +
	myTheme()+ylim(c(-7,12))+xlim(c(-1.3,2.3))+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)
Figure3_c_78 = ggplot() +
	geom_ribbon(data=subset(conformity_exponent_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, ymin=conformity_exponent_p25,ymax=conformity_exponent_p75), fill = "grey80")+
	geom_hline(yintercept=1, linetype='solid', size=0.7, alpha=0.7)+
	geom_hline(yintercept=-1, linetype='solid', size=0.7, alpha=0.7)+
	geom_point(data=subset(performanceSummary4, environment=='e1+e2'&frequency_dependence_4!='xPos-freq-dep'&taskDifficulty=='78%'),mapping=aes(groupSize_standerized2, theta, colour=frequency_dependence_4, shape=frequency_dependence_4))+
	geom_line(data=subset(conformity_exponent_pred_data,taskDifficulty=='78%'&groupSize_standerized2>=-1.3&groupSize_standerized2<=2.2), aes(x=groupSize_standerized2, y=conformity_exponent_p50), linetype='dashed', size=0.6,color='grey20')+
	scale_colour_manual(values=c("Single-condition"="black","Asocial"="grey50","Random-copying"="black",'Neg-freq-dep'="blue", "Pos-freq-dep"="red", "Random choice"="grey80"), name="Categories of\nlearning strategies") +
	scale_shape_manual(values=c("Single-condition"=1,"Asocial"=1,"Random-copying"=1,'Neg-freq-dep'=18, "Pos-freq-dep"=2, "Random choice"=1), name="Categories of\nlearning strategies") +
	labs(x='Group size\n (standardized)', y=bquote(atop('Conformity', 'exponent' ~ theta[i])), title='') +
	myTheme()+ylim(c(-7,12))+xlim(c(-1.3,2.3))+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

Figure3_c = plot_grid(Figure3_c_25, Figure3_c_50, Figure3_c_78, labels = c('','',''), ncol = 3, align = 'v')
#Figure3_c

## summary
#Figure3 = plot_grid(Figure3_top, Figure3_b, Figure3_c, labels=c("","",""), ncol = 1, align = 'v')
Figure3 = plot_grid(GS_proportionSL_25,GS_proportionSL_50,GS_proportionSL_78,Figure3_b_25, Figure3_b_50, Figure3_b_78,Figure3_c_25, Figure3_c_50, Figure3_c_78, labels=c("","",""), ncol = 3, align = 'v')
#Figure3
ggsave(file = "Figure3.pdf", plot = Figure3, dpi = 700, width = 10, height = 10)
ggsave(file = "Figure3.tiff", plot = Figure3, dpi = 700, width = 10, height = 10)



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
saveRDS(conformity_exponent_fitting_posOnly, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/conformity_exponent_fitting_posOnly.rds')
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
## Other parameters????
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
saveRDS(alpha_fitting, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/alpha_fitting.rds')
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
saveRDS(average_invTemp_fitting, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/average_invTemp_fitting.rds')
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
low_sigma = 0.4
high_sigma = 0.9
low_theta = 1
high_theta = 3

## To get the tidy-data frame for the time series data
UNC_Moderate_otherParam1_timeSeries_sup1 = subset(UNC_Moderate_otherParam1_timeSeries_sup10, theta==low_theta|theta==high_theta)
UNC_Moderate_otherParam1_timeSeries_sup12 = data.frame(
  round = rep(1:70, 2*2*3),
  groupSize = c(rep(3, 70*2*2), rep(10, 70*2*2), rep(30, 70*2*2)),
  groupSize_factor = c(rep('n = 3', 70*2*2), rep('n = 10', 70*2*2), rep('n = 30', 70*2*2)),
  theta = rep(c(rep(low_theta, 70*2), rep(high_theta,70*2)), 3),
  #theta = rep(c(rep(1, 70*2), rep(3,70*2)), 3),
  lambda = rep(rep(c(rep(low_sigma,70), rep(high_sigma,70)), 2), 3)
)
UNC_Moderate_otherParam1_timeSeries_sup12$groupSize_factor = factor(UNC_Moderate_otherParam1_timeSeries_sup12$groupSize_factor, levels=c('n = 3','n = 10','n = 30'))
UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate = NA
for(i in 1:nrow(UNC_Moderate_otherParam1_timeSeries_sup1)) {
  g = as.numeric(as.character(UNC_Moderate_otherParam1_timeSeries_sup1$groupSize[i]))
  the = as.numeric(as.character(UNC_Moderate_otherParam1_timeSeries_sup1$theta[i]))
  la = as.numeric(as.character(UNC_Moderate_otherParam1_timeSeries_sup1$lambda[i]))
  location = which(UNC_Moderate_otherParam1_timeSeries_sup12$groupSize==g&UNC_Moderate_otherParam1_timeSeries_sup12$theta==the&UNC_Moderate_otherParam1_timeSeries_sup12$lambda==la)
  if(length(location)>0) {
  UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate[location] = UNC_Moderate_otherParam1_timeSeries_sup1[i,4:73]
  }
}
UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate = as.numeric(as.character(UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate))

for(t in 1:70) {
  UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate[which(UNC_Moderate_otherParam1_timeSeries_sup12$groupSize==1&UNC_Moderate_otherParam1_timeSeries_sup12$round==t)] = mean(UNC_Moderate_otherParam1_timeSeries_sup12$optimChoiceRate[which(UNC_Moderate_otherParam1_timeSeries_sup12$groupSize==1&UNC_Moderate_otherParam1_timeSeries_sup12$round==t)])
}

Figure1_otherParam1_a =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, lambda==low_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize_factor, colour=groupSize_factor, linetype=groupSize_factor),size=1.5)+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.4\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_manual(name = 'Group size', values = c("black", "orange", "red"))+
  scale_linetype_manual(name = 'Group size', values = c("dotted", "dashed", "solid"))+
  #scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()+
  theme(legend.position = c(0.06, 0.2), legend.key.width = unit(15, "mm"))
Figure1_otherParam1_b =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==3 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==10 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==30 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.4\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()
Figure1_otherParam1_e =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==3 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==10 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==30 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()
Figure1_otherParam1_f =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==3 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==10 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam1_timeSeries_sup12, groupSize==30 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()

Figure1_otherParam1 = plot_grid(Figure1_otherParam1_a, Figure1_otherParam1_b, Figure1_otherParam1_e, Figure1_otherParam1_f, labels = c('','','',''), ncol = 2, align = 'v')
ggsave(file = "Figure1_otherParam1.pdf", plot = Figure1_otherParam1, dpi = 600, width = 9, height = 6)
#ggsave(file = "Figure1_otherParam1.tiff", plot = Figure1_otherParam1, dpi = 600, width = 9, height = 6)

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
low_sigma = 0.4
high_sigma = 0.9
low_theta = 1
high_theta = 3

## To get the tidy-data frame for the time series data
UNC_Moderate_otherParam2_timeSeries_sup1 = subset(UNC_Moderate_otherParam2_timeSeries_sup10, theta==low_theta|theta==high_theta)
UNC_Moderate_otherParam2_timeSeries_sup12 = data.frame(
  round = rep(1:70, 2*2*3),
  groupSize = c(rep(3, 70*2*2), rep(10, 70*2*2), rep(30, 70*2*2)),
  groupSize_factor = c(rep('n = 3', 70*2*2), rep('n = 10', 70*2*2), rep('n = 30', 70*2*2)),
  theta = rep(c(rep(low_theta, 70*2), rep(high_theta,70*2)), 3),
  #theta = rep(c(rep(1, 70*2), rep(3,70*2)), 3),
  lambda = rep(rep(c(rep(low_sigma,70), rep(high_sigma,70)), 2), 3)
)
UNC_Moderate_otherParam2_timeSeries_sup12$groupSize_factor = factor(UNC_Moderate_otherParam2_timeSeries_sup12$groupSize_factor, levels=c('n = 3','n = 10','n = 30'))
UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate = NA
for(i in 1:nrow(UNC_Moderate_otherParam2_timeSeries_sup1)) {
  g = as.numeric(as.character(UNC_Moderate_otherParam2_timeSeries_sup1$groupSize[i]))
  the = as.numeric(as.character(UNC_Moderate_otherParam2_timeSeries_sup1$theta[i]))
  la = as.numeric(as.character(UNC_Moderate_otherParam2_timeSeries_sup1$lambda[i]))
  location = which(UNC_Moderate_otherParam2_timeSeries_sup12$groupSize==g&UNC_Moderate_otherParam2_timeSeries_sup12$theta==the&UNC_Moderate_otherParam2_timeSeries_sup12$lambda==la)
  if(length(location)>0) {
  UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate[location] = UNC_Moderate_otherParam2_timeSeries_sup1[i,4:73]
  }
}
UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate = as.numeric(as.character(UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate))

for(t in 1:70) {
  UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate[which(UNC_Moderate_otherParam2_timeSeries_sup12$groupSize==1&UNC_Moderate_otherParam2_timeSeries_sup12$round==t)] = mean(UNC_Moderate_otherParam2_timeSeries_sup12$optimChoiceRate[which(UNC_Moderate_otherParam2_timeSeries_sup12$groupSize==1&UNC_Moderate_otherParam2_timeSeries_sup12$round==t)])
}

Figure1_otherParam2_a =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, lambda==low_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize_factor, colour=groupSize_factor, linetype=groupSize_factor),size=1.5)+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.4\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_manual(name = 'Group size', values = c("black", "orange", "red"))+
  scale_linetype_manual(name = 'Group size', values = c("dotted", "dashed", "solid"))+
  #scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()+
  theme(legend.position = c(0.06, 0.2), legend.key.width = unit(15, "mm"))
Figure1_otherParam2_b =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==3 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==10 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==30 & lambda==low_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.4\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()
Figure1_otherParam2_e =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==3 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==10 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==30 & lambda==high_sigma & theta==low_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 1\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()
Figure1_otherParam2_f =
  ggplot() +
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==3 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='black',size=1.5,linetype='dotted')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==10 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='orange',size=1.5,linetype='dashed')+
  geom_line(data=subset(UNC_Moderate_otherParam2_timeSeries_sup12, groupSize==30 & lambda==high_sigma & theta==high_theta), mapping = aes(round, optimChoiceRate,group=groupSize),colour='red',size=1.5,linetype='solid')+
  annotate("text", x = 20, y = 0.9, label = "paste(bar(italic(sigma)), \" = 0.9\", \", \", bar(italic(theta)), \" = 3\" )", parse = TRUE, size=6, family="Times")+
  labs(x='Rounds', y='Proportion of choosing\n the best option')+
  scale_color_distiller(name = 'Group size', palette = "Spectral", direction = -1)+
  ylim(c(0,1))+
  geom_vline(xintercept = 40.5, linetype="dashed")+
  myTheme_inner()

Figure1_otherParam2 = plot_grid(Figure1_otherParam2_a, Figure1_otherParam2_b, Figure1_otherParam2_e, Figure1_otherParam2_f, labels = c('','','',''), ncol = 2, align = 'v')
ggsave(file = "Figure1_otherParam2.pdf", plot = Figure1_otherParam2, dpi = 600, width = 9, height = 6)
#ggsave(file = "Figure1_otherParam2.tiff", plot = Figure1_otherParam2, dpi = 600, width = 9, height = 6)

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



