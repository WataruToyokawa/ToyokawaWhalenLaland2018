# RESCOLA-WAGNER + SOFTMAX-ANNEALING CHOICE MODEL (model_changingBST.stan)
# PARAMETERS:
# 	alpha = learning rate (lower=0, upper=1)
# 	beta = inv_temp (lower=0)
# ASSUMPTIONS:
# 	Both alpha and beta are generated from normal distributions (hierarchical approach)

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
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(2, "lines")
  )
}
myTheme_legend = function() {
  theme(
    legend.position = 'right',
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Times" ), #"Times"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=15, family="Times" ,colour='black'),
    axis.text.y = element_text(size=15, family="Times" ,colour='black'),
    axis.title.x=element_text(size=16, family="Times" ),
    axis.title.y=element_text(size=16, family="Times" ),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white")
  )
}

# Static condition
library(rstan)
library(ggmcmc)
library(readr)
library(cowplot)
library(readr)
	## parallel coputing
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
	## WAIC function
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

	# WBIC when all subjects are included
WBIC = function(fit) {
	log_lik = rstan::extract(fit)$log_lik
	wbic = - mean(rowSums(log_lik))
	return(wbic)
}
	# WBIC for each individual subject
WBIC_indv = function(fit) {
	log_lik = rstan::extract(fit)$log_lik
	wbic = - apply(log_lik, 2, mean)
	return(wbic)
}

## ==========  Data reading and cleaning

allBehaviouralData_25 <- read_csv("allBehaviouralData_25.csv")
allBehaviouralData_50 <- read_csv("allBehaviouralData_50.csv")
allBehaviouralData_78 <- read_csv("allBehaviouralData_78.csv")

allBehaviouralData_all = rbind(rbind(allBehaviouralData_25, allBehaviouralData_50), allBehaviouralData_78[,1:22])

allBehaviouralData_all$taskDifficulty_num = as.numeric(as.factor(allBehaviouralData_all$taskDifficulty))
#testBehaviourData = allBehaviouralData_all
testBehaviourData = subset(allBehaviouralData_all, sizeCategory > 1) # N >= 2
# Using only the subjects who played the task at least for 31 rounds
testBehaviourData = subset(testBehaviourData, amazonID %in% names(which(table(testBehaviourData$amazonID)>30)))

# debug #
#testNames = head(names(table(subset(allBehaviouralData_all, taskDifficulty_num==1&sizeCategory>5)$amazonID)),10)
#testNames <- append(testNames, head(names(table(subset(allBehaviouralData_all, taskDifficulty_num==2&sizeCategory>5)$amazonID)), 10))
#testNames <- append(testNames, head(names(table(subset(allBehaviouralData_all, taskDifficulty_num==3&sizeCategory>5)$amazonID)), 10))
#testBehaviourData = subset(allBehaviouralData_all, amazonID %in% testNames)
# end - debug #



testBehaviourData$sub = as.numeric(as.factor(testBehaviourData$amazonID))

    # insert missed trials
for(i in 1:length(table(testBehaviourData$amazonID))) {
    thisSubject = subset(testBehaviourData,sub==i)
    for(t in 1:70) {
        if(nrow(thisSubject[which(thisSubject$round==t),])==0) {
            testBehaviourData <- rbind(testBehaviourData, rep(NA, ncol(testBehaviourData)))
            testBehaviourData$machine[nrow(testBehaviourData)] = -1
            testBehaviourData$result[nrow(testBehaviourData)] = -1
            testBehaviourData$round[nrow(testBehaviourData)] = t
            testBehaviourData$amazonID[nrow(testBehaviourData)] = testBehaviourData$amazonID[which(testBehaviourData$sub==i)][1]
            testBehaviourData$sub[nrow(testBehaviourData)] = testBehaviourData$sub[which(testBehaviourData$sub==i)][1]
            testBehaviourData$taskDifficulty_num[nrow(testBehaviourData)] = testBehaviourData$taskDifficulty_num[which(testBehaviourData$sub==i)][1]
        }
    }
}
testBehaviourData = testBehaviourData[order(testBehaviourData$round),] # Sorting by round

ReinforcementLearningStanData_group = list(
    All = nrow(testBehaviourData),
    Nsub = length(table(testBehaviourData$amazonID)),
    Ncue = 3, # number of options
    Ntrial = 70,
    Ncon = 3,
    sub = testBehaviourData$sub,
    Y = testBehaviourData$machine + 1,
    trial = testBehaviourData$round,
    outcome = testBehaviourData$result * 100
    )
con = numeric(length(table(testBehaviourData$amazonID)))
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    con[i] <- testBehaviourData$taskDifficulty_num[which(testBehaviourData$sub==i)][1]
}
ReinforcementLearningStanData_group$condition = con
startT = c()
MaxTrial = c()
for(i in 1:length(table(testBehaviourData$amazonID))){
    startT = append(startT, min(subset(testBehaviourData, sub==i)$round))
    MaxTrial = append(MaxTrial, length(which(subset(testBehaviourData, sub==i)$machine>=0)))
}
ReinforcementLearningStanData_group$startT = startT
ReinforcementLearningStanData_group$MaxTrial = MaxTrial

## compfun() で高速化
library(compiler)
F_calculation = function(array, Nsub, Ntrial, data)
{
    F <- array
    # F (Simulation data とは t+1 ずれてるから注意！)
    for(i in 1:Nsub) {
        F[i,1,1] = 0; F[i,2,1] = 0; F[i,3,1] = 0;
        for(t in 2:Ntrial) {
            lastChoice = 0
            if(subset(data,sub==i&round==(t-1))$machine+1>0) {
                lastChoice = subset(data,sub==i&round==(t-1))$machine + 1
            }
            F[i,1,t] = subset(data,sub==i&round==t)$socialFreq0
            F[i,2,t] = subset(data,sub==i&round==t)$socialFreq1
            F[i,3,t] = subset(data,sub==i&round==t)$socialFreq2
            if(lastChoice>0){
                F[i,lastChoice,t] = F[i,lastChoice,t] - 1
            }
            if(length(which(F[i,,t]>=0))==0){
                F[i,,t] <- c(-1,-1,-1)
            }
        }
    }
    return(F)
}

F_calculation.compiled <- cmpfun(F_calculation)
F0 = array(rep(NA,nrow(testBehaviourData)), c(ReinforcementLearningStanData_group$Nsub, ReinforcementLearningStanData_group$Ncue, ReinforcementLearningStanData_group$Ntrial))
F = F_calculation.compiled(F0,ReinforcementLearningStanData_group$Nsub,ReinforcementLearningStanData_group$Ntrial,testBehaviourData)

ReinforcementLearningStanData_group$F = F

## ==========  Data reading and cleaning END


################################################################################
##
##                   FITTING: FULL MODEL
##                  beta, sigma, and theta are all time dependent variavles (alpha is time independent)
##
################################################################################

    ## FITTING
model_changingBST = stan_model(file="~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model_changingBST.stan") # debug
fit_changingBST = sampling(
    model_changingBST, data=ReinforcementLearningStanData_group, seed=77,
    #pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_soc','s_soc','mu_theta','s_theta','mu_soc_reduction','s_soc_reduction','mu_annealing','s_annealing','alpha','beta','annealing','theta','soc_raw','soc_reduction','soc','netBeta','log_lik'),
    pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_annealing','s_annealing','mu_soc','s_soc','mu_soc_reduction','s_soc_reduction','mu_theta','s_theta','mu_theta_slope','s_theta_slope','alpha','beta','annealing','theta','theta_slope','soc_raw','soc_reduction','soc','netTheta','netBeta','log_lik'),
    init=function() {
        list(
            mu_beta = runif(3,0,4),# because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            mu_annealing = runif(3,0,4),# because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_alpha = runif(3,1,3),
            s_beta = runif(3,1,3), # because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_annealing = runif(3,1,3), # because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_alpha = runif(3,1,3),
            s_theta = runif(3,1,3),
            s_theta_slope = runif(3,1,3),
            s_soc = runif(3,1,3),
            s_soc_reduction = runif(3,1,3)#rep(2,3)#runif(1, 1, 4)
            )
    },
    #chains=8, iter=1200, warmup=600, thin=5
    chains=8, iter=2000, warmup=1000, thin=5
    #chains=8, iter=120, warmup=60, thin=1
)
fit_changingBST
saveRDS(fit_changingBST, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/fit_changingBST.rds')


ms_changingBST <- rstan::extract(fit_changingBST)
mcmc_results_plot_changingBST = plot(fit_changingBST, pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_soc','s_soc','mu_theta','s_theta','mu_theta_slope','s_theta_slope','mu_soc_reduction','s_soc_reduction','mu_annealing','s_annealing'))
trace_1_changingBST = traceplot(fit_changingBST, pars=c('mu_alpha','mu_beta','mu_annealing'), inc_warmup=FALSE)
trace_2_changingBST = traceplot(fit_changingBST, pars=c('mu_theta','mu_theta_slope','mu_soc','mu_soc_reduction'), inc_warmup=FALSE)
trace_3_changingBST = traceplot(fit_changingBST, pars=c('s_alpha','s_beta','s_theta','s_theta_slope','s_soc','s_soc_reduction'), inc_warmup=FALSE)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/mcmc_results_plot_changingBST.pdf", plot = mcmc_results_plot_changingBST, dpi = 600, width = 9, height = 9)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/trace_1_changingBST.pdf", plot = trace_1_changingBST, dpi = 600, width = 9, height = 9)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/trace_2_changingBST.pdf", plot = trace_2_changingBST, dpi = 600, width = 9, height = 9)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/trace_3_changingBST.pdf", plot = trace_3_changingBST, dpi = 600, width = 9, height = 9)


temp_changingBST = data.frame(WAICi_changingBST = WAIC_indv(fit_changingBST)$waic)

temp_changingBST <- cbind(temp_changingBST, as.data.frame(t(apply(ms_changingBST$alpha[,1:ncol(ms_changingBST$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST <- cbind(temp_changingBST, as.data.frame(t(apply(ms_changingBST$beta[,1:ncol(ms_changingBST$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST <- cbind(temp_changingBST, as.data.frame(t(apply(ms_changingBST$annealing[,1:ncol(ms_changingBST$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST <- cbind(temp_changingBST, as.data.frame(t(apply(ms_changingBST$theta[,1:ncol(ms_changingBST$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST <- cbind(temp_changingBST, as.data.frame(t(apply(ms_changingBST$theta_slope[,1:ncol(ms_changingBST$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST <- cbind(temp_changingBST, as.data.frame(t(apply(ms_changingBST$soc_reduction[,1:ncol(ms_changingBST$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST <- cbind(temp_changingBST, as.data.frame(t(apply(ms_changingBST$soc_raw[,1:ncol(ms_changingBST$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(temp_changingBST) <- c("changingBSTi",
                            "changingBST_alpha_p2.5", "changingBST_alpha_p25", "changingBST_alpha_p50", "changingBST_alpha_p75", "changingBST_alpha_p97.5",
                            "changingBST_beta_p2.5", "changingBST_beta_p25", "changingBST_beta_p50", "changingBST_beta_p75", "changingBST_beta_p97.5",
                            "changingBST_annealing_p2.5", "changingBST_annealing_p25", "changingBST_annealing_p50", "changingBST_annealing_p75", "changingBST_annealing_p97.5",
                            "changingBST_theta_p2.5", "changingBST_theta_p25", "changingBST_theta_p50", "changingBST_theta_p75", "changingBST_theta_p97.5",
                            "changingBST_theta_slope_p2.5", "changingBST_theta_slope_p25", "changingBST_theta_slope_p50", "changingBST_theta_slope_p75", "changingBST_theta_slope_p97.5",
                            "changingBST_soc_reduction_p2.5", "changingBST_soc_reduction_p25", "changingBST_soc_reduction_p50", "changingBST_soc_reduction_p75", "changingBST_soc_reduction_p97.5",
                            "changingBST_soc_raw_p2.5", "changingBST_soc_raw_p25", "changingBST_soc_raw_p50", "changingBST_soc_raw_p75", "changingBST_soc_raw_p97.5")
max_log_lik = apply(ms_changingBST$log_lik[,1:ncol(ms_changingBST$log_lik)], 2, which.max)
max_changingBST_alpha=rep(NA, ncol(ms_changingBST$log_lik))
max_changingBST_beta=rep(NA, ncol(ms_changingBST$log_lik))
max_changingBST_annealing=rep(NA, ncol(ms_changingBST$log_lik))
max_changingBST_theta=rep(NA, ncol(ms_changingBST$log_lik))
max_changingBST_theta_slope=rep(NA, ncol(ms_changingBST$log_lik))
max_changingBST_soc_reduction=rep(NA, ncol(ms_changingBST$log_lik))
max_changingBST_soc_raw=rep(NA, ncol(ms_changingBST$log_lik))
for(i in 1:ncol(ms_changingBST$log_lik)){
    max_changingBST_alpha[i] <- ms_changingBST$alpha[max_log_lik[i],i]
    max_changingBST_beta[i] <- ms_changingBST$beta[max_log_lik[i],i]
    max_changingBST_annealing[i] <- ms_changingBST$annealing[max_log_lik[i],i]
    max_changingBST_theta[i] <- ms_changingBST$theta[max_log_lik[i],i]
    max_changingBST_theta_slope[i] <- ms_changingBST$theta_slope[max_log_lik[i],i]
    max_changingBST_soc_reduction[i] <- ms_changingBST$soc_reduction[max_log_lik[i],i]
    max_changingBST_soc_raw[i] <- ms_changingBST$soc_raw[max_log_lik[i],i]
}
temp_changingBST$max_changingBST_alpha = max_changingBST_alpha
temp_changingBST$max_changingBST_beta = max_changingBST_beta
temp_changingBST$max_changingBST_annealing = max_changingBST_annealing
temp_changingBST$max_changingBST_theta = max_changingBST_theta
temp_changingBST$max_changingBST_theta_slope = max_changingBST_theta_slope
temp_changingBST$max_changingBST_soc_reduction = max_changingBST_soc_reduction
temp_changingBST$max_changingBST_soc_raw = max_changingBST_soc_raw
temp_changingBST$max_changingBST_log_lik = apply(ms_changingBST$log_lik[,1:ncol(ms_changingBST$log_lik)], 2, max)
numParam = 7
temp_changingBST$AIC_changingBST = -temp_changingBST$max_changingBST_log_lik + numParam ## 正確には AIC/2
temp_changingBST$BIC_changingBST = -temp_changingBST$max_changingBST_log_lik + (numParam/2)*log(70) ## more precisely, BIC/2
temp_changingBST$WAICi_changingBST = WAIC_indv(fit_changingBST)$waic

temp_changingBST$amazonID = NA
temp_changingBST$taskDifficulty = NA
temp_changingBST$totalRound = NA
temp_changingBST$correctNum = NA
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    temp_changingBST$amazonID[i] = testBehaviourData$amazonID[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)][1]
    temp_changingBST$taskDifficulty[i] = names(table(testBehaviourData$taskDifficulty[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)]))
    temp_changingBST$totalRound[i] = length(which(as.numeric(as.factor(subset(testBehaviourData, machine>=0)$amazonID))==i))
    temp_changingBST$correctNum[i] = sum(testBehaviourData$optimalChoice[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)], na.rm=TRUE)
}
temp_changingBST$correctProp = temp_changingBST$correctNum/temp_changingBST$totalRound

temp_changingBST$log_lik_n_eff = NA
temp_changingBST$log_lik_Rhat = NA
temp_changingBST$alpha_n_eff = NA
temp_changingBST$alpha_Rhat = NA
temp_changingBST$beta_n_eff = NA
temp_changingBST$beta_Rhat = NA
temp_changingBST$soc_raw_n_eff = NA
temp_changingBST$soc_raw_Rhat = NA
temp_changingBST$soc_reduction_n_eff = NA
temp_changingBST$soc_reduction_Rhat = NA
temp_changingBST$theta_n_eff = NA
temp_changingBST$theta_Rhat = NA
temp_changingBST$theta_slope_n_eff = NA
temp_changingBST$theta_slope_Rhat = NA
temp_changingBST$annealing_n_eff = NA
temp_changingBST$annealing_Rhat = NA
for(i in 1:nrow(temp_changingBST)){
    temp_changingBST$log_lik_n_eff[i] <- summary(fit_changingBST)$summary[paste('log_lik[',i,']',sep=''),'n_eff']
    temp_changingBST$log_lik_Rhat[i] <- summary(fit_changingBST)$summary[paste('log_lik[',i,']',sep=''),'Rhat']
    temp_changingBST$alpha_n_eff[i] <- summary(fit_changingBST)$summary[paste('alpha[',i,']',sep=''),'n_eff']
    temp_changingBST$alpha_Rhat[i] <- summary(fit_changingBST)$summary[paste('alpha[',i,']',sep=''),'Rhat']
    temp_changingBST$beta_n_eff[i] <- summary(fit_changingBST)$summary[paste('beta[',i,']',sep=''),'n_eff']
    temp_changingBST$beta_Rhat[i] <- summary(fit_changingBST)$summary[paste('beta[',i,']',sep=''),'Rhat']
    temp_changingBST$soc_raw_n_eff[i] <- summary(fit_changingBST)$summary[paste('soc_raw[',i,']',sep=''),'n_eff']
    temp_changingBST$soc_raw_Rhat[i] <- summary(fit_changingBST)$summary[paste('soc_raw[',i,']',sep=''),'Rhat']
    temp_changingBST$soc_reduction_n_eff[i] <- summary(fit_changingBST)$summary[paste('soc_reduction[',i,']',sep=''),'n_eff']
    temp_changingBST$soc_reduction_Rhat[i] <- summary(fit_changingBST)$summary[paste('soc_reduction[',i,']',sep=''),'Rhat']
    temp_changingBST$theta_n_eff[i] <- summary(fit_changingBST)$summary[paste('theta[',i,']',sep=''),'n_eff']
    temp_changingBST$theta_Rhat[i] <- summary(fit_changingBST)$summary[paste('theta[',i,']',sep=''),'Rhat']
    temp_changingBST$theta_slope_n_eff[i] <- summary(fit_changingBST)$summary[paste('theta_slope[',i,']',sep=''),'n_eff']
    temp_changingBST$theta_slope_Rhat[i] <- summary(fit_changingBST)$summary[paste('theta_slope[',i,']',sep=''),'Rhat']
    temp_changingBST$annealing_n_eff[i] <- summary(fit_changingBST)$summary[paste('annealing[',i,']',sep=''),'n_eff']
    temp_changingBST$annealing_Rhat[i] <- summary(fit_changingBST)$summary[paste('annealing[',i,']',sep=''),'Rhat']
}

write.csv(temp_changingBST,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/temp_changingBST.csv",
            row.names=FALSE)

    ## Estimated copying rate at each round
soc_table_changingBST = data.frame(subject = 1:ncol(ms_changingBST$log_lik), round = rep(1, ncol(ms_changingBST$log_lik)))
soc_table_changingBST = cbind(soc_table_changingBST, as.data.frame(t(apply(ms_changingBST$soc[,1:ncol(ms_changingBST$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(soc_table_changingBST) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
    soc_table_changingBST = rbind(soc_table_changingBST, cbind(cbind(1:ncol(ms_changingBST$log_lik), rep(t,ncol(ms_changingBST$log_lik))), as.data.frame(t(apply(ms_changingBST$soc[,1:ncol(ms_changingBST$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(soc_table_changingBST) <- c("subject", "round", "changingBST_soc_p2.5", "changingBST_soc_p25", "changingBST_soc_p50", "changingBST_soc_p75", "changingBST_soc_p97.5")
soc_table_changingBST$amazonID = NA
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    soc_table_changingBST$amazonID[which(soc_table_changingBST$subject==i)] = testBehaviourData$amazonID[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)][1]
}

write.csv(soc_table_changingBST,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/soc_table_changingBST.csv",
            row.names=FALSE)

    ## Estimated netBeta at each round
netBeta_table_changingBST = data.frame(subject = 1:ncol(ms_changingBST$log_lik), round = rep(1, ncol(ms_changingBST$log_lik)))
netBeta_table_changingBST = cbind(netBeta_table_changingBST, as.data.frame(t(apply(ms_changingBST$netBeta[,1:ncol(ms_changingBST$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(netBeta_table_changingBST) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
    netBeta_table_changingBST = rbind(netBeta_table_changingBST, cbind(cbind(1:ncol(ms_changingBST$log_lik), rep(t,ncol(ms_changingBST$log_lik))), as.data.frame(t(apply(ms_changingBST$netBeta[,1:ncol(ms_changingBST$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(netBeta_table_changingBST) <- c("subject", "round", "changingBST_netBeta_p2.5", "changingBST_netBeta_p25", "changingBST_netBeta_p50", "changingBST_netBeta_p75", "changingBST_netBeta_p97.5")
netBeta_table_changingBST$amazonID = NA
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    netBeta_table_changingBST$amazonID[which(netBeta_table_changingBST$subject==i)] = testBehaviourData$amazonID[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)][1]
}

write.csv(netBeta_table_changingBST,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/netBeta_table_changingBST.csv",
            row.names=FALSE)

    ## Estimated netTheta at each round
netTheta_table_changingBST = data.frame(subject = 1:ncol(ms_changingBST$log_lik), round = rep(1, ncol(ms_changingBST$log_lik)))
netTheta_table_changingBST = cbind(netTheta_table_changingBST, as.data.frame(t(apply(ms_changingBST$netTheta[,1:ncol(ms_changingBST$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(netTheta_table_changingBST) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
    netTheta_table_changingBST = rbind(netTheta_table_changingBST, cbind(cbind(1:ncol(ms_changingBST$log_lik), rep(t,ncol(ms_changingBST$log_lik))), as.data.frame(t(apply(ms_changingBST$netTheta[,1:ncol(ms_changingBST$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(netTheta_table_changingBST) <- c("subject", "round", "changingBST_netTheta_p2.5", "changingBST_netTheta_p25", "changingBST_netTheta_p50", "changingBST_netTheta_p75", "changingBST_netTheta_p97.5")
netTheta_table_changingBST$amazonID = NA
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    netTheta_table_changingBST$amazonID[which(netTheta_table_changingBST$subject==i)] = testBehaviourData$amazonID[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)][1]
}

write.csv(netTheta_table_changingBST,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/netTheta_table_changingBST.csv",
            row.names=FALSE)

################################################################################
##
##                         FITTING: model_changingBT
##
################################################################################

    ## FITTING
model_changingBT = stan_model(file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model_changingBT.stan') # debug
fit_changingBT = sampling(
    model_changingBT, data=ReinforcementLearningStanData_group, seed=77,
    #control = list(adapt_delta = 0.9999, max_treedepth =20, stepsize = 0.001),
    pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_annealing','s_annealing','mu_soc','s_soc','mu_theta','s_theta','mu_theta_slope','s_theta_slope','alpha','beta','annealing','theta','theta_slope','soc','netBeta','netTheta','log_lik'),
    init=function() {
        list(
            mu_beta = runif(3,0,4),# because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            mu_annealing = runif(3,0,4),# because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_alpha = runif(3,0,3),
            s_beta = runif(3,0,3), # because exp(beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_annealing = runif(3,0,3), # because exp(beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_alpha = runif(3,0,3),
            s_theta = runif(3,0,3),
            s_theta_slope = runif(3,0,3),
            s_soc = runif(3,0,3)
            #s_soc_reduction = runif(3,0,3)#rep(2,3)#runif(1, 1, 4)
            )
    },
    #chains=8, iter=1200, warmup=600, thin=5
    chains=8, iter=2000, warmup=1000, thin=5
    #chains=8, iter=120, warmup=60, thin=1
    #sample_file = "stan_output/fit_changingBT_"
    #sample_file = "/Volumes/LaCie/stan_fit/st-andrews-onlineexp2/fit_changingBT_"
)
fit_changingBT
saveRDS(fit_changingBT, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/fit_changingBT.rds')
ms_changingBT <- rstan::extract(fit_changingBT)
plot(fit_changingBT, pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_annealing','s_annealing','mu_soc','s_soc','mu_theta','s_theta','mu_theta_slope','s_theta_slope'))
#plot(fit_changingBT, pars=c('log_lik'))
traceplot(fit_changingBT, pars=c('mu_alpha','mu_beta','mu_annealing'), inc_warmup=TRUE)
traceplot(fit_changingBT, pars=c('mu_theta','mu_theta_slope','mu_soc'), inc_warmup=TRUE)
traceplot(fit_changingBT, pars=c('s_alpha','s_beta','s_annealing','s_theta','s_theta_slope','s_soc'), inc_warmup=TRUE)



temp_changingBT = data.frame(WAICi_changingBT = WAIC_indv(fit_changingBT)$waic)

temp_changingBT <- cbind(temp_changingBT, as.data.frame(t(apply(ms_changingBT$alpha[,1:ncol(ms_changingBT$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBT <- cbind(temp_changingBT, as.data.frame(t(apply(ms_changingBT$beta[,1:ncol(ms_changingBT$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBT <- cbind(temp_changingBT, as.data.frame(t(apply(ms_changingBT$annealing[,1:ncol(ms_changingBT$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBT <- cbind(temp_changingBT, as.data.frame(t(apply(ms_changingBT$theta[,1:ncol(ms_changingBT$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
#temp_changingBT <- cbind(temp_changingBT, as.data.frame(t(apply(ms_changingBT$soc_reduction[,1:ncol(ms_changingBT$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBT <- cbind(temp_changingBT, as.data.frame(t(apply(ms_changingBT$soc[,1:ncol(ms_changingBT$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(temp_changingBT) <- c("changingBTi",
                            "changingBT_alpha_p2.5", "changingBT_alpha_p25", "changingBT_alpha_p50", "changingBT_alpha_p75", "changingBT_alpha_p97.5",
                            "changingBT_beta_p2.5", "changingBT_beta_p25", "changingBT_beta_p50", "changingBT_beta_p75", "changingBT_beta_p97.5",
                            "changingBT_annealing_p2.5", "changingBT_annealing_p25", "changingBT_annealing_p50", "changingBT_annealing_p75", "changingBT_annealing_p97.5",
                            "changingBT_theta_p2.5", "changingBT_theta_p25", "changingBT_theta_p50", "changingBT_theta_p75", "changingBT_theta_p97.5",
                            #"changingBT_soc_reduction_p2.5", "changingBT_soc_reduction_p25", "changingBT_soc_reduction_p50", "changingBT_soc_reduction_p75", "changingBT_soc_reduction_p97.5",
                            "changingBT_soc_p2.5", "changingBT_soc_p25", "changingBT_soc_p50", "changingBT_soc_p75", "changingBT_soc_p97.5")
max_log_lik = apply(ms_changingBT$log_lik[,1:ncol(ms_changingBT$log_lik)], 2, which.max)
max_changingBT_alpha=rep(NA, ncol(ms_changingBT$log_lik))
max_changingBT_beta=rep(NA, ncol(ms_changingBT$log_lik))
max_changingBT_annealing=rep(NA, ncol(ms_changingBT$log_lik))
max_changingBT_theta=rep(NA, ncol(ms_changingBT$log_lik))
#max_changingBT_soc_reduction=rep(NA, ncol(ms_changingBT$log_lik))
max_changingBT_soc=rep(NA, ncol(ms_changingBT$log_lik))
for(i in 1:ncol(ms_changingBT$log_lik)){
    max_changingBT_alpha[i] <- ms_changingBT$alpha[max_log_lik[i],i]
    max_changingBT_beta[i] <- ms_changingBT$beta[max_log_lik[i],i]
    max_changingBT_annealing[i] <- ms_changingBT$annealing[max_log_lik[i],i]
    max_changingBT_theta[i] <- ms_changingBT$theta[max_log_lik[i],i]
    #max_changingBT_soc_reduction[i] <- ms_changingBT$soc_reduction[max_log_lik[i],i]
    max_changingBT_soc[i] <- ms_changingBT$soc[max_log_lik[i],i]
}
temp_changingBT$max_changingBT_alpha = max_changingBT_alpha
temp_changingBT$max_changingBT_beta = max_changingBT_beta
temp_changingBT$max_changingBT_annealing = max_changingBT_annealing
temp_changingBT$max_changingBT_theta = max_changingBT_theta
#temp_changingBT$max_changingBT_soc_reduction = max_changingBT_soc_reduction
temp_changingBT$max_changingBT_soc = max_changingBT_soc
temp_changingBT$max_changingBT_log_lik = apply(ms_changingBT$log_lik[,1:ncol(ms_changingBT$log_lik)], 2, max)
numParam = 5
temp_changingBT$AIC_changingBT = -temp_changingBT$max_changingBT_log_lik + numParam ## 正確には AIC/2
temp_changingBT$BIC_changingBT = -temp_changingBT$max_changingBT_log_lik + (numParam/2)*log(70) ## more precisely, BIC/2
temp_changingBT$WAICi_changingBT = WAIC_indv(fit_changingBT)$waic

temp_changingBT$amazonID = NA
temp_changingBT$taskDifficulty = NA
temp_changingBT$totalRound = NA
temp_changingBT$correctNum = NA
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    temp_changingBT$amazonID[i] = testBehaviourData$amazonID[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)][1]
    temp_changingBT$taskDifficulty[i] = names(table(testBehaviourData$taskDifficulty[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)]))
    temp_changingBT$totalRound[i] = length(which(as.numeric(as.factor(subset(testBehaviourData, machine>=0)$amazonID))==i))
    temp_changingBT$correctNum[i] = sum(testBehaviourData$optimalChoice[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)], na.rm=TRUE)
}
temp_changingBT$correctProp = temp_changingBT$correctNum/temp_changingBT$totalRound

temp_changingBT$log_lik_n_eff = NA
temp_changingBT$log_lik_Rhat = NA
temp_changingBT$alpha_n_eff = NA
temp_changingBT$alpha_Rhat = NA
temp_changingBT$beta_n_eff = NA
temp_changingBT$beta_Rhat = NA
temp_changingBT$soc_n_eff = NA
temp_changingBT$soc_Rhat = NA
#temp_changingBT$soc_reduction_n_eff = NA
#temp_changingBT$soc_reduction_Rhat = NA
temp_changingBT$theta_n_eff = NA
temp_changingBT$theta_Rhat = NA
temp_changingBT$annealing_n_eff = NA
temp_changingBT$annealing_Rhat = NA
for(i in 1:nrow(temp_changingBT)){
    temp_changingBT$log_lik_n_eff[i] <- summary(fit_changingBT)$summary[paste('log_lik[',i,']',sep=''),'n_eff']
    temp_changingBT$log_lik_Rhat[i] <- summary(fit_changingBT)$summary[paste('log_lik[',i,']',sep=''),'Rhat']
    temp_changingBT$alpha_n_eff[i] <- summary(fit_changingBT)$summary[paste('alpha[',i,']',sep=''),'n_eff']
    temp_changingBT$alpha_Rhat[i] <- summary(fit_changingBT)$summary[paste('alpha[',i,']',sep=''),'Rhat']
    temp_changingBT$beta_n_eff[i] <- summary(fit_changingBT)$summary[paste('beta[',i,']',sep=''),'n_eff']
    temp_changingBT$beta_Rhat[i] <- summary(fit_changingBT)$summary[paste('beta[',i,']',sep=''),'Rhat']
    temp_changingBT$soc_n_eff[i] <- summary(fit_changingBT)$summary[paste('soc[',i,']',sep=''),'n_eff']
    temp_changingBT$soc_Rhat[i] <- summary(fit_changingBT)$summary[paste('soc[',i,']',sep=''),'Rhat']
    #temp_changingBT$soc_reduction_n_eff[i] <- summary(fit_changingBT)$summary[paste('soc_reduction[',i,']',sep=''),'n_eff']
    #temp_changingBT$soc_reduction_Rhat[i] <- summary(fit_changingBT)$summary[paste('soc_reduction[',i,']',sep=''),'Rhat']
    temp_changingBT$theta_n_eff[i] <- summary(fit_changingBT)$summary[paste('theta[',i,']',sep=''),'n_eff']
    temp_changingBT$theta_Rhat[i] <- summary(fit_changingBT)$summary[paste('theta[',i,']',sep=''),'Rhat']
    temp_changingBT$annealing_n_eff[i] <- summary(fit_changingBT)$summary[paste('annealing[',i,']',sep=''),'n_eff']
    temp_changingBT$annealing_Rhat[i] <- summary(fit_changingBT)$summary[paste('annealing[',i,']',sep=''),'Rhat']
}

write.csv(temp_changingBT,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/temp_changingBT.csv",
            #"/Users/watall/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/temp_changingBT.csv",
            row.names=FALSE)


    ## Estimated netBeta at each round
netBeta_table_changingBT = data.frame(subject = 1:ncol(ms_changingBT$log_lik), round = rep(1, ncol(ms_changingBT$log_lik)))
netBeta_table_changingBT = cbind(netBeta_table_changingBT, as.data.frame(t(apply(ms_changingBT$netBeta[,1:ncol(ms_changingBT$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(netBeta_table_changingBT) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
    netBeta_table_changingBT = rbind(netBeta_table_changingBT, cbind(cbind(1:ncol(ms_changingBT$log_lik), rep(t,ncol(ms_changingBT$log_lik))), as.data.frame(t(apply(ms_changingBT$netBeta[,1:ncol(ms_changingBT$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(netBeta_table_changingBT) <- c("subject", "round", "changingBT_netBeta_p2.5", "changingBT_netBeta_p25", "changingBT_netBeta_p50", "changingBT_netBeta_p75", "changingBT_netBeta_p97.5")
netBeta_table_changingBT$amazonID = NA
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    netBeta_table_changingBT$amazonID[which(netBeta_table_changingBT$subject==i)] = testBehaviourData$amazonID[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)][1]
}

write.csv(netBeta_table_changingBT,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/netBeta_table_changingBT.csv",
            #"/Users/watall/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/netBeta_table_changingBT.csv",
            row.names=FALSE)



