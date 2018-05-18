# RESCOLA-WAGNER + SOFTMAX-ANNEALING CHOICE MODEL (model_UNC6_sReduc_annealing.stan)
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
#testNames = head(names(table(subset(allBehaviouralData_all, taskDifficulty_num==1&sizeCategory>5)$amazonID)),5)
#testNames <- append(testNames, head(names(table(subset(allBehaviouralData_all, taskDifficulty_num==2&sizeCategory>5)$amazonID)), 5))
#testNames <- append(testNames, head(names(table(subset(allBehaviouralData_all, taskDifficulty_num==3&sizeCategory>5)$amazonID)), 5))
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
##                   FITTING: UNC6_sReduc_annealing model
##
################################################################################

    ## FITTING
model_UNC6_sReduc_annealing = stan_model(file="~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model_UNC6_sReduc_annealing.stan") # debug
fit_UNC6_sReduc_annealing = sampling(
    model_UNC6_sReduc_annealing, data=ReinforcementLearningStanData_group, seed=69,
    pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_soc','s_soc','mu_theta','s_theta','mu_soc_reduction','s_soc_reduction','mu_annealing','s_annealing','alpha','beta','annealing','theta','soc_raw','soc_reduction','soc','netBeta','log_lik'),
    init=function() {
        list(
            mu_beta = runif(3,0,4),# because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            mu_annealing = runif(3,0,4),# because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_alpha = runif(3,1,3),
            s_beta = runif(3,1,3), # because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_annealing = runif(3,1,3), # because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_alpha = runif(3,1,3),
            s_theta = runif(3,1,3),
            s_soc = runif(3,1,3),
            s_soc_reduction = runif(3,1,3)#rep(2,3)#runif(1, 1, 4)
            )
    },
    chains=8, iter=2000, warmup=1000, thin=5
    #chains=8, iter=1200, warmup=600, thin=3
)
fit_UNC6_sReduc_annealing
saveRDS(fit_UNC6_sReduc_annealing, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/fit_UNC6_sReduc_annealing.rds')

ms_UNC6_sReduc_annealing <- rstan::extract(fit_UNC6_sReduc_annealing)
mcmc_results_plot = plot(fit_UNC6_sReduc_annealing, pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_soc','s_soc','mu_theta','s_theta','mu_soc_reduction','s_soc_reduction','mu_annealing','s_annealing'))
trace_1 = traceplot(fit_UNC6_sReduc_annealing, pars=c('mu_alpha','mu_beta','mu_annealing'), inc_warmup=FALSE)
trace_2 = traceplot(fit_UNC6_sReduc_annealing, pars=c('mu_theta','mu_soc','mu_soc_reduction'), inc_warmup=FALSE)
trace_3 = traceplot(fit_UNC6_sReduc_annealing, pars=c('s_alpha','s_beta','s_theta','s_soc','s_soc_reduction'), inc_warmup=FALSE)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/mcmc_results_plot.pdf", plot = mcmc_results_plot, dpi = 600, width = 9, height = 9)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/trace_1.pdf", plot = trace_1, dpi = 600, width = 9, height = 9)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/trace_2.pdf", plot = trace_2, dpi = 600, width = 9, height = 9)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/trace_3.pdf", plot = trace_3, dpi = 600, width = 9, height = 9)


temp_UNC6_sReduc_annealing = data.frame(WAICi_UNC6_sReduc_annealing = WAIC_indv(fit_UNC6_sReduc_annealing)$waic)

temp_UNC6_sReduc_annealing <- cbind(temp_UNC6_sReduc_annealing, as.data.frame(t(apply(ms_UNC6_sReduc_annealing$alpha[,1:ncol(ms_UNC6_sReduc_annealing$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing <- cbind(temp_UNC6_sReduc_annealing, as.data.frame(t(apply(ms_UNC6_sReduc_annealing$beta[,1:ncol(ms_UNC6_sReduc_annealing$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing <- cbind(temp_UNC6_sReduc_annealing, as.data.frame(t(apply(ms_UNC6_sReduc_annealing$annealing[,1:ncol(ms_UNC6_sReduc_annealing$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing <- cbind(temp_UNC6_sReduc_annealing, as.data.frame(t(apply(ms_UNC6_sReduc_annealing$theta[,1:ncol(ms_UNC6_sReduc_annealing$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing <- cbind(temp_UNC6_sReduc_annealing, as.data.frame(t(apply(ms_UNC6_sReduc_annealing$soc_reduction[,1:ncol(ms_UNC6_sReduc_annealing$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing <- cbind(temp_UNC6_sReduc_annealing, as.data.frame(t(apply(ms_UNC6_sReduc_annealing$soc_raw[,1:ncol(ms_UNC6_sReduc_annealing$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(temp_UNC6_sReduc_annealing) <- c("UNC6_sReduc_annealingi",
                            "UNC6_sReduc_annealing_alpha_p2.5", "UNC6_sReduc_annealing_alpha_p25", "UNC6_sReduc_annealing_alpha_p50", "UNC6_sReduc_annealing_alpha_p75", "UNC6_sReduc_annealing_alpha_p97.5",
                            "UNC6_sReduc_annealing_beta_p2.5", "UNC6_sReduc_annealing_beta_p25", "UNC6_sReduc_annealing_beta_p50", "UNC6_sReduc_annealing_beta_p75", "UNC6_sReduc_annealing_beta_p97.5",
                            "UNC6_sReduc_annealing_annealing_p2.5", "UNC6_sReduc_annealing_annealing_p25", "UNC6_sReduc_annealing_annealing_p50", "UNC6_sReduc_annealing_annealing_p75", "UNC6_sReduc_annealing_annealing_p97.5",
                            "UNC6_sReduc_annealing_theta_p2.5", "UNC6_sReduc_annealing_theta_p25", "UNC6_sReduc_annealing_theta_p50", "UNC6_sReduc_annealing_theta_p75", "UNC6_sReduc_annealing_theta_p97.5",
                            "UNC6_sReduc_annealing_soc_reduction_p2.5", "UNC6_sReduc_annealing_soc_reduction_p25", "UNC6_sReduc_annealing_soc_reduction_p50", "UNC6_sReduc_annealing_soc_reduction_p75", "UNC6_sReduc_annealing_soc_reduction_p97.5",
                            "UNC6_sReduc_annealing_soc_raw_p2.5", "UNC6_sReduc_annealing_soc_raw_p25", "UNC6_sReduc_annealing_soc_raw_p50", "UNC6_sReduc_annealing_soc_raw_p75", "UNC6_sReduc_annealing_soc_raw_p97.5")
max_log_lik = apply(ms_UNC6_sReduc_annealing$log_lik[,1:ncol(ms_UNC6_sReduc_annealing$log_lik)], 2, which.max)
max_UNC6_sReduc_annealing_alpha=rep(NA, ncol(ms_UNC6_sReduc_annealing$log_lik))
max_UNC6_sReduc_annealing_beta=rep(NA, ncol(ms_UNC6_sReduc_annealing$log_lik))
max_UNC6_sReduc_annealing_annealing=rep(NA, ncol(ms_UNC6_sReduc_annealing$log_lik))
max_UNC6_sReduc_annealing_theta=rep(NA, ncol(ms_UNC6_sReduc_annealing$log_lik))
max_UNC6_sReduc_annealing_soc_reduction=rep(NA, ncol(ms_UNC6_sReduc_annealing$log_lik))
max_UNC6_sReduc_annealing_soc_raw=rep(NA, ncol(ms_UNC6_sReduc_annealing$log_lik))
for(i in 1:ncol(ms_UNC6_sReduc_annealing$log_lik)){
    max_UNC6_sReduc_annealing_alpha[i] <- ms_UNC6_sReduc_annealing$alpha[max_log_lik[i],i]
    max_UNC6_sReduc_annealing_beta[i] <- ms_UNC6_sReduc_annealing$beta[max_log_lik[i],i]
    max_UNC6_sReduc_annealing_annealing[i] <- ms_UNC6_sReduc_annealing$annealing[max_log_lik[i],i]
    max_UNC6_sReduc_annealing_theta[i] <- ms_UNC6_sReduc_annealing$theta[max_log_lik[i],i]
    max_UNC6_sReduc_annealing_soc_reduction[i] <- ms_UNC6_sReduc_annealing$soc_reduction[max_log_lik[i],i]
    max_UNC6_sReduc_annealing_soc_raw[i] <- ms_UNC6_sReduc_annealing$soc_raw[max_log_lik[i],i]
}
temp_UNC6_sReduc_annealing$max_UNC6_sReduc_annealing_alpha = max_UNC6_sReduc_annealing_alpha
temp_UNC6_sReduc_annealing$max_UNC6_sReduc_annealing_beta = max_UNC6_sReduc_annealing_beta
temp_UNC6_sReduc_annealing$max_UNC6_sReduc_annealing_annealing = max_UNC6_sReduc_annealing_annealing
temp_UNC6_sReduc_annealing$max_UNC6_sReduc_annealing_theta = max_UNC6_sReduc_annealing_theta
temp_UNC6_sReduc_annealing$max_UNC6_sReduc_annealing_soc_reduction = max_UNC6_sReduc_annealing_soc_reduction
temp_UNC6_sReduc_annealing$max_UNC6_sReduc_annealing_soc_raw = max_UNC6_sReduc_annealing_soc_raw
temp_UNC6_sReduc_annealing$max_UNC6_sReduc_annealing_log_lik = apply(ms_UNC6_sReduc_annealing$log_lik[,1:ncol(ms_UNC6_sReduc_annealing$log_lik)], 2, max)
numParam = 6
temp_UNC6_sReduc_annealing$AIC_UNC6_sReduc_annealing = -temp_UNC6_sReduc_annealing$max_UNC6_sReduc_annealing_log_lik + numParam ## 正確には AIC/2
temp_UNC6_sReduc_annealing$BIC_UNC6_sReduc_annealing = -temp_UNC6_sReduc_annealing$max_UNC6_sReduc_annealing_log_lik + (numParam/2)*log(70) ## more precisely, BIC/2
temp_UNC6_sReduc_annealing$WAICi_UNC6_sReduc_annealing = WAIC_indv(fit_UNC6_sReduc_annealing)$waic

temp_UNC6_sReduc_annealing$amazonID = NA
temp_UNC6_sReduc_annealing$taskDifficulty = NA
temp_UNC6_sReduc_annealing$totalRound = NA
temp_UNC6_sReduc_annealing$correctNum = NA
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    temp_UNC6_sReduc_annealing$amazonID[i] = testBehaviourData$amazonID[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)][1]
    temp_UNC6_sReduc_annealing$taskDifficulty[i] = names(table(testBehaviourData$taskDifficulty[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)]))
    temp_UNC6_sReduc_annealing$totalRound[i] = length(which(as.numeric(as.factor(subset(testBehaviourData, machine>=0)$amazonID))==i))
    temp_UNC6_sReduc_annealing$correctNum[i] = sum(testBehaviourData$optimalChoice[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)], na.rm=TRUE)
}
temp_UNC6_sReduc_annealing$correctProp = temp_UNC6_sReduc_annealing$correctNum/temp_UNC6_sReduc_annealing$totalRound

temp_UNC6_sReduc_annealing$log_lik_n_eff = NA
temp_UNC6_sReduc_annealing$log_lik_Rhat = NA
temp_UNC6_sReduc_annealing$alpha_n_eff = NA
temp_UNC6_sReduc_annealing$alpha_Rhat = NA
temp_UNC6_sReduc_annealing$beta_n_eff = NA
temp_UNC6_sReduc_annealing$beta_Rhat = NA
temp_UNC6_sReduc_annealing$soc_raw_n_eff = NA
temp_UNC6_sReduc_annealing$soc_raw_Rhat = NA
temp_UNC6_sReduc_annealing$soc_reduction_n_eff = NA
temp_UNC6_sReduc_annealing$soc_reduction_Rhat = NA
temp_UNC6_sReduc_annealing$theta_n_eff = NA
temp_UNC6_sReduc_annealing$theta_Rhat = NA
temp_UNC6_sReduc_annealing$annealing_n_eff = NA
temp_UNC6_sReduc_annealing$annealing_Rhat = NA
for(i in 1:nrow(temp_UNC6_sReduc_annealing)){
    temp_UNC6_sReduc_annealing$log_lik_n_eff[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('log_lik[',i,']',sep=''),'n_eff']
    temp_UNC6_sReduc_annealing$log_lik_Rhat[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('log_lik[',i,']',sep=''),'Rhat']
    temp_UNC6_sReduc_annealing$alpha_n_eff[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('alpha[',i,']',sep=''),'n_eff']
    temp_UNC6_sReduc_annealing$alpha_Rhat[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('alpha[',i,']',sep=''),'Rhat']
    temp_UNC6_sReduc_annealing$beta_n_eff[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('beta[',i,']',sep=''),'n_eff']
    temp_UNC6_sReduc_annealing$beta_Rhat[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('beta[',i,']',sep=''),'Rhat']
    temp_UNC6_sReduc_annealing$soc_raw_n_eff[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('soc_raw[',i,']',sep=''),'n_eff']
    temp_UNC6_sReduc_annealing$soc_raw_Rhat[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('soc_raw[',i,']',sep=''),'Rhat']
    temp_UNC6_sReduc_annealing$soc_reduction_n_eff[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('soc_reduction[',i,']',sep=''),'n_eff']
    temp_UNC6_sReduc_annealing$soc_reduction_Rhat[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('soc_reduction[',i,']',sep=''),'Rhat']
    temp_UNC6_sReduc_annealing$theta_n_eff[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('theta[',i,']',sep=''),'n_eff']
    temp_UNC6_sReduc_annealing$theta_Rhat[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('theta[',i,']',sep=''),'Rhat']
    temp_UNC6_sReduc_annealing$annealing_n_eff[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('annealing[',i,']',sep=''),'n_eff']
    temp_UNC6_sReduc_annealing$annealing_Rhat[i] <- summary(fit_UNC6_sReduc_annealing)$summary[paste('annealing[',i,']',sep=''),'Rhat']
}

write.csv(temp_UNC6_sReduc_annealing,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/temp_UNC6_sReduc_annealing.csv",
            row.names=FALSE)

    ## Estimated copying rate at each round
soc_table_UNC6_sReduc_annealing = data.frame(subject = 1:ncol(ms_UNC6_sReduc_annealing$log_lik), round = rep(1, ncol(ms_UNC6_sReduc_annealing$log_lik)))
soc_table_UNC6_sReduc_annealing = cbind(soc_table_UNC6_sReduc_annealing, as.data.frame(t(apply(ms_UNC6_sReduc_annealing$soc[,1:ncol(ms_UNC6_sReduc_annealing$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(soc_table_UNC6_sReduc_annealing) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
    soc_table_UNC6_sReduc_annealing = rbind(soc_table_UNC6_sReduc_annealing, cbind(cbind(1:ncol(ms_UNC6_sReduc_annealing$log_lik), rep(t,ncol(ms_UNC6_sReduc_annealing$log_lik))), as.data.frame(t(apply(ms_UNC6_sReduc_annealing$soc[,1:ncol(ms_UNC6_sReduc_annealing$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(soc_table_UNC6_sReduc_annealing) <- c("subject", "round", "UNC6_sReduc_annealing_soc_p2.5", "UNC6_sReduc_annealing_soc_p25", "UNC6_sReduc_annealing_soc_p50", "UNC6_sReduc_annealing_soc_p75", "UNC6_sReduc_annealing_soc_p97.5")
soc_table_UNC6_sReduc_annealing$amazonID = NA
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    soc_table_UNC6_sReduc_annealing$amazonID[which(soc_table_UNC6_sReduc_annealing$subject==i)] = testBehaviourData$amazonID[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)][1]
}

write.csv(soc_table_UNC6_sReduc_annealing,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/soc_table_UNC6_sReduc_annealing.csv",
            row.names=FALSE)

    ## Estimated netBeta at each round
netBeta_table_UNC6_sReduc_annealing = data.frame(subject = 1:ncol(ms_UNC6_sReduc_annealing$log_lik), round = rep(1, ncol(ms_UNC6_sReduc_annealing$log_lik)))
netBeta_table_UNC6_sReduc_annealing = cbind(netBeta_table_UNC6_sReduc_annealing, as.data.frame(t(apply(ms_UNC6_sReduc_annealing$netBeta[,1:ncol(ms_UNC6_sReduc_annealing$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(netBeta_table_UNC6_sReduc_annealing) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
    netBeta_table_UNC6_sReduc_annealing = rbind(netBeta_table_UNC6_sReduc_annealing, cbind(cbind(1:ncol(ms_UNC6_sReduc_annealing$log_lik), rep(t,ncol(ms_UNC6_sReduc_annealing$log_lik))), as.data.frame(t(apply(ms_UNC6_sReduc_annealing$netBeta[,1:ncol(ms_UNC6_sReduc_annealing$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(netBeta_table_UNC6_sReduc_annealing) <- c("subject", "round", "UNC6_sReduc_annealing_netBeta_p2.5", "UNC6_sReduc_annealing_netBeta_p25", "UNC6_sReduc_annealing_netBeta_p50", "UNC6_sReduc_annealing_netBeta_p75", "UNC6_sReduc_annealing_netBeta_p97.5")
netBeta_table_UNC6_sReduc_annealing$amazonID = NA
for(i in 1:ReinforcementLearningStanData_group$Nsub){
    netBeta_table_UNC6_sReduc_annealing$amazonID[which(netBeta_table_UNC6_sReduc_annealing$subject==i)] = testBehaviourData$amazonID[which(as.numeric(as.factor(testBehaviourData$amazonID))==i)][1]
}

write.csv(netBeta_table_UNC6_sReduc_annealing,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/netBeta_table_UNC6_sReduc_annealing.csv",
            row.names=FALSE)



################################################################################
##
##                   FITTING: AL6_annealing model for Solitary condition
##
################################################################################


# Making ReinforcementLearningStanData_indiv
testBehaviourData_indiv = subset(allBehaviouralData_all, sizeCategory == 1) # N >= 2
# Using only the subjects who played the task at least for 31 rounds
testBehaviourData_indiv = subset(testBehaviourData_indiv, amazonID %in% names(which(table(testBehaviourData_indiv$amazonID)>30)))

# For debug #
#testNames = head(names(table(subset(allBehaviouralData_all, taskDifficulty_num==1&sizeCategory==1)$amazonID)),5)
#testNames <- append(testNames, head(names(table(subset(allBehaviouralData_all, taskDifficulty_num==2&sizeCategory==1)$amazonID)), 5))
#testNames <- append(testNames, head(names(table(subset(allBehaviouralData_all, taskDifficulty_num==3&sizeCategory==1)$amazonID)), 5))
#testBehaviourData_indiv = subset(allBehaviouralData_all, amazonID %in% testNames)
# end - debug #

testBehaviourData_indiv$sub = as.numeric(as.factor(testBehaviourData_indiv$amazonID))

    # insert missed trials
for(i in 1:length(table(testBehaviourData_indiv$amazonID))) {
    thisSubject = subset(testBehaviourData_indiv,sub==i)
    for(t in 1:70) {
        if(nrow(thisSubject[which(thisSubject$round==t),])==0) {
            testBehaviourData_indiv <- rbind(testBehaviourData_indiv, rep(NA, ncol(testBehaviourData_indiv)))
            testBehaviourData_indiv$machine[nrow(testBehaviourData_indiv)] = -1
            testBehaviourData_indiv$result[nrow(testBehaviourData_indiv)] = -1
            testBehaviourData_indiv$round[nrow(testBehaviourData_indiv)] = t
            testBehaviourData_indiv$amazonID[nrow(testBehaviourData_indiv)] = testBehaviourData_indiv$amazonID[which(testBehaviourData_indiv$sub==i)][1]
            testBehaviourData_indiv$sub[nrow(testBehaviourData_indiv)] = testBehaviourData_indiv$sub[which(testBehaviourData_indiv$sub==i)][1]
            testBehaviourData_indiv$taskDifficulty_num[nrow(testBehaviourData_indiv)] = testBehaviourData_indiv$taskDifficulty_num[which(testBehaviourData_indiv$sub==i)][1]
        }
    }
}
testBehaviourData_indiv = testBehaviourData_indiv[order(testBehaviourData_indiv$round),] # Sorting by round

ReinforcementLearningStanData_indiv = list(
    All = nrow(testBehaviourData_indiv),
    Nsub = length(table(testBehaviourData_indiv$amazonID)),
    Ncue = 3, # number of options
    Ntrial = 70,
    Ncon = 3,
    sub = testBehaviourData_indiv$sub,
    Y = testBehaviourData_indiv$machine + 1,
    trial = testBehaviourData_indiv$round,
    outcome = testBehaviourData_indiv$result * 100
    )
con = numeric(length(table(testBehaviourData_indiv$amazonID)))
for(i in 1:ReinforcementLearningStanData_indiv$Nsub){
    con[i] <- testBehaviourData_indiv$taskDifficulty_num[which(testBehaviourData_indiv$sub==i)][1]
}
ReinforcementLearningStanData_indiv$condition = con
startT = c()
MaxTrial = c()
for(i in 1:length(table(testBehaviourData_indiv$amazonID))){
    startT = append(startT, min(subset(testBehaviourData_indiv, sub==i)$round))
    MaxTrial = append(MaxTrial, length(which(subset(testBehaviourData_indiv, sub==i)$machine>=0)))
}
ReinforcementLearningStanData_indiv$startT = startT
ReinforcementLearningStanData_indiv$MaxTrial = MaxTrial


    ## FITTING
model_AL6_annealing = stan_model(file="~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model_AL6_annealing.stan") # debug
fit_AL6_annealing_indiv = sampling(
    model_AL6_annealing, data=ReinforcementLearningStanData_indiv, seed=69,
    #control = list(adapt_delta = 0.90, max_treedepth=13),
    pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_annealing','s_annealing','alpha','beta','annealing','netBeta','log_lik'),
    init=function() {
        list(
            mu_beta = runif(3,0,4),# because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            mu_annealing = runif(3,0,4),# because (beta + annealing*t/70), both beta and annealing shouldn't be too large!
            s_alpha = runif(3,0,3),#rep(2,3),
            s_beta = runif(3,0,3),#rep(2,3),
            s_annealing = runif(3,0,3)
            )
    },
    chains=8, iter=2000, warmup=1000, thin=5
    #chains=8, iter=1200, warmup=600, thin=3
)
fit_AL6_annealing_indiv
saveRDS(fit_AL6_annealing_indiv, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/fit_AL6_annealing_indiv.rds')

ms_AL6_annealing_indiv <- rstan::extract(fit_AL6_annealing_indiv)
mcmc_results_plot_AL = plot(fit_AL6_annealing_indiv, pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_annealing','s_annealing'))
#plot(fit_AL6_annealing_indiv, pars=c('log_lik'))

    ## 事後診断
trace_1_AL = traceplot(fit_AL6_annealing_indiv, pars=c('mu_alpha','mu_beta','mu_annealing'), inc_warmup=FALSE)
trace_2_AL = traceplot(fit_AL6_annealing_indiv, pars=c('s_alpha','s_beta','s_annealing'), inc_warmup=FALSE)

ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/mcmc_results_plot_AL.pdf", plot = mcmc_results_plot_AL, dpi = 600, width = 9, height = 9)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/trace_1_AL.pdf", plot = trace_1_AL, dpi = 600, width = 9, height = 9)
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/trace_2_AL.pdf", plot = trace_2_AL, dpi = 600, width = 9, height = 9)

AL6_annealing_indiv_mcmc <- data.frame(alpha=ms_AL6_annealing_indiv$alpha[,25], beta=ms_AL6_annealing_indiv$beta[,25], log_lik=ms_AL6_annealing_indiv$log_lik[,25])
ggplot(AL6_annealing_indiv_mcmc,aes(alpha,beta, colour=log_lik))+geom_point(alpha=1/5)+scale_color_distiller(palette="Spectral")


temp_AL6_annealing_indiv = data.frame(WAICi = WAIC_indv(fit_AL6_annealing_indiv)$waic)
temp_AL6_annealing_indiv$amazonID = NA
temp_AL6_annealing_indiv$taskDifficulty = NA
temp_AL6_annealing_indiv$totalRound = NA
temp_AL6_annealing_indiv$correctNum = NA
for(i in 1:ReinforcementLearningStanData_indiv$Nsub){
    temp_AL6_annealing_indiv$amazonID[i] = testBehaviourData_indiv$amazonID[which(as.numeric(as.factor(testBehaviourData_indiv$amazonID))==i)][1]
    temp_AL6_annealing_indiv$taskDifficulty[i] = names(table(testBehaviourData_indiv$taskDifficulty[which(as.numeric(as.factor(testBehaviourData_indiv$amazonID))==i)]))
    temp_AL6_annealing_indiv$totalRound[i] = length(which(as.numeric(as.factor(subset(testBehaviourData_indiv, machine>=0)$amazonID))==i))
    temp_AL6_annealing_indiv$correctNum[i] = sum(testBehaviourData_indiv$optimalChoice[which(as.numeric(as.factor(testBehaviourData_indiv$amazonID))==i)], na.rm=TRUE)
}
temp_AL6_annealing_indiv$correctProp = temp_AL6_annealing_indiv$correctNum/temp_AL6_annealing_indiv$totalRound

temp_AL6_annealing_indiv <- cbind(temp_AL6_annealing_indiv, as.data.frame(t(apply(ms_AL6_annealing_indiv$log_lik[,1:ncol(ms_AL6_annealing_indiv$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_AL6_annealing_indiv <- cbind(temp_AL6_annealing_indiv, as.data.frame(t(apply(ms_AL6_annealing_indiv$alpha[,1:ncol(ms_AL6_annealing_indiv$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_AL6_annealing_indiv <- cbind(temp_AL6_annealing_indiv, as.data.frame(t(apply(ms_AL6_annealing_indiv$beta[,1:ncol(ms_AL6_annealing_indiv$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_AL6_annealing_indiv <- cbind(temp_AL6_annealing_indiv, as.data.frame(t(apply(ms_AL6_annealing_indiv$annealing[,1:ncol(ms_AL6_annealing_indiv$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(temp_AL6_annealing_indiv) <- c("WAICi","amazonID","taskDifficulty","totalRound","correctNum","correctProp",
                            "AL6_annealing_indiv_log_lik_p2.5", "AL6_annealing_indiv_log_lik_p25", "AL6_annealing_indiv_log_lik_p50", "AL6_annealing_indiv_log_lik_p75", "AL6_annealing_indiv_log_lik_p97.5",
                            "AL6_annealing_indiv_alpha_p2.5", "AL6_annealing_indiv_alpha_p25", "AL6_annealing_indiv_alpha_p50", "AL6_annealing_indiv_alpha_p75", "AL6_annealing_indiv_alpha_p97.5",
                            "AL6_annealing_indiv_beta_p2.5", "AL6_annealing_indiv_beta_p25", "AL6_annealing_indiv_beta_p50", "AL6_annealing_indiv_beta_p75", "AL6_annealing_indiv_beta_p97.5",
                            "AL6_annealing_indiv_annealing_p2.5", "AL6_annealing_indiv_annealing_p25", "AL6_annealing_indiv_annealing_p50", "AL6_annealing_indiv_annealing_p75", "AL6_annealing_indiv_annealing_p97.5")
max_log_lik = apply(ms_AL6_annealing_indiv$log_lik[,1:ncol(ms_AL6_annealing_indiv$log_lik)], 2, which.max)
max_AL6_annealing_indiv_alpha=rep(NA, ncol(ms_AL6_annealing_indiv$log_lik))
max_AL6_annealing_indiv_beta=rep(NA, ncol(ms_AL6_annealing_indiv$log_lik))
max_AL6_annealing_indiv_annealing=rep(NA, ncol(ms_AL6_annealing_indiv$log_lik))
for(i in 1:ncol(ms_AL6_annealing_indiv$log_lik)){
    max_AL6_annealing_indiv_alpha[i] <- ms_AL6_annealing_indiv$alpha[max_log_lik[i],i]
    max_AL6_annealing_indiv_beta[i] <- ms_AL6_annealing_indiv$beta[max_log_lik[i],i]
    max_AL6_annealing_indiv_annealing[i] <- ms_AL6_annealing_indiv$annealing[max_log_lik[i],i]
}
temp_AL6_annealing_indiv$max_AL6_annealing_indiv_alpha = max_AL6_annealing_indiv_alpha
temp_AL6_annealing_indiv$max_AL6_annealing_indiv_beta = max_AL6_annealing_indiv_beta
temp_AL6_annealing_indiv$max_AL6_annealing_indiv_annealing = max_AL6_annealing_indiv_annealing
temp_AL6_annealing_indiv$max_AL6_annealing_indiv_log_lik = apply(ms_AL6_annealing_indiv$log_lik[,1:ncol(ms_AL6_annealing_indiv$log_lik)], 2, max)
numParam = 3
temp_AL6_annealing_indiv$AIC_AL6_annealing_indiv = -temp_AL6_annealing_indiv$max_AL6_annealing_indiv_log_lik + numParam ## 正確には AIC/2
temp_AL6_annealing_indiv$BIC_AL6_annealing_indiv = -temp_AL6_annealing_indiv$max_AL6_annealing_indiv_log_lik + (numParam/2)*log(temp_AL6_annealing_indiv$totalRound) ## more precisely, BIC/2
temp_AL6_annealing_indiv$WAICi_AL6_annealing_indiv = WAIC_indv(fit_AL6_annealing_indiv)$waic

AL6_annealing_indiv_waic = WAIC(fit_AL6_annealing_indiv)$waic

temp_AL6_annealing_indiv$log_lik_n_eff = NA
temp_AL6_annealing_indiv$log_lik_Rhat = NA
temp_AL6_annealing_indiv$alpha_n_eff = NA
temp_AL6_annealing_indiv$alpha_Rhat = NA
temp_AL6_annealing_indiv$beta_n_eff = NA
temp_AL6_annealing_indiv$beta_Rhat = NA
temp_AL6_annealing_indiv$annealing_n_eff = NA
temp_AL6_annealing_indiv$annealing_Rhat = NA

for(i in 1:nrow(temp_AL6_annealing_indiv)){
    temp_AL6_annealing_indiv$log_lik_n_eff[i] <- summary(fit_AL6_annealing_indiv)$summary[paste('log_lik[',i,']',sep=''),'n_eff']
    temp_AL6_annealing_indiv$log_lik_Rhat[i] <- summary(fit_AL6_annealing_indiv)$summary[paste('log_lik[',i,']',sep=''),'Rhat']
    temp_AL6_annealing_indiv$alpha_n_eff[i] <- summary(fit_AL6_annealing_indiv)$summary[paste('alpha[',i,']',sep=''),'n_eff']
    temp_AL6_annealing_indiv$alpha_Rhat[i] <- summary(fit_AL6_annealing_indiv)$summary[paste('alpha[',i,']',sep=''),'Rhat']
    temp_AL6_annealing_indiv$beta_n_eff[i] <- summary(fit_AL6_annealing_indiv)$summary[paste('beta[',i,']',sep=''),'n_eff']
    temp_AL6_annealing_indiv$beta_Rhat[i] <- summary(fit_AL6_annealing_indiv)$summary[paste('beta[',i,']',sep=''),'Rhat']
    temp_AL6_annealing_indiv$annealing_n_eff[i] <- summary(fit_AL6_annealing_indiv)$summary[paste('annealing[',i,']',sep=''),'n_eff']
    temp_AL6_annealing_indiv$annealing_Rhat[i] <- summary(fit_AL6_annealing_indiv)$summary[paste('annealing[',i,']',sep=''),'Rhat']
}

write.csv(temp_AL6_annealing_indiv,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/temp_AL6_annealing_indiv.csv",
            row.names=FALSE)

    ## Estimated netBeta at each round
netBeta_table_AL6_annealing_indiv = data.frame(subject = 1:ncol(ms_AL6_annealing_indiv$log_lik), round = rep(1, ncol(ms_AL6_annealing_indiv$log_lik)))
netBeta_table_AL6_annealing_indiv = cbind(netBeta_table_AL6_annealing_indiv, as.data.frame(t(apply(ms_AL6_annealing_indiv$netBeta[,1:ncol(ms_AL6_annealing_indiv$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(netBeta_table_AL6_annealing_indiv) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
    netBeta_table_AL6_annealing_indiv = rbind(netBeta_table_AL6_annealing_indiv, cbind(cbind(1:ncol(ms_AL6_annealing_indiv$log_lik), rep(t,ncol(ms_AL6_annealing_indiv$log_lik))), as.data.frame(t(apply(ms_AL6_annealing_indiv$netBeta[,1:ncol(ms_AL6_annealing_indiv$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(netBeta_table_AL6_annealing_indiv) <- c("subject", "round", "AL6_annealing_indiv_netBeta_p2.5", "AL6_annealing_indiv_netBeta_p25", "AL6_annealing_indiv_netBeta_p50", "AL6_annealing_indiv_netBeta_p75", "AL6_annealing_indiv_netBeta_p97.5")
netBeta_table_AL6_annealing_indiv$amazonID = NA
for(i in 1:ReinforcementLearningStanData_indiv$Nsub){
    netBeta_table_AL6_annealing_indiv$amazonID[which(netBeta_table_AL6_annealing_indiv$subject==i)] = testBehaviourData_indiv$amazonID[which(as.numeric(as.factor(testBehaviourData_indiv$amazonID))==i)][1]
}

write.csv(netBeta_table_AL6_annealing_indiv,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/netBeta_table_AL6_annealing_indiv.csv",
            row.names=FALSE)













################################################################################
##
##                           Parameter Analyses
##
################################################################################




## Only considering UNC6_sReduc_annealing and AL6_annealing_indiv
temp_UNC6_sReduc_annealing <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/temp_UNC6_sReduc_annealing.csv")
temp_AL6_annealing_indiv <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/temp_AL6_annealing_indiv.csv")
soc_table_UNC6_sReduc_annealing <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/soc_table_UNC6_sReduc_annealing.csv")
netBeta_table_UNC6_sReduc_annealing <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/netBeta_table_UNC6_sReduc_annealing.csv")
netBeta_table_AL6_annealing_indiv <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/netBeta_table_AL6_annealing_indiv.csv")

performanceSummary4 = performanceSummary2
performanceSummary4$amazonID = performanceSummary4$subjectID

    ## Parameter values
performanceSummary4$alpha = NA
performanceSummary4$beta = NA
performanceSummary4$soc_raw = NA
performanceSummary4$soc_change = NA
performanceSummary4$soc2 = NA
performanceSummary4$soc70 = NA
performanceSummary4$theta = NA
performanceSummary4$low_alpha = NA
performanceSummary4$low_beta = NA
performanceSummary4$low_soc_raw = NA
performanceSummary4$low_soc_change = NA
performanceSummary4$low_soc2 = NA
performanceSummary4$low_theta = NA
performanceSummary4$high_alpha = NA
performanceSummary4$high_beta = NA
performanceSummary4$high_soc_raw = NA
performanceSummary4$high_soc_change = NA
performanceSummary4$high_soc2 = NA
performanceSummary4$high_theta = NA
performanceSummary4$alpha_25 = NA
performanceSummary4$beta_25 = NA
performanceSummary4$soc_raw_25 = NA
performanceSummary4$soc_change_25 = NA
performanceSummary4$soc2_25 = NA
performanceSummary4$theta_25 = NA
performanceSummary4$alpha_75 = NA
performanceSummary4$beta_75 = NA
performanceSummary4$soc_raw_75 = NA
performanceSummary4$soc_change_75 = NA
performanceSummary4$soc2_75 = NA
performanceSummary4$theta_75 = NA
performanceSummary4$annealing = NA
performanceSummary4$annealing_25 = NA
performanceSummary4$annealing_75 = NA
performanceSummary4$high_annealing = NA
performanceSummary4$low_annealing = NA

for(i in names(table(performanceSummary4$subjectID))){
    # If there is no subject mutching i, this means he was in the single condition
    if (length(temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_alpha_p50[which(temp_UNC6_sReduc_annealing$amazonID==i)])>0) {
        #thisSubjectModel = performanceSummary4$smallestWaicModel[which(performanceSummary4$subjectID==i)][1]
        thisSubjectRow = which(performanceSummary4$subjectID==i)
        #if(thisSubjectModel != '---') {
            # find this subject in temp file
            thisSubjectRow_in_temp = which(temp_UNC6_sReduc_annealing$amazonID==i)
            # median of posterior
            performanceSummary4$alpha[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_alpha_p50[thisSubjectRow_in_temp]
            performanceSummary4$beta[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_beta_p50[thisSubjectRow_in_temp]
            performanceSummary4$annealing[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_annealing_p50[thisSubjectRow_in_temp]
            performanceSummary4$soc_raw[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_raw_p50[thisSubjectRow_in_temp]
            performanceSummary4$soc_change[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_reduction_p50[thisSubjectRow_in_temp]
            performanceSummary4$soc2[thisSubjectRow] <- soc_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_p50[which(soc_table_UNC6_sReduc_annealing$round==2&soc_table_UNC6_sReduc_annealing$amazonID==i)]
            performanceSummary4$soc70[thisSubjectRow] <- soc_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_p50[which(soc_table_UNC6_sReduc_annealing$round==70&soc_table_UNC6_sReduc_annealing$amazonID==i)]
            performanceSummary4$theta[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_theta_p50[thisSubjectRow_in_temp]
            # low
            performanceSummary4$low_alpha[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_alpha_p2.5[thisSubjectRow_in_temp]
            performanceSummary4$low_beta[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_beta_p2.5[thisSubjectRow_in_temp]
            performanceSummary4$low_annealing[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_annealing_p2.5[thisSubjectRow_in_temp]
            performanceSummary4$low_soc_raw[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_raw_p2.5[thisSubjectRow_in_temp]
            performanceSummary4$low_soc_change[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_reduction_p2.5[thisSubjectRow_in_temp]
            performanceSummary4$low_soc2[thisSubjectRow] <- soc_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_p2.5[which(soc_table_UNC6_sReduc_annealing$round==2&soc_table_UNC6_sReduc_annealing$amazonID==i)]
            performanceSummary4$low_theta[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_theta_p2.5[thisSubjectRow_in_temp]
            # high
            performanceSummary4$high_alpha[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_alpha_p97.5[thisSubjectRow_in_temp]
            performanceSummary4$high_beta[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_beta_p97.5[thisSubjectRow_in_temp]
            performanceSummary4$high_annealing[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_annealing_p97.5[thisSubjectRow_in_temp]
            performanceSummary4$high_soc_raw[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_raw_p97.5[thisSubjectRow_in_temp]
            performanceSummary4$high_soc_change[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_reduction_p97.5[thisSubjectRow_in_temp]
            performanceSummary4$high_soc2[thisSubjectRow] <- soc_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_p97.5[which(soc_table_UNC6_sReduc_annealing$round==2&soc_table_UNC6_sReduc_annealing$amazonID==i)]
            performanceSummary4$high_theta[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_theta_p97.5[thisSubjectRow_in_temp]
            # p25
            performanceSummary4$alpha_25[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_alpha_p25[thisSubjectRow_in_temp]
            performanceSummary4$beta_25[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_beta_p25[thisSubjectRow_in_temp]
            performanceSummary4$annealing_25[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_annealing_p25[thisSubjectRow_in_temp]
            performanceSummary4$soc_raw_25[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_raw_p25[thisSubjectRow_in_temp]
            performanceSummary4$soc_change_25[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_reduction_p25[thisSubjectRow_in_temp]
            performanceSummary4$soc2_25[thisSubjectRow] <- soc_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_p25[which(soc_table_UNC6_sReduc_annealing$round==2&soc_table_UNC6_sReduc_annealing$amazonID==i)]
            performanceSummary4$theta_25[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_theta_p25[thisSubjectRow_in_temp]
            # p75
            performanceSummary4$alpha_75[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_alpha_p75[thisSubjectRow_in_temp]
            performanceSummary4$beta_75[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_beta_p75[thisSubjectRow_in_temp]
            performanceSummary4$annealing_75[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_annealing_p75[thisSubjectRow_in_temp]
            performanceSummary4$soc_raw_75[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_raw_p75[thisSubjectRow_in_temp]
            performanceSummary4$soc_change_75[thisSubjectRow] = temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_reduction_p75[thisSubjectRow_in_temp]
            performanceSummary4$soc2_75[thisSubjectRow] <- soc_table_UNC6_sReduc_annealing$UNC6_sReduc_annealing_soc_p75[which(soc_table_UNC6_sReduc_annealing$round==2&soc_table_UNC6_sReduc_annealing$amazonID==i)]
            performanceSummary4$theta_75[thisSubjectRow] <- temp_UNC6_sReduc_annealing$UNC6_sReduc_annealing_theta_p75[thisSubjectRow_in_temp]
        #}
    }
}

for(i in names(table(performanceSummary4$subjectID))){
    # If there is no subject mutching i, this means he was in the single condition
    if (length(temp_AL6_annealing_indiv$AL6_annealing_indiv_alpha_p50[which(temp_AL6_annealing_indiv$amazonID==i)])>0) {
        #thisSubjectModel = modelComparison_WAIC_indivCondition$smallestWaicModel[which(modelComparison_WAIC_indivCondition$amazonID==i)]
        #performanceSummary4$smallestWaicModel[which(performanceSummary4$subjectID==i)] = as.character(thisSubjectModel)
        thisSubjectRow = which(performanceSummary4$subjectID==i)
        #if(thisSubjectModel == 'AL_ann') {
            # find this subject in temp file
            thisSubjectRow_in_temp = which(temp_AL6_annealing_indiv$amazonID==i)
            # median of posterior
            performanceSummary4$alpha[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_alpha_p50[thisSubjectRow_in_temp]
            performanceSummary4$beta[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_beta_p50[thisSubjectRow_in_temp]
            performanceSummary4$annealing[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_annealing_p50[thisSubjectRow_in_temp]
            #performanceSummary4$soc_raw[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_raw_p50[thisSubjectRow_in_temp]
            #performanceSummary4$soc_change[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_reduction_p50[thisSubjectRow_in_temp]
            #performanceSummary4$soc2[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_p50[thisSubjectRow_in_temp]
            #performanceSummary4$soc70[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_p50[thisSubjectRow_in_temp]
            #performanceSummary4$theta[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_theta_p50[thisSubjectRow_in_temp]
            # low
            performanceSummary4$low_alpha[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_alpha_p2.5[thisSubjectRow_in_temp]
            performanceSummary4$low_beta[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_beta_p2.5[thisSubjectRow_in_temp]
            performanceSummary4$low_annealing[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_annealing_p2.5[thisSubjectRow_in_temp]
            #performanceSummary4$low_soc_raw[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_raw_p2.5[thisSubjectRow_in_temp]
            #performanceSummary4$low_soc_change[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_reduction_p2.5[thisSubjectRow_in_temp]
            #performanceSummary4$low_soc2[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_p2.5[thisSubjectRow_in_temp]
            #performanceSummary4$low_theta[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_theta_p2.5[thisSubjectRow_in_temp]
            # high
            performanceSummary4$high_alpha[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_alpha_p97.5[thisSubjectRow_in_temp]
            performanceSummary4$high_beta[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_beta_p97.5[thisSubjectRow_in_temp]
            performanceSummary4$high_annealing[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_annealing_p97.5[thisSubjectRow_in_temp]
            #performanceSummary4$high_soc_raw[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_raw_p97.5[thisSubjectRow_in_temp]
            #performanceSummary4$high_soc_change[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_reduction_p97.5[thisSubjectRow_in_temp]
            #performanceSummary4$high_soc2[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_p97.5[thisSubjectRow_in_temp]
            #performanceSummary4$high_theta[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_theta_p97.5[thisSubjectRow_in_temp]
            # p25
            performanceSummary4$alpha_25[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_alpha_p25[thisSubjectRow_in_temp]
            performanceSummary4$beta_25[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_beta_p25[thisSubjectRow_in_temp]
            performanceSummary4$annealing_25[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_annealing_p25[thisSubjectRow_in_temp]
            #performanceSummary4$soc_raw_25[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_raw_p25[thisSubjectRow_in_temp]
            #performanceSummary4$soc_change_25[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_reduction_p25[thisSubjectRow_in_temp]
            #performanceSummary4$soc2_25[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_p25[thisSubjectRow_in_temp]
            #performanceSummary4$theta_25[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_theta_p25[thisSubjectRow_in_temp]
            # p75
            performanceSummary4$alpha_75[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_alpha_p75[thisSubjectRow_in_temp]
            performanceSummary4$beta_75[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_beta_p75[thisSubjectRow_in_temp]
            performanceSummary4$annealing_75[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_annealing_p75[thisSubjectRow_in_temp]
            #performanceSummary4$soc_raw_75[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_raw_p75[thisSubjectRow_in_temp]
            #performanceSummary4$soc_change_75[thisSubjectRow] = temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_reduction_p75[thisSubjectRow_in_temp]
            #performanceSummary4$soc2_75[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_soc_p75[thisSubjectRow_in_temp]
            #performanceSummary4$theta_75[thisSubjectRow] <- temp_AL6_annealing_indiv$AL6_annealing_indiv_theta_p75[thisSubjectRow_in_temp]
        #}
    }
}

performanceSummary4$smallestWaicModel = "AL_ann"
performanceSummary4$smallestWaicModel[which(performanceSummary4$groupSize>1)] = "UNC6_sReduc_ann"

write.csv(performanceSummary4,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/performanceSummary4.csv",
            row.names=FALSE)


    ## Copying rate
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
copyingProbReduction_data_fullModelOnly$frequency_dependence_4 = 'random-copying'
copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$theta_75 < 0)] = 'neg-freq-dep'
copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$theta_25 > 0)] = 'pos-freq-dep'
#copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$best_model == 'AL' | copyingProbReduction_data_fullModelOnly$best_model == 'AL_ann')] = 'AL'
#copyingProbReduction_data_fullModelOnly$frequency_dependence_4[which(copyingProbReduction_data_fullModelOnly$best_model == 'Random')] = 'Random'
copyingProbReduction_data_fullModelOnly$frequency_dependence_4 = factor(copyingProbReduction_data_fullModelOnly$frequency_dependence_4, levels=c("Random","AL","random-copying",'neg-freq-dep', "pos-freq-dep"))

## Task difficulty 2
copyingProbReduction_data_fullModelOnly$taskDifficulty2 = NA
copyingProbReduction_data_fullModelOnly$taskDifficulty2[which(copyingProbReduction_data_fullModelOnly$taskDifficulty=='25%')] <- 'Easy'
copyingProbReduction_data_fullModelOnly$taskDifficulty2[which(copyingProbReduction_data_fullModelOnly$taskDifficulty=='50%')] <- 'Moderate'
copyingProbReduction_data_fullModelOnly$taskDifficulty2[which(copyingProbReduction_data_fullModelOnly$taskDifficulty=='78%')] <- 'Difficult'
copyingProbReduction_data_fullModelOnly$taskDifficulty2 = factor(copyingProbReduction_data_fullModelOnly$taskDifficulty2, levels=c('Easy','Moderate','Difficult'))

# Using 50% CI for categorisation
# Three categories: random-copying, neg-freq-dep, pos-freq-dep
copyingRatePlot_fullModel = ggplot() +
                    geom_line(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'), mapping=aes(round,copyingProb,group=amazonID,colour=frequency_dependence_4),alpha=2/4)+
                    stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,copyingProb),fun.y = "median",geom="line", colour = "black", size = 1)+
                    #stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,copyingProb_low),fun.y = "median",geom="line", colour = "grey50", size = 1)+
                    #stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,copyingProb_high),fun.y = "median",geom="line", colour = "grey50", size = 1)+
                    facet_grid(taskDifficulty2~frequency_dependence_4)+
                    labs(x='Round', y='Copying probability', title='') +
                    ylim(c(0,1))+
                    xlim(c(2,70))+
                    scale_colour_manual(values=c("Random"="#999999","AL"="#999999","random-copying"="#999999",'neg-freq-dep'="blue", "pos-freq-dep"="red"), name="Strength of \nfrequency dependence") +
                    #scale_colour_distiller(name="Frequency \ndependence", palette = "Spectral", direction = -1)+
                    myTheme()+
                    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/round_copyingRate_theta_fullModel.pdf", plot = copyingRatePlot_fullModel, dpi = 600, width = 9, height = 9)


## changing beta
netBetaPlot_fullModel = ggplot() +
                    geom_line(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'), mapping=aes(round,beta,group=amazonID,colour=frequency_dependence_4),alpha=2/4)+
                    stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,beta),fun.y = "median",geom="line", colour = "black", size = 1)+
                    #stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,beta_low),fun.y = "median",geom="line", colour = "grey50", size = 1)+
                    #stat_summary(data=subset(copyingProbReduction_data_fullModelOnly, best_model!='Random'&best_model!='AL'&best_model!='AL_ann'),aes(round,beta_high),fun.y = "median",geom="line", colour = "grey50", size = 1)+
                    facet_grid(taskDifficulty2~frequency_dependence_4)+
                    labs(x='Round', y='Inverse temperature', title='') +
                    #ylim(c(0,15))+
                    xlim(c(1,70))+
                    scale_colour_manual(values=c("Random"="#999999","AL"="#999999","random-copying"="#999999",'neg-freq-dep'="blue", "pos-freq-dep"="red"), name="Strength of \nfrequency dependence") +
                    #scale_colour_distiller(name="Frequency \ndependence", palette = "Spectral", direction = -1)+
                    myTheme()+
                    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)

ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/round_beta_fullModel.pdf", plot = netBetaPlot_fullModel, dpi = 600, width = 9, height = 9)


# conservative categories - rand, neg, pos (Considering Asocial model and using 50% CI bound for freq-dep category)
performanceSummary4$frequency_dependence_4 = '--'
performanceSummary4$frequency_dependence_4[which(performanceSummary4$groupSize == 1)] = 'Single-condition'
performanceSummary4$frequency_dependence_4[which(performanceSummary4$groupSize > 1 & performanceSummary4$frequency_dependence != "--")] = 'Random-copying'
performanceSummary4$frequency_dependence_4[which(performanceSummary4$theta_75 > 0&performanceSummary4$theta_25 < 0)] = 'Random-copying'
performanceSummary4$frequency_dependence_4[which(performanceSummary4$theta_75 < 0)] = 'Neg-freq-dep'
performanceSummary4$frequency_dependence_4[which(performanceSummary4$theta_25 > 0)] = 'Pos-freq-dep'
#performanceSummary4$frequency_dependence_4[which(performanceSummary4$smallestWaicModel == "AL_ann"|performanceSummary4$smallestWaicModel == "AL")] = 'Asocial'
#performanceSummary4$frequency_dependence_4[which(performanceSummary4$smallestWaicModel == "Random")] = 'Random'
performanceSummary4$frequency_dependence_4 = factor(performanceSummary4$frequency_dependence_4, levels=c("--","Single-condition","Random","Asocial","Random-copying",'Neg-freq-dep', "Pos-freq-dep"))


## Read the questionnaire data
data_all <- read.csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/data_all.csv")
data_all_additional_indiv <- read.csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/results/data_all_additional_indiv.csv")

performanceSummary4$male_female_other = NA
performanceSummary4$age = NA
for(i in names(table(performanceSummary4$amazonID))){
    if(length(data_all$male_female_other[which(data_all$amazonID==i)])>0) {
        performanceSummary4$male_female_other[which(performanceSummary4$amazonID==i)] <- data_all$male_female_other[which(data_all$amazonID==i)]
    }
    if(length(data_all$age[which(data_all$amazonID==i)])>0) {
        performanceSummary4$age[which(performanceSummary4$amazonID==i)] <- data_all$age[which(data_all$amazonID==i)]
    }
    if(length(data_all_additional_indiv$male_female_other[which(data_all_additional_indiv$amazonID==i)])>0) {
        performanceSummary4$male_female_other[which(performanceSummary4$amazonID==i)] <- data_all_additional_indiv$male_female_other[which(data_all_additional_indiv$amazonID==i)]
    }
    if(length(data_all_additional_indiv$age[which(data_all_additional_indiv$amazonID==i)])>0) {
        performanceSummary4$age[which(performanceSummary4$amazonID==i)] <- data_all_additional_indiv$age[which(data_all_additional_indiv$amazonID==i)]
    }
}

performanceSummary4$average_copy_rate = NA
for(i in names(table(performanceSummary4$amazonID))){
    performanceSummary4$average_copy_rate[which(performanceSummary4$amazonID==i)] <- mean(subset(copyingProbReduction_data_fullModelOnly,amazonID==i)$copyingProb, na.rm=TRUE)
}

performanceSummary4$average_invTemp = NA
for(i in names(table(performanceSummary4$amazonID))){
    performanceSummary4$average_invTemp[which(performanceSummary4$amazonID==i)] <- mean(subset(copyingProbReduction_data_fullModelOnly,amazonID==i)$beta, na.rm=TRUE)
}

write.csv(performanceSummary4,
            "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/performanceSummary4.csv",
            row.names=FALSE)
