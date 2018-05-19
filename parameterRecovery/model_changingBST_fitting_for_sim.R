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
  	text = element_text(family="Times"),
    legend.position = 'right',
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Times" ), #"Times"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=15, family="Times" ,colour='black'),
    axis.text.y = element_text(size=15, family="Times" ,colour='black'),
    legend.text = element_text(size=15, family="Times" ,colour='black'),
    legend.title = element_text(size=16, family="Times" ,colour='black'),
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

################################################################################
##
##                         changingBST model
##
################################################################################
modelSimulation_changingBST_25 <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/modelSimulation_changingBST_25.csv", col_names = FALSE)
modelSimulation_changingBST_50 <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/modelSimulation_changingBST_50.csv", col_names = FALSE)
modelSimulation_changingBST_78 <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/modelSimulation_changingBST_78.csv", col_names = FALSE)

modelSimulation_changingBST = rbind(rbind(modelSimulation_changingBST_25,modelSimulation_changingBST_50),modelSimulation_changingBST_78)
colnames(modelSimulation_changingBST) <- c('Rep','groupSize','alpha','beta','annealing','netBeta','soc_raw','soc_reduc','soc','theta0','theta_slope','theta','Y','outcome','performance','socialFreq0','socialFreq1','socialFreq2')
modelSimulation_changingBST$taskDifficulty = c(rep('25%',nrow(modelSimulation_changingBST_25)),rep('50%',nrow(modelSimulation_changingBST_50)),rep('78%',nrow(modelSimulation_changingBST_78)))
modelSimulation_changingBST$round = rep(1:70, nrow(modelSimulation_changingBST)/70)
modelSimulation_changingBST$subject = NA
for(t in 1:70){
	modelSimulation_changingBST$subject[which(modelSimulation_changingBST$round==t)] = 1:(nrow(modelSimulation_changingBST)/70)
}
modelSimulation_changingBST$optimalChoice = modelSimulation_changingBST$performance
modelSimulation_changingBST$optimalChoice[which(modelSimulation_changingBST$optimalChoice==1&modelSimulation_changingBST$round>40)] <- 0
modelSimulation_changingBST$optimalChoice[which(modelSimulation_changingBST$optimalChoice==2)] <- 1
modelSimulation_changingBST$con = c(rep(1,nrow(modelSimulation_changingBST_25)),rep(2,nrow(modelSimulation_changingBST_50)),rep(3,nrow(modelSimulation_changingBST_78)))

ggplot(modelSimulation_changingBST,aes(outcome,fill=as.factor(Y)))+geom_histogram()+facet_grid(.~taskDifficulty)
#ggplot(modelSimulation_changingBST,aes(round,performance,colour=groupSize)) +
	stat_summary(aes(group=groupSize), fun.y = mean, geom="line") + ylim(c(0,2))+
	facet_grid(.~taskDifficulty)
ggplot(modelSimulation_changingBST,aes(round,optimalChoice,colour=groupSize)) +
	stat_summary(aes(group=groupSize), fun.y = mean, geom="line") + ylim(c(0,1))+
	facet_grid(.~taskDifficulty)

#ggplot(modelSimulation_changingBST)+geom_line(aes(round,Y))+geom_point(aes(round,Y),alpha=1/4)+facet_grid(groupSize+Rep~taskDifficulty)

ReinforcementLearningStanData_simulation = list(
	All = nrow(modelSimulation_changingBST),
	Nsub = length(table(modelSimulation_changingBST$subject)),
	Ncue = 3, # number of options
	Ntrial = 70,
	Ncon = 3,
	sub = modelSimulation_changingBST$subject,
	Y = modelSimulation_changingBST$Y,
	trial = modelSimulation_changingBST$round,
	outcome = modelSimulation_changingBST$outcome,
	trueAlpha=modelSimulation_changingBST$alpha,
	trueBeta = modelSimulation_changingBST$beta,
	trueAnnealing = modelSimulation_changingBST$annealing,
	trueNetBeta=modelSimulation_changingBST$netBeta,
	trueSoc=modelSimulation_changingBST$soc,
	trueSocRaw=modelSimulation_changingBST$soc_raw,
	trueSocReduc=modelSimulation_changingBST$soc_reduc,
	trueTheta0=modelSimulation_changingBST$theta0,
	trueThetaSlope=modelSimulation_changingBST$theta_slope,
	trueTheta=modelSimulation_changingBST$theta
	)
con = numeric(length(table(modelSimulation_changingBST$subject)))
for(i in 1:ReinforcementLearningStanData_simulation$Nsub){
	con[i] <- modelSimulation_changingBST$con[which(modelSimulation_changingBST$subject==i)][1]
}
ReinforcementLearningStanData_simulation$condition = con
startT = c()
MaxTrial = c()
for(i in 1:length(table(modelSimulation_changingBST$subject))){
	startT = append(startT, min(subset(modelSimulation_changingBST, subject==i)$round))
	MaxTrial = append(MaxTrial, length(which(subset(modelSimulation_changingBST, subject==i)$Y>=1)))
}
ReinforcementLearningStanData_simulation$startT = startT
ReinforcementLearningStanData_simulation$MaxTrial = MaxTrial

## compfun() で高速化
library(compiler)
F_calculation_simulation = function(array, Nsub, Ntrial, data)
{
	F <- array
	# F (Experimental data とは t+1 ずれてるから注意！)
	for(i in 1:Nsub) {
		F[i,1,1] = 0; F[i,2,1] = 0; F[i,3,1] = 0;
		for(t in 1:(Ntrial-1)) {
			lastChoice = 0
			if(subset(data,subject==i&round==t)$Y>0) {
				lastChoice = subset(data,subject==i&round==t)$Y
			}
			F[i,1,(t+1)] = subset(data,subject==i&round==t)$socialFreq0
			F[i,2,(t+1)] = subset(data,subject==i&round==t)$socialFreq1
			F[i,3,(t+1)] = subset(data,subject==i&round==t)$socialFreq2
			if(lastChoice>0){
				F[i,lastChoice,(t+1)] = F[i,lastChoice,(t+1)] - 1
			}
			if(length(which(F[i,,(t+1)]>=0))==0){
				F[i,,t] <- c(-1,-1,-1)
			}
		}
	}
	return(F)
}

F_calculation_simulation.compiled <- cmpfun(F_calculation_simulation)
F0_sim = array(rep(NA,nrow(modelSimulation_changingBST)), c(ReinforcementLearningStanData_simulation$Nsub, ReinforcementLearningStanData_simulation$Ncue, ReinforcementLearningStanData_simulation$Ntrial))
F = F_calculation_simulation.compiled(F0_sim,ReinforcementLearningStanData_simulation$Nsub,ReinforcementLearningStanData_simulation$Ntrial,modelSimulation_changingBST)
ReinforcementLearningStanData_simulation$F = round(F)


################################################################################
##
##                FITTING: changingBST model -- STANDARD
##
################################################################################

	## FITTING
model_changingBST = stan_model(file="~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model_changingBST.stan")
fit_changingBST_simulation = sampling(
	model_changingBST, data=ReinforcementLearningStanData_simulation, seed=77,
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
    chains=8, iter=2000, warmup=1000, thin=5
    #chains=8, iter=1200, warmup=600, thin=3
)
fit_changingBST_simulation
saveRDS(fit_changingBST_simulation, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/fit_changingBST_simulation.rds')
ms_changingBST_simulation <- rstan::extract(fit_changingBST_simulation)
plot(fit_changingBST_simulation, pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_soc','s_soc','mu_theta','s_theta','mu_theta_slope','s_theta_slope','mu_soc_reduction','s_soc_reduction','mu_annealing','s_annealing'))

	## 事後診断
traceplot(fit_changingBST_simulation, pars=c('mu_alpha','mu_beta','mu_annealing'), inc_warmup=TRUE)
traceplot(fit_changingBST_simulation, pars=c('mu_theta','mu_theta_slope','mu_soc','mu_soc_reduction'), inc_warmup=TRUE)
traceplot(fit_changingBST_simulation, pars=c('s_alpha','s_beta','s_theta','s_theta_slope','s_soc','s_soc_reduction'), inc_warmup=TRUE)

temp_changingBST_simulation = data.frame(WAICi_changingBST = WAIC_indv(fit_changingBST_simulation)$waic)

temp_changingBST_simulation <- cbind(temp_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$alpha[,1:ncol(ms_changingBST_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST_simulation <- cbind(temp_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$beta[,1:ncol(ms_changingBST_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST_simulation <- cbind(temp_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$annealing[,1:ncol(ms_changingBST_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST_simulation <- cbind(temp_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$theta[,1:ncol(ms_changingBST_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST_simulation <- cbind(temp_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$theta_slope[,1:ncol(ms_changingBST_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST_simulation <- cbind(temp_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$soc_reduction[,1:ncol(ms_changingBST_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_changingBST_simulation <- cbind(temp_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$soc_raw[,1:ncol(ms_changingBST_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(temp_changingBST_simulation) <- c("changingBST_simulation_WAICi",
							"changingBST_simulation_alpha_p2.5", "changingBST_simulation_alpha_p25", "changingBST_simulation_alpha_p50", "changingBST_simulation_alpha_p75", "changingBST_simulation_alpha_p97.5",
							"changingBST_simulation_beta_p2.5", "changingBST_simulation_beta_p25", "changingBST_simulation_beta_p50", "changingBST_simulation_beta_p75", "changingBST_simulation_beta_p97.5",
							"changingBST_simulation_annealing_p2.5", "changingBST_simulation_annealing_p25", "changingBST_simulation_annealing_p50", "changingBST_simulation_annealing_p75", "changingBST_simulation_annealing_p97.5",
							"changingBST_simulation_theta_p2.5", "changingBST_simulation_theta_p25", "changingBST_simulation_theta_p50", "changingBST_simulation_theta_p75", "changingBST_simulation_theta_p97.5",
							"changingBST_simulation_theta_slope_p2.5", "changingBST_simulation_theta_slope_p25", "changingBST_simulation_theta_slope_p50", "changingBST_simulation_theta_slope_p75", "changingBST_simulation_theta_slope_p97.5",
							"changingBST_simulation_soc_reduction_p2.5", "changingBST_simulation_soc_reduction_p25", "changingBST_simulation_soc_reduction_p50", "changingBST_simulation_soc_reduction_p75", "changingBST_simulation_soc_reduction_p97.5",
							"changingBST_simulation_soc_raw_p2.5", "changingBST_simulation_soc_raw_p25", "changingBST_simulation_soc_raw_p50", "changingBST_simulation_soc_raw_p75", "changingBST_simulation_soc_raw_p97.5")
max_log_lik = apply(ms_changingBST_simulation$log_lik[,1:ncol(ms_changingBST_simulation$log_lik)], 2, which.max)
max_changingBST_simulation_alpha=rep(NA, ncol(ms_changingBST_simulation$log_lik))
max_changingBST_simulation_beta=rep(NA, ncol(ms_changingBST_simulation$log_lik))
max_changingBST_simulation_annealing=rep(NA, ncol(ms_changingBST_simulation$log_lik))
max_changingBST_simulation_theta=rep(NA, ncol(ms_changingBST_simulation$log_lik))
max_changingBST_simulation_theta_slope=rep(NA, ncol(ms_changingBST_simulation$log_lik))
max_changingBST_simulation_soc_reduction=rep(NA, ncol(ms_changingBST_simulation$log_lik))
max_changingBST_simulation_soc_raw=rep(NA, ncol(ms_changingBST_simulation$log_lik))
for(i in 1:ncol(ms_changingBST_simulation$log_lik)){
	max_changingBST_simulation_alpha[i] <- ms_changingBST_simulation$alpha[max_log_lik[i],i]
	max_changingBST_simulation_beta[i] <- ms_changingBST_simulation$beta[max_log_lik[i],i]
	max_changingBST_simulation_annealing[i] <- ms_changingBST_simulation$annealing[max_log_lik[i],i]
	max_changingBST_simulation_theta[i] <- ms_changingBST_simulation$theta[max_log_lik[i],i]
	max_changingBST_simulation_theta_slope[i] <- ms_changingBST_simulation$theta_slope[max_log_lik[i],i]
	max_changingBST_simulation_soc_reduction[i] <- ms_changingBST_simulation$soc_reduction[max_log_lik[i],i]
	max_changingBST_simulation_soc_raw[i] <- ms_changingBST_simulation$soc_raw[max_log_lik[i],i]
}
temp_changingBST_simulation$max_changingBST_simulation_alpha = max_changingBST_simulation_alpha
temp_changingBST_simulation$max_changingBST_simulation_beta = max_changingBST_simulation_beta
temp_changingBST_simulation$max_changingBST_simulation_annealing = max_changingBST_simulation_annealing
temp_changingBST_simulation$max_changingBST_simulation_theta = max_changingBST_simulation_theta
temp_changingBST_simulation$max_changingBST_simulation_theta_slope = max_changingBST_simulation_theta_slope
temp_changingBST_simulation$max_changingBST_simulation_soc_reduction = max_changingBST_simulation_soc_reduction
temp_changingBST_simulation$max_changingBST_simulation_soc_raw = max_changingBST_simulation_soc_raw
temp_changingBST_simulation$max_changingBST_simulation_log_lik = apply(ms_changingBST_simulation$log_lik[,1:ncol(ms_changingBST_simulation$log_lik)], 2, max)
numParam = 7
temp_changingBST_simulation$AIC_changingBST = -temp_changingBST_simulation$max_changingBST_simulation_log_lik + numParam ## 正確には AIC/2
temp_changingBST_simulation$BIC_changingBST = -temp_changingBST_simulation$max_changingBST_simulation_log_lik + (numParam/2)*log(70) ## more precisely, BIC/2
temp_changingBST_simulation$WAICi_changingBST = WAIC_indv(fit_changingBST_simulation)$waic


temp_changingBST_simulation$log_lik_n_eff = NA
temp_changingBST_simulation$log_lik_Rhat = NA
temp_changingBST_simulation$alpha_n_eff = NA
temp_changingBST_simulation$alpha_Rhat = NA
temp_changingBST_simulation$beta_n_eff = NA
temp_changingBST_simulation$beta_Rhat = NA
temp_changingBST_simulation$soc_raw_n_eff = NA
temp_changingBST_simulation$soc_raw_Rhat = NA
temp_changingBST_simulation$soc_reduction_n_eff = NA
temp_changingBST_simulation$soc_reduction_Rhat = NA
temp_changingBST_simulation$theta_n_eff = NA
temp_changingBST_simulation$theta_Rhat = NA
temp_changingBST_simulation$theta_slope_n_eff = NA
temp_changingBST_simulation$theta_slope_Rhat = NA
temp_changingBST_simulation$annealing_n_eff = NA
temp_changingBST_simulation$annealing_Rhat = NA


	## Estimated copying rate at each round
soc_table_changingBST_simulation = data.frame(subject = 1:ncol(ms_changingBST_simulation$log_lik), round = rep(1, ncol(ms_changingBST_simulation$log_lik)))
soc_table_changingBST_simulation = cbind(soc_table_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$soc[,1:ncol(ms_changingBST_simulation$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(soc_table_changingBST_simulation) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
	soc_table_changingBST_simulation = rbind(soc_table_changingBST_simulation, cbind(cbind(1:ncol(ms_changingBST_simulation$log_lik), rep(t,ncol(ms_changingBST_simulation$log_lik))), as.data.frame(t(apply(ms_changingBST_simulation$soc[,1:ncol(ms_changingBST_simulation$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(soc_table_changingBST_simulation) <- c("subject", "round", "changingBST_soc_p2.5", "changingBST_soc_p25", "changingBST_soc_p50", "changingBST_soc_p75", "changingBST_soc_p97.5")

write.csv(soc_table_changingBST_simulation,
			"~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/soc_table_changingBST_simulation.csv",
			row.names=FALSE)

	## Estimated netBeta at each round
netBeta_table_changingBST_simulation = data.frame(subject = 1:ncol(ms_changingBST_simulation$log_lik), round = rep(1, ncol(ms_changingBST_simulation$log_lik)))
netBeta_table_changingBST_simulation = cbind(netBeta_table_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$netBeta[,1:ncol(ms_changingBST_simulation$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(netBeta_table_changingBST_simulation) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
	netBeta_table_changingBST_simulation = rbind(netBeta_table_changingBST_simulation, cbind(cbind(1:ncol(ms_changingBST_simulation$log_lik), rep(t,ncol(ms_changingBST_simulation$log_lik))), as.data.frame(t(apply(ms_changingBST_simulation$netBeta[,1:ncol(ms_changingBST_simulation$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(netBeta_table_changingBST_simulation) <- c("subject", "round", "changingBST_netBeta_p2.5", "changingBST_netBeta_p25", "changingBST_netBeta_p50", "changingBST_netBeta_p75", "changingBST_netBeta_p97.5")

write.csv(netBeta_table_changingBST_simulation,
			"~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/netBeta_table_changingBST_simulation.csv", row.names=FALSE)

	## Estimated netTheta at each round
netTheta_table_changingBST_simulation = data.frame(subject = 1:ncol(ms_changingBST_simulation$log_lik), round = rep(1, ncol(ms_changingBST_simulation$log_lik)))
netTheta_table_changingBST_simulation = cbind(netTheta_table_changingBST_simulation, as.data.frame(t(apply(ms_changingBST_simulation$netTheta[,1:ncol(ms_changingBST_simulation$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(netTheta_table_changingBST_simulation) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
	netTheta_table_changingBST_simulation = rbind(netTheta_table_changingBST_simulation, cbind(cbind(1:ncol(ms_changingBST_simulation$log_lik), rep(t,ncol(ms_changingBST_simulation$log_lik))), as.data.frame(t(apply(ms_changingBST_simulation$netTheta[,1:ncol(ms_changingBST_simulation$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(netTheta_table_changingBST_simulation) <- c("subject", "round", "changingBST_netTheta_p2.5", "changingBST_netTheta_p25", "changingBST_netTheta_p50", "changingBST_netTheta_p75", "changingBST_netTheta_p97.5")

write.csv(netTheta_table_changingBST_simulation,
			"~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/netTheta_table_changingBST_simulation.csv", row.names=FALSE)


## Comparing with the TRUE values
## Were parameter values correctly estimated?
modelSimulation_changingBST$estimated_alpha = NA
modelSimulation_changingBST$estimated_alpha_high = NA
modelSimulation_changingBST$estimated_alpha_low = NA
modelSimulation_changingBST$estimated_alpha_25 = NA
modelSimulation_changingBST$estimated_alpha_75 = NA

modelSimulation_changingBST$estimated_beta = NA
modelSimulation_changingBST$estimated_beta_high = NA
modelSimulation_changingBST$estimated_beta_low = NA
modelSimulation_changingBST$estimated_beta_25 = NA
modelSimulation_changingBST$estimated_beta_75 = NA

modelSimulation_changingBST$estimated_annealing = NA
modelSimulation_changingBST$estimated_annealing_high = NA
modelSimulation_changingBST$estimated_annealing_low = NA
modelSimulation_changingBST$estimated_annealing_25 = NA
modelSimulation_changingBST$estimated_annealing_75 = NA

modelSimulation_changingBST$estimated_netBeta = NA
modelSimulation_changingBST$estimated_netBeta_high = NA
modelSimulation_changingBST$estimated_netBeta_low = NA
modelSimulation_changingBST$estimated_netBeta_25 = NA
modelSimulation_changingBST$estimated_netBeta_75 = NA

modelSimulation_changingBST$estimated_theta = NA
modelSimulation_changingBST$estimated_theta_high = NA
modelSimulation_changingBST$estimated_theta_low = NA
modelSimulation_changingBST$estimated_theta_25 = NA
modelSimulation_changingBST$estimated_theta_75 = NA

modelSimulation_changingBST$estimated_theta_slope = NA
modelSimulation_changingBST$estimated_theta_slope_high = NA
modelSimulation_changingBST$estimated_theta_slope_low = NA
modelSimulation_changingBST$estimated_theta_slope_25 = NA
modelSimulation_changingBST$estimated_theta_slope_75 = NA

modelSimulation_changingBST$estimated_netTheta = NA
modelSimulation_changingBST$estimated_netTheta_high = NA
modelSimulation_changingBST$estimated_netTheta_low = NA
modelSimulation_changingBST$estimated_netTheta_25 = NA
modelSimulation_changingBST$estimated_netTheta_75 = NA

modelSimulation_changingBST$estimated_soc_raw = NA
modelSimulation_changingBST$estimated_soc_raw_high = NA
modelSimulation_changingBST$estimated_soc_raw_low = NA
modelSimulation_changingBST$estimated_soc_raw_25 = NA
modelSimulation_changingBST$estimated_soc_raw_75 = NA

modelSimulation_changingBST$estimated_soc_reduction = NA
modelSimulation_changingBST$estimated_soc_reduction_high = NA
modelSimulation_changingBST$estimated_soc_reduction_low = NA
modelSimulation_changingBST$estimated_soc_reduction_25 = NA
modelSimulation_changingBST$estimated_soc_reduction_75 = NA

modelSimulation_changingBST$estimated_soc = NA
modelSimulation_changingBST$estimated_soc_high = NA
modelSimulation_changingBST$estimated_soc_low = NA
modelSimulation_changingBST$estimated_soc_25 = NA
modelSimulation_changingBST$estimated_soc_75 = NA

modelSimulation_changingBST$bestfit_alpha = NA
modelSimulation_changingBST$bestfit_beta = NA
modelSimulation_changingBST$bestfit_annealing = NA
modelSimulation_changingBST$bestfit_netBeta = NA
modelSimulation_changingBST$bestfit_theta = NA
modelSimulation_changingBST$bestfit_theta_slope = NA
modelSimulation_changingBST$bestfit_netTheta = NA
modelSimulation_changingBST$bestfit_soc_raw = NA
modelSimulation_changingBST$bestfit_soc_reduction = NA
modelSimulation_changingBST$bestfit_soc = NA

for (i in 1:nrow(temp_changingBST_simulation)) {
        # alpha
    modelSimulation_changingBST$bestfit_alpha[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$max_changingBST_simulation_alpha[i]
    modelSimulation_changingBST$estimated_alpha[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_alpha_p50[i]
    modelSimulation_changingBST$estimated_alpha_low[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_alpha_p2.5[i]
    modelSimulation_changingBST$estimated_alpha_high[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_alpha_p97.5[i]
    modelSimulation_changingBST$estimated_alpha_25[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_alpha_p25[i]
    modelSimulation_changingBST$estimated_alpha_75[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_alpha_p75[i]
        # beta
    modelSimulation_changingBST$estimated_beta[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_beta_p50[i]
    modelSimulation_changingBST$estimated_beta_low[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_beta_p2.5[i]
    modelSimulation_changingBST$estimated_beta_high[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_beta_p97.5[i]
    modelSimulation_changingBST$bestfit_beta[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$max_changingBST_simulation_beta[i]
    modelSimulation_changingBST$estimated_beta_25[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_beta_p25[i]
    modelSimulation_changingBST$estimated_beta_75[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_beta_p75[i]
        # annealing
    modelSimulation_changingBST$estimated_annealing[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_annealing_p50[i]
    modelSimulation_changingBST$estimated_annealing_low[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_annealing_p2.5[i]
    modelSimulation_changingBST$estimated_annealing_high[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_annealing_p97.5[i]
    modelSimulation_changingBST$bestfit_annealing[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$max_changingBST_simulation_annealing[i]
    modelSimulation_changingBST$estimated_annealing_25[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_annealing_p25[i]
    modelSimulation_changingBST$estimated_annealing_75[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_annealing_p75[i]
        # theta
    modelSimulation_changingBST$estimated_theta[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_p50[i]
    modelSimulation_changingBST$estimated_theta_low[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_p2.5[i]
    modelSimulation_changingBST$estimated_theta_high[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_p97.5[i]
    modelSimulation_changingBST$bestfit_theta[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$max_changingBST_simulation_theta[i]
    modelSimulation_changingBST$estimated_theta_25[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_p25[i]
    modelSimulation_changingBST$estimated_theta_75[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_p75[i]
        # theta_slope
    modelSimulation_changingBST$estimated_theta_slope[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_slope_p50[i]
    modelSimulation_changingBST$estimated_theta_slope_low[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_slope_p2.5[i]
    modelSimulation_changingBST$estimated_theta_slope_high[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_slope_p97.5[i]
    modelSimulation_changingBST$bestfit_theta_slope[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$max_changingBST_simulation_theta_slope[i]
    modelSimulation_changingBST$estimated_theta_slope_25[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_slope_p25[i]
    modelSimulation_changingBST$estimated_theta_slope_75[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_theta_slope_p75[i]
        # soc_raw
    modelSimulation_changingBST$estimated_soc_raw[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_raw_p50[i]
    modelSimulation_changingBST$estimated_soc_raw_low[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_raw_p2.5[i]
    modelSimulation_changingBST$estimated_soc_raw_high[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_raw_p97.5[i]
    modelSimulation_changingBST$bestfit_soc_raw[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$max_changingBST_simulation_soc_raw[i]
    modelSimulation_changingBST$estimated_soc_raw_25[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_raw_p25[i]
    modelSimulation_changingBST$estimated_soc_raw_75[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_raw_p75[i]
        # soc_reduction
    modelSimulation_changingBST$estimated_soc_reduction[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_reduction_p50[i]
    modelSimulation_changingBST$estimated_soc_reduction_low[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_reduction_p2.5[i]
    modelSimulation_changingBST$estimated_soc_reduction_high[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_reduction_p97.5[i]
    modelSimulation_changingBST$bestfit_soc_reduction[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$max_changingBST_simulation_soc_reduction[i]
    modelSimulation_changingBST$estimated_soc_reduction_25[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_reduction_p25[i]
    modelSimulation_changingBST$estimated_soc_reduction_75[which(modelSimulation_changingBST$subject==i)] <- temp_changingBST_simulation$changingBST_simulation_soc_reduction_p75[i]
}

for (i in 1:nrow(temp_changingBST_simulation)) {
	for (t in 1:70) {
		# soc
		modelSimulation_changingBST$estimated_soc[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =soc_table_changingBST_simulation$changingBST_soc_p50[which(soc_table_changingBST_simulation$subject==i&soc_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_soc_low[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =soc_table_changingBST_simulation$changingBST_soc_p2.5[which(soc_table_changingBST_simulation$subject==i&soc_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_soc_high[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =soc_table_changingBST_simulation$changingBST_soc_p97.5[which(soc_table_changingBST_simulation$subject==i&soc_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_soc_25[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =soc_table_changingBST_simulation$changingBST_soc_p25[which(soc_table_changingBST_simulation$subject==i&soc_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_soc_75[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =soc_table_changingBST_simulation$changingBST_soc_p75[which(soc_table_changingBST_simulation$subject==i&soc_table_changingBST_simulation$round==t)]
		# netBeta
		modelSimulation_changingBST$estimated_netBeta[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netBeta_table_changingBST_simulation$changingBST_netBeta_p50[which(netBeta_table_changingBST_simulation$subject==i&netBeta_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_netBeta_low[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netBeta_table_changingBST_simulation$changingBST_netBeta_p2.5[which(netBeta_table_changingBST_simulation$subject==i&netBeta_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_netBeta_high[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netBeta_table_changingBST_simulation$changingBST_netBeta_p97.5[which(netBeta_table_changingBST_simulation$subject==i&netBeta_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_netBeta_25[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netBeta_table_changingBST_simulation$changingBST_netBeta_p25[which(netBeta_table_changingBST_simulation$subject==i&netBeta_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_netBeta_75[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netBeta_table_changingBST_simulation$changingBST_netBeta_p75[which(netBeta_table_changingBST_simulation$subject==i&netBeta_table_changingBST_simulation$round==t)]
		# netTheta
		modelSimulation_changingBST$estimated_netTheta[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netTheta_table_changingBST_simulation$changingBST_netTheta_p50[which(netTheta_table_changingBST_simulation$subject==i&netTheta_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_netTheta_low[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netTheta_table_changingBST_simulation$changingBST_netTheta_p2.5[which(netTheta_table_changingBST_simulation$subject==i&netTheta_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_netTheta_high[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netTheta_table_changingBST_simulation$changingBST_netTheta_p97.5[which(netTheta_table_changingBST_simulation$subject==i&netTheta_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_netTheta_25[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netTheta_table_changingBST_simulation$changingBST_netTheta_p25[which(netTheta_table_changingBST_simulation$subject==i&netTheta_table_changingBST_simulation$round==t)]
		modelSimulation_changingBST$estimated_netTheta_75[which(modelSimulation_changingBST$subject==i&modelSimulation_changingBST$round==t)] =netTheta_table_changingBST_simulation$changingBST_netTheta_p75[which(netTheta_table_changingBST_simulation$subject==i&netTheta_table_changingBST_simulation$round==t)]
	}
}











###################################################################
##
## Parameter recovery test
##
###################################################################

## Alpha
alpha_recovery_data = data.frame(
								true_alpha = modelSimulation_changingBST$alpha[which(modelSimulation_changingBST$round==1)],
								estimated_alpha = modelSimulation_changingBST$estimated_alpha[which(modelSimulation_changingBST$round==1)],
								bestfit_alpha = modelSimulation_changingBST$bestfit_alpha[which(modelSimulation_changingBST$round==1)],
								estimated_alpha_low = modelSimulation_changingBST$estimated_alpha_low[which(modelSimulation_changingBST$round==1)],
								taskDifficulty = modelSimulation_changingBST$taskDifficulty[which(modelSimulation_changingBST$round==1)],
								estimated_alpha_25 = modelSimulation_changingBST$estimated_alpha_25[which(modelSimulation_changingBST$round==1)],
								estimated_alpha_75 = modelSimulation_changingBST$estimated_alpha_75[which(modelSimulation_changingBST$round==1)],
								estimated_alpha_high = modelSimulation_changingBST$estimated_alpha_high[which(modelSimulation_changingBST$round==1)]
								)
subject_order_alpha = order(modelSimulation_changingBST$alpha[which(modelSimulation_changingBST$round==1)])
alpha_recovery_data = alpha_recovery_data[subject_order_alpha,]
alpha_recovery_data$subject = 1:nrow(alpha_recovery_data)
alpha_plot = ggplot(data=alpha_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_alpha_low, ymax=estimated_alpha_high), width=.001) +
	geom_point(aes(subject, estimated_alpha), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_alpha), colour='orange',size=.5)+
	geom_point(aes(subject, true_alpha), colour='red',size=.5) +
	myTheme()+
	ylim(c(0,1)) +
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ',alpha,sep="")))

## beta
beta_recovery_data = data.frame(
								true_beta = modelSimulation_changingBST$beta[which(modelSimulation_changingBST$round==1)],
								estimated_beta = modelSimulation_changingBST$estimated_beta[which(modelSimulation_changingBST$round==1)],
								bestfit_beta = modelSimulation_changingBST$bestfit_beta[which(modelSimulation_changingBST$round==1)],
								estimated_beta_low = modelSimulation_changingBST$estimated_beta_low[which(modelSimulation_changingBST$round==1)],
								taskDifficulty = modelSimulation_changingBST$taskDifficulty[which(modelSimulation_changingBST$round==1)],
								estimated_beta_25 = modelSimulation_changingBST$estimated_beta_25[which(modelSimulation_changingBST$round==1)],
								estimated_beta_75 = modelSimulation_changingBST$estimated_beta_75[which(modelSimulation_changingBST$round==1)],
								estimated_beta_high = modelSimulation_changingBST$estimated_beta_high[which(modelSimulation_changingBST$round==1)]
								)
subject_order_beta = order(modelSimulation_changingBST$beta[which(modelSimulation_changingBST$round==1)])
beta_recovery_data = beta_recovery_data[subject_order_beta,]
beta_recovery_data$subject = 1:nrow(beta_recovery_data)
beta_plot = ggplot(data=beta_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_beta_low, ymax=estimated_beta_high), width=.001) +
	geom_point(aes(subject, estimated_beta), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_beta), colour='orange',size=.5)+
	geom_point(aes(subject, true_beta), colour='red',size=.5)+
	myTheme()+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ',beta[0],sep="")))

## annealing (epsilon)
annealing_recovery_data = data.frame(
								true_annealing = modelSimulation_changingBST$annealing[which(modelSimulation_changingBST$round==1)],
								estimated_annealing = modelSimulation_changingBST$estimated_annealing[which(modelSimulation_changingBST$round==1)],
								bestfit_annealing = modelSimulation_changingBST$bestfit_annealing[which(modelSimulation_changingBST$round==1)],
								estimated_annealing_low = modelSimulation_changingBST$estimated_annealing_low[which(modelSimulation_changingBST$round==1)],
								taskDifficulty = modelSimulation_changingBST$taskDifficulty[which(modelSimulation_changingBST$round==1)],
								estimated_annealing_25 = modelSimulation_changingBST$estimated_annealing_25[which(modelSimulation_changingBST$round==1)],
								estimated_annealing_75 = modelSimulation_changingBST$estimated_annealing_75[which(modelSimulation_changingBST$round==1)],
								estimated_annealing_high = modelSimulation_changingBST$estimated_annealing_high[which(modelSimulation_changingBST$round==1)]
								)
subject_order_annealing = order(modelSimulation_changingBST$annealing[which(modelSimulation_changingBST$round==1)])
annealing_recovery_data = annealing_recovery_data[subject_order_annealing,]
annealing_recovery_data$subject = 1:nrow(annealing_recovery_data)
annealing_plot = ggplot(data=annealing_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_annealing_low, ymax=estimated_annealing_high), width=.001) +
	geom_point(aes(subject, estimated_annealing), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_annealing), colour='orange',size=.5)+
	geom_point(aes(subject, true_annealing), colour='red',size=.5)+
	myTheme()+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ', epsilon, sep="")))

## theta_slope (gamma)
theta_slope_recovery_data = data.frame(
								true_theta_slope = modelSimulation_changingBST$theta_slope[which(modelSimulation_changingBST$round==1)],
								estimated_theta_slope = modelSimulation_changingBST$estimated_theta_slope[which(modelSimulation_changingBST$round==1)],
								bestfit_theta_slope = modelSimulation_changingBST$bestfit_theta_slope[which(modelSimulation_changingBST$round==1)],
								estimated_theta_slope_low = modelSimulation_changingBST$estimated_theta_slope_low[which(modelSimulation_changingBST$round==1)],
								taskDifficulty = modelSimulation_changingBST$taskDifficulty[which(modelSimulation_changingBST$round==1)],
								estimated_theta_slope_25 = modelSimulation_changingBST$estimated_theta_slope_25[which(modelSimulation_changingBST$round==1)],
								estimated_theta_slope_75 = modelSimulation_changingBST$estimated_theta_slope_75[which(modelSimulation_changingBST$round==1)],
								estimated_theta_slope_high = modelSimulation_changingBST$estimated_theta_slope_high[which(modelSimulation_changingBST$round==1)]
								)
subject_order_theta_slope = order(modelSimulation_changingBST$theta_slope[which(modelSimulation_changingBST$round==1)])
theta_slope_recovery_data = theta_slope_recovery_data[subject_order_theta_slope,]
theta_slope_recovery_data$subject = 1:nrow(theta_slope_recovery_data)
theta_slope_plot = ggplot(data=theta_slope_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_theta_slope_low, ymax=estimated_theta_slope_high), width=.001) +
	geom_point(aes(subject, estimated_theta_slope), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_theta_slope), colour='orange',size=.5)+
	geom_point(aes(subject, true_theta_slope), colour='red',size=.5)+
	myTheme()+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ', gamma, sep="")))

## theta
theta_recovery_data = data.frame(
								true_theta = modelSimulation_changingBST$theta[which(modelSimulation_changingBST$round==1)],
								estimated_theta = modelSimulation_changingBST$estimated_theta[which(modelSimulation_changingBST$round==1)],
								bestfit_theta = modelSimulation_changingBST$bestfit_theta[which(modelSimulation_changingBST$round==1)],
								estimated_theta_low = modelSimulation_changingBST$estimated_theta_low[which(modelSimulation_changingBST$round==1)],
								taskDifficulty = modelSimulation_changingBST$taskDifficulty[which(modelSimulation_changingBST$round==1)],
								estimated_theta_25 = modelSimulation_changingBST$estimated_theta_25[which(modelSimulation_changingBST$round==1)],
								estimated_theta_75 = modelSimulation_changingBST$estimated_theta_75[which(modelSimulation_changingBST$round==1)],
								estimated_theta_high = modelSimulation_changingBST$estimated_theta_high[which(modelSimulation_changingBST$round==1)]
								)

theta_recovery_data$true_mean_netTheta = tapply(modelSimulation_changingBST$theta, modelSimulation_changingBST$subject, mean)

theta_recovery_data$estimated_frequency_dependence = 'random-copying'
theta_recovery_data$estimated_frequency_dependence[which(theta_recovery_data$estimated_theta_high<0)] = 'neg-freq-dep'
theta_recovery_data$estimated_frequency_dependence[which(theta_recovery_data$estimated_theta_low>=0&theta_recovery_data$estimated_theta_low<=1)] = 'weak-pos-freq-dep'
theta_recovery_data$estimated_frequency_dependence[which(theta_recovery_data$estimated_theta_low>1)] = 'conformist'
theta_recovery_data$estimated_frequency_dependence2 = 'random-copying'
theta_recovery_data$estimated_frequency_dependence2[which(theta_recovery_data$estimated_theta_75<0)] = 'neg-freq-dep'
theta_recovery_data$estimated_frequency_dependence2[which(theta_recovery_data$estimated_theta_25>=0&theta_recovery_data$estimated_theta_25<=1)] = 'weak-pos-freq-dep'
theta_recovery_data$estimated_frequency_dependence2[which(theta_recovery_data$estimated_theta_25>1)] = 'conformist'
	# rand, neg, pos
theta_recovery_data$estimated_frequency_dependence31 = 'random-copying'
theta_recovery_data$estimated_frequency_dependence31[which(theta_recovery_data$estimated_theta_high<0)] = 'neg-freq-dep'
theta_recovery_data$estimated_frequency_dependence31[which(theta_recovery_data$estimated_theta_low>0)] = 'pos-freq-dep'
theta_recovery_data$estimated_frequency_dependence32 = 'random-copying'
theta_recovery_data$estimated_frequency_dependence32[which(theta_recovery_data$estimated_theta_75<0)] = 'neg-freq-dep'
theta_recovery_data$estimated_frequency_dependence32[which(theta_recovery_data$estimated_theta_25>0)] = 'pos-freq-dep'


subject_order_theta = order(modelSimulation_changingBST$theta[which(modelSimulation_changingBST$round==1)])
theta_recovery_data = theta_recovery_data[subject_order_theta,]
theta_recovery_data$subject = 1:nrow(theta_recovery_data)

theta_plot = ggplot(data=theta_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_theta_low, ymax=estimated_theta_high), width=.001, alpha=1/1) +
	geom_point(aes(subject, estimated_theta), colour='blue',size=.5)+
	geom_point(aes(subject, true_theta), colour='red',size=.5) +
	myTheme()+
	ylim(c(-33, 30))+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ',theta,sep="")))

theta_plot_95 = ggplot(data=theta_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_theta_low, ymax=estimated_theta_high, colour=estimated_frequency_dependence31), width=.001, alpha=1/1) +
	geom_point(aes(subject, estimated_theta, colour=estimated_frequency_dependence31),size=.5)+
	#geom_point(aes(subject, bestfit_theta), colour='orange',size=.5)+
	geom_point(aes(subject, true_theta),size=.5) +
	geom_hline(yintercept=0, colour='black') +
	#geom_hline(yintercept=1, colour='red') +
	scale_colour_manual(values=c("random-copying"="#999999",'neg-freq-dep'="blue", "weak-pos-freq-dep"="orange", "conformist"="red", "RANDOM"="#999999", "AL"="#999999", "pos-freq-dep"="red"), name="Strength of \nfrequency dependence") +
	myTheme()+
	ylim(c(-33, 30))+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ',theta,sep="")))
theta_plot_50 = ggplot(data=theta_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_theta_25, ymax=estimated_theta_75, colour=estimated_frequency_dependence32), width=.001, alpha=1/1) +
	geom_point(aes(subject, estimated_theta, colour=estimated_frequency_dependence32),size=.5)+
	#geom_point(aes(subject, bestfit_theta), colour='orange',size=.5)+
	geom_point(aes(subject, true_theta),size=.5) +
	geom_hline(yintercept=0, colour='black') +
	#geom_hline(yintercept=1, colour='red') +
	myTheme()+
	ylim(c(-33, 30))+
	scale_colour_manual(values=c("random-copying"="#999999",'neg-freq-dep'="blue", "weak-pos-freq-dep"="red", "conformist"="red", "RANDOM"="#999999", "AL"="#999999", "pos-freq-dep"="red"), name="Strength of \nfrequency dependence") +
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste(theta,' (using 50% CI)', sep="")))


## soc_raw (sigma)
soc_raw_recovery_data = data.frame(
								true_soc_raw = modelSimulation_changingBST$soc_raw[which(modelSimulation_changingBST$round==1)],
								estimated_soc_raw = modelSimulation_changingBST$estimated_soc_raw[which(modelSimulation_changingBST$round==1)],
								bestfit_soc_raw = modelSimulation_changingBST$bestfit_soc_raw[which(modelSimulation_changingBST$round==1)],
								taskDifficulty = modelSimulation_changingBST$taskDifficulty[which(modelSimulation_changingBST$round==1)],
								estimated_soc_raw_25 = modelSimulation_changingBST$estimated_soc_raw_25[which(modelSimulation_changingBST$round==1)],
								estimated_soc_raw_75 = modelSimulation_changingBST$estimated_soc_raw_75[which(modelSimulation_changingBST$round==1)],
								estimated_soc_raw_low = modelSimulation_changingBST$estimated_soc_raw_low[which(modelSimulation_changingBST$round==1)],
								estimated_soc_raw_high = modelSimulation_changingBST$estimated_soc_raw_high[which(modelSimulation_changingBST$round==1)]
								)
subject_order_soc_raw = order(modelSimulation_changingBST$soc_raw[which(modelSimulation_changingBST$round==1)])
soc_raw_recovery_data = soc_raw_recovery_data[subject_order_soc_raw,]
soc_raw_recovery_data$subject = 1:nrow(soc_raw_recovery_data)
soc_raw_plot =ggplot(data=soc_raw_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_soc_raw_low, ymax=estimated_soc_raw_high), width=.001) +
	geom_point(aes(subject, estimated_soc_raw), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_soc_raw), colour='orange',size=.5)+
	geom_point(aes(subject, true_soc_raw), colour='red',size=.5) +
	myTheme()+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ',sigma[0],sep="")))

## soc_reduction
modelSimulation_changingBST$soc_reduction = modelSimulation_changingBST$soc_reduc
soc_reduction_recovery_data = data.frame(
								true_soc_reduction = modelSimulation_changingBST$soc_reduction[which(modelSimulation_changingBST$round==1)],
								estimated_soc_reduction = modelSimulation_changingBST$estimated_soc_reduction[which(modelSimulation_changingBST$round==1)],
								bestfit_soc_reduction = modelSimulation_changingBST$bestfit_soc_reduction[which(modelSimulation_changingBST$round==1)],
								estimated_soc_reduction_25 = modelSimulation_changingBST$estimated_soc_reduction_25[which(modelSimulation_changingBST$round==1)],
								estimated_soc_reduction_75 = modelSimulation_changingBST$estimated_soc_reduction_75[which(modelSimulation_changingBST$round==1)],
								estimated_soc_reduction_low = modelSimulation_changingBST$estimated_soc_reduction_low[which(modelSimulation_changingBST$round==1)],
								taskDifficulty = modelSimulation_changingBST$taskDifficulty[which(modelSimulation_changingBST$round==1)],
								estimated_soc_reduction_high = modelSimulation_changingBST$estimated_soc_reduction_high[which(modelSimulation_changingBST$round==1)]
								)
subject_order_soc_reduction = order(modelSimulation_changingBST$soc_reduction[which(modelSimulation_changingBST$round==1)])
soc_reduction_recovery_data = soc_reduction_recovery_data[subject_order_soc_reduction,]
soc_reduction_recovery_data$subject = 1:nrow(soc_reduction_recovery_data)

soc_reduction_recovery_data$true_soc_change = NA
soc_reduction_recovery_data$true_soc_change[which(soc_reduction_recovery_data$true_soc_reduction<0)] = 'Decreasing'
soc_reduction_recovery_data$true_soc_change[which(soc_reduction_recovery_data$true_soc_reduction==0)] = 'Unchanged'
soc_reduction_recovery_data$true_soc_change[which(soc_reduction_recovery_data$true_soc_reduction>0)] = 'Increasing'
soc_reduction_recovery_data$estimated_soc_change = NA
soc_reduction_recovery_data$estimated_soc_change[which(soc_reduction_recovery_data$estimated_soc_reduction<0)] = 'Decreasing'
soc_reduction_recovery_data$estimated_soc_change[which(soc_reduction_recovery_data$estimated_soc_reduction==0)] = 'Unchanged'
soc_reduction_recovery_data$estimated_soc_change[which(soc_reduction_recovery_data$estimated_soc_reduction>0)] = 'Increasing'


soc_reduction_plot = ggplot(data=soc_reduction_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_soc_reduction_low, ymax=estimated_soc_reduction_high, colour=estimated_soc_change), width=.001) +
	geom_point(aes(subject, estimated_soc_reduction, colour=estimated_soc_change),size=.5)+
	#geom_point(aes(subject, bestfit_soc_reduction), colour='orange',size=.5)+
	geom_point(aes(subject, true_soc_reduction),size=.5) +
	geom_hline(yintercept=0, colour='black') +
	myTheme()+
	scale_colour_manual(values=c("Unchanged"="#999999",'Decreasing'="blue", "Increasing"="red"), name="Strength of \nfrequency dependence") +
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Agents', y=expression(paste('Individual ',delta,sep="")))

soc_reduction_plot2 = ggplot(data=soc_reduction_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_soc_reduction_low, ymax=estimated_soc_reduction_high), width=.001) +
	geom_point(aes(subject, estimated_soc_reduction), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_soc_reduction), colour='orange',size=.5)+
	geom_point(aes(subject, true_soc_reduction), colour='red',size=.5) +
	geom_hline(yintercept=0, colour='black') +
	myTheme()+
	#scale_colour_manual(values=c("Unchanged"="#999999",'Decreasing'="blue", "Increasing"="red"), name="Strength of \nfrequency dependence") +
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	ylim(c(-25,25))+
	labs(x='Agents', y=expression(paste('Individual ',delta,sep="")))


library(cowplot)
#plot_grid(alpha_plot, beta_plot, theta_plot, soc_raw_plot, soc_reduction_plot, soc_plot,  labels = c("","","","","",""), ncol = 2, align = 'v')
param_recov_plot = plot_grid(alpha_plot, beta_plot, annealing_plot, NULL, theta_plot, theta_slope_plot, soc_raw_plot, soc_reduction_plot2, labels = c("a","b","c","","d","e","f","g"), ncol = 2, align = 'v')
#param_recov_plot

ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs/param_recov_plot.pdf", plot = param_recov_plot, dpi = 600, width = 12, height = 9)

# soc (sigma)
modelSimulation_changingBST$taskDifficulty2 = 'Low Uncertainty'
modelSimulation_changingBST$taskDifficulty2[which(modelSimulation_changingBST$taskDifficulty=='50%')] = 'Moderate Uncertainty'
modelSimulation_changingBST$taskDifficulty2[which(modelSimulation_changingBST$taskDifficulty=='78%')] = 'High Uncertainty'
modelSimulation_changingBST$taskDifficulty2 = factor(modelSimulation_changingBST$taskDifficulty2, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))
true_sigma = ggplot(modelSimulation_changingBST,aes(round,soc,colour=theta))+
	geom_line(aes(group=subject)) + ylim(c(0,1))+
	facet_grid(.~taskDifficulty2)+
	labs(x='Round',y=bquote(atop('True social learning', 'weight' ~ sigma[i][t])))+
	scale_color_distiller(name=bquote(atop('True conformity', 'exponent' ~ theta[i][t])),palette='RdYlBu',direction=-1)+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
fit_sigma = ggplot(modelSimulation_changingBST,aes(round,estimated_soc,colour=theta))+
    geom_line(aes(group=subject)) + ylim(c(0,1))+ facet_grid(.~taskDifficulty2)+
    scale_color_distiller(name=bquote(atop('Fitted conformity', 'exponent' ~ theta[i][t])),palette='RdYlBu',direction=-1)+
    labs(x='Round',y=bquote(atop('Fitted social learning', 'weight' ~ sigma[i][t])))+
    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
sigma_plot = plot_grid(true_sigma, fit_sigma, labels = c("a","b"), ncol = 1, align = 'v')
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs/sigma_plot.pdf", plot = sigma_plot, dpi = 600, width = 9, height = 6)

# netBeta
true_netBeta = ggplot(modelSimulation_changingBST,aes(round,netBeta,colour=annealing))+
	geom_line(aes(group=subject)) + #ylim(c(0,6))+
	facet_grid(.~taskDifficulty2)+
	labs(x='Round',y=expression(paste('True inverse temperature ',beta,sep="")))+
	scale_color_distiller(name=expression(paste('',epsilon,sep="")),palette='RdYlBu',direction=-1)+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
fit_netBeta = ggplot(modelSimulation_changingBST,aes(round,estimated_netBeta,colour=estimated_annealing))+
    geom_line(aes(group=subject)) + #ylim(c(0,6))+
    facet_grid(.~taskDifficulty2)+
    scale_color_distiller(name=expression(paste('Fitted ',epsilon,sep="")),palette='RdYlBu',direction=-1)+
    labs(x='Round',y=expression(paste('Fitted inverse temperature ',beta,sep="")))+
    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
netBeta_plot = plot_grid(true_netBeta, fit_netBeta, labels = c("a","b"), ncol = 1, align = 'v')
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs/netBeta_plot.pdf", plot = netBeta_plot, dpi = 600, width = 9, height = 6)

# netTheta
modelSimulation_changingBST$netTheta = modelSimulation_changingBST$theta
true_netTheta = ggplot(modelSimulation_changingBST,aes(round,netTheta,colour=soc))+
	geom_line(aes(group=subject), alpha=2/3) + #ylim(c(0,6))+
	facet_grid(.~taskDifficulty2)+
	labs(x='Round',y=bquote(atop('True', 'conformity exponent' ~ theta[i][t])))+
	#scale_color_distiller(name=bquote(atop('True social learning weight' ~ sigma[i][t])), palette='RdYlBu',direction=-1)+
	scale_color_viridis(name=bquote(atop('True social learning', 'weight' ~ sigma[i][t])), option="magma", direction = -1)+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
fit_netTheta = ggplot(modelSimulation_changingBST,aes(round,estimated_netTheta,colour=estimated_soc))+
    geom_line(aes(group=subject), alpha=2/3) + #ylim(c(0,6))+
    facet_grid(.~taskDifficulty2)+
    #scale_color_distiller(name=bquote(atop('Fitted social learning weight' ~ sigma[i][t])), palette='RdYlBu',direction=-1)+
    scale_color_viridis(name=bquote(atop('True social learning', 'weight' ~ sigma[i][t])), option="magma", direction = -1)+
    labs(x='Round',y=bquote(atop('Fitted', 'conformity exponent' ~ theta[i][t])))+
    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
netTheta_plot = plot_grid(true_netTheta, fit_netTheta, labels = c("a","b"), ncol = 1, align = 'v')
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs/netTheta_plot.pdf", plot = netTheta_plot, dpi = 600, width = 9, height = 6)

## alpha
alpha_scatter = ggplot(subset(modelSimulation_changingBST, round==1),aes(alpha,estimated_alpha,colour=((alpha - estimated_alpha)^2)^(1/2)))+
    geom_point() + ylim(c(0,1))+ xlim(c(0,1))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',alpha,sep="")),y=expression(paste('Fitted ',alpha,sep="")))+
    annotate("text", x=0.15, y=0.90, label='italic(r) == 0.84', parse = TRUE, size = 5)+
    myTheme_legend()
## beta
beta_scatter = ggplot(subset(modelSimulation_changingBST, round==1),aes(beta,estimated_beta,colour=((beta - estimated_beta)^2)^(1/2)))+
    geom_point() + #ylim(c(0,5)) + xlim(c(0,5))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',beta[0],sep="")),y=expression(paste('Fitted ',beta[0],sep="")))+
    annotate("text", x=-2, y=3.5, label='italic(r) == 0.68', parse = TRUE, size = 5)+
    myTheme_legend()
## annealing
annealing_scatter = ggplot(subset(modelSimulation_changingBST, round==1),aes(annealing,estimated_annealing,colour=((annealing - estimated_annealing)^2)^(1/2)))+
    geom_point() + #ylim(c(0,5)) + xlim(c(0,5))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',epsilon,sep="")),y=expression(paste('Fitted ',epsilon,sep="")))+
    annotate("text", x=-4, y=5, label='italic(r) == 0.67', parse = TRUE, size = 5)+
    myTheme_legend()
## theta
theta_scatter = ggplot(subset(modelSimulation_changingBST, round==1),aes(theta0,estimated_theta,colour=((theta0 - estimated_theta)^2)^(1/2)))+
    geom_point() + #ylim(c(0,5)) + xlim(c(0,5))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',theta[0],sep="")),y=expression(paste('Fitted ',theta[0],sep="")))+
    annotate("text", x=-6, y=2.7, label='italic(r) == 0.45', parse = TRUE, size = 5)+
    myTheme_legend()
## theta_slope
theta_slope_scatter = ggplot(subset(modelSimulation_changingBST, round==1),aes(theta_slope,estimated_theta_slope,colour=((theta_slope - estimated_theta_slope)^2)^(1/2)))+
    geom_point() + #ylim(c(0,5)) + xlim(c(0,5))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',gamma,sep="")),y=expression(paste('Fitted ',gamma,sep="")))+
    annotate("text", x=-15, y=7, label='italic(r) == 0.51', parse = TRUE, size = 5)+
    myTheme_legend()
## soc_raw
soc_raw_scatter = ggplot(subset(modelSimulation_changingBST, round==1),aes(soc_raw,estimated_soc_raw,colour=((soc_raw - estimated_soc_raw)^2)^(1/2)))+
    geom_point() + #ylim(c(-3,3)) + xlim(c(-3,3))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',sigma[0],sep="")),y=expression(paste('Fitted ',sigma[0],sep="")))+
    annotate("text", x=-4, y=1, label='italic(r) == 0.61', parse = TRUE, size = 5)+
    myTheme_legend()
## soc_reduction
soc_reduction_scatter = ggplot(subset(modelSimulation_changingBST, round==1),aes(soc_reduc,estimated_soc_reduction,colour=((soc_reduc - estimated_soc_reduction)^2)^(1/2)))+
    geom_point() + #ylim(c(-6,0)) + xlim(c(-6,0))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',delta,sep="")),y=expression(paste('Fitted ',delta,sep="")))+
    annotate("text", x=-10, y=4, label='italic(r) == 0.70', parse = TRUE, size = 5)+
    myTheme_legend()


param_recov_scatter_plot = plot_grid(alpha_scatter, NULL, beta_scatter, annealing_scatter, theta_scatter, theta_slope_scatter, soc_raw_scatter, soc_reduction_scatter, labels = c("a","","b","c","d","e","f","g"), ncol = 2, align = 'v')
#param_recov_scatter_plot
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs/param_recov_scatter_plot.pdf", plot = param_recov_scatter_plot, dpi = 600, width = 9, height = 7)

## Correlation coefficients
cor.test(modelSimulation_changingBST$alpha, modelSimulation_changingBST$estimated_alpha)
cor.test(modelSimulation_changingBST$beta, modelSimulation_changingBST$estimated_beta)
cor.test(modelSimulation_changingBST$annealing, modelSimulation_changingBST$estimated_annealing)
cor.test(modelSimulation_changingBST$theta0, modelSimulation_changingBST$estimated_theta)
cor.test(modelSimulation_changingBST$theta_slope, modelSimulation_changingBST$estimated_theta_slope)
cor.test(modelSimulation_changingBST$soc_raw, modelSimulation_changingBST$estimated_soc_raw)
cor.test(modelSimulation_changingBST$soc_reduction, modelSimulation_changingBST$estimated_soc_reduction)


### MCMC diagnostics
for(i in 1:nrow(temp_changingBST_simulation)){
	temp_changingBST_simulation$log_lik_n_eff[i] <- summary(fit_changingBST_simulation)$summary[paste('log_lik[',i,']',sep=''),'n_eff']
	temp_changingBST_simulation$log_lik_Rhat[i] <- summary(fit_changingBST_simulation)$summary[paste('log_lik[',i,']',sep=''),'Rhat']
	temp_changingBST_simulation$alpha_n_eff[i] <- summary(fit_changingBST_simulation)$summary[paste('alpha[',i,']',sep=''),'n_eff']
	temp_changingBST_simulation$alpha_Rhat[i] <- summary(fit_changingBST_simulation)$summary[paste('alpha[',i,']',sep=''),'Rhat']
	temp_changingBST_simulation$beta_n_eff[i] <- summary(fit_changingBST_simulation)$summary[paste('beta[',i,']',sep=''),'n_eff']
	temp_changingBST_simulation$beta_Rhat[i] <- summary(fit_changingBST_simulation)$summary[paste('beta[',i,']',sep=''),'Rhat']
	temp_changingBST_simulation$soc_raw_n_eff[i] <- summary(fit_changingBST_simulation)$summary[paste('soc_raw[',i,']',sep=''),'n_eff']
	temp_changingBST_simulation$soc_raw_Rhat[i] <- summary(fit_changingBST_simulation)$summary[paste('soc_raw[',i,']',sep=''),'Rhat']
	temp_changingBST_simulation$soc_reduction_n_eff[i] <- summary(fit_changingBST_simulation)$summary[paste('soc_reduction[',i,']',sep=''),'n_eff']
	temp_changingBST_simulation$soc_reduction_Rhat[i] <- summary(fit_changingBST_simulation)$summary[paste('soc_reduction[',i,']',sep=''),'Rhat']
	temp_changingBST_simulation$theta_n_eff[i] <- summary(fit_changingBST_simulation)$summary[paste('theta[',i,']',sep=''),'n_eff']
	temp_changingBST_simulation$theta_Rhat[i] <- summary(fit_changingBST_simulation)$summary[paste('theta[',i,']',sep=''),'Rhat']
	temp_changingBST_simulation$theta_slope_n_eff[i] <- summary(fit_changingBST_simulation)$summary[paste('theta_slope[',i,']',sep=''),'n_eff']
	temp_changingBST_simulation$theta_slope_Rhat[i] <- summary(fit_changingBST_simulation)$summary[paste('theta_slope[',i,']',sep=''),'Rhat']
	temp_changingBST_simulation$annealing_n_eff[i] <- summary(fit_changingBST_simulation)$summary[paste('annealing[',i,']',sep=''),'n_eff']
	temp_changingBST_simulation$annealing_Rhat[i] <- summary(fit_changingBST_simulation)$summary[paste('annealing[',i,']',sep=''),'Rhat']
}
### Save the file again
write.csv(temp_changingBST_simulation,
			"~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/temp_changingBST_simulation.csv",
			row.names=FALSE)


### 25 - 75 ###
alpha_plot_50 = ggplot(data=alpha_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_alpha_25, ymax=estimated_alpha_75), width=.05) +
	geom_point(aes(subject, estimated_alpha), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_alpha), colour='orange',size=.5)+
	geom_point(aes(subject, true_alpha), colour='red',size=.5) +
	myTheme()+
	ylim(c(0,1)) +labs(x='Agents', y='alpha') #+ facet_grid(.~taskDifficulty)
beta_plot_50 = ggplot(data=beta_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_beta_25, ymax=estimated_beta_75), width=.05) +
	geom_point(aes(subject, estimated_beta), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_beta), colour='orange',size=.5)+
	geom_point(aes(subject, true_beta), colour='red',size=.5)+
	myTheme()+
	labs(x='Agents', y='beta')#+ facet_grid(.~taskDifficulty)
annealing_plot_50 = ggplot(data=annealing_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_annealing_25, ymax=estimated_annealing_75), width=.05) +
	geom_point(aes(subject, estimated_annealing), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_annealing), colour='orange',size=.5)+
	geom_point(aes(subject, true_annealing), colour='red',size=.5)+
	myTheme()+
	labs(x='Agents', y='annealing')#+ facet_grid(.~taskDifficulty)
theta_plot_50 = ggplot(data=theta_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_theta_25, ymax=estimated_theta_75), width=.05) +
	geom_point(aes(subject, estimated_theta), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_theta), colour='orange',size=.5)+
	geom_point(aes(subject, true_theta), colour='red',size=.5)+
	myTheme()+
	labs(x='Agents', y='theta')#+ facet_grid(.~taskDifficulty)
theta_slope_plot_50 = ggplot(data=theta_slope_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_theta_slope_25, ymax=estimated_theta_slope_75), width=.05) +
	geom_point(aes(subject, estimated_theta_slope), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_theta_slope), colour='orange',size=.5)+
	geom_point(aes(subject, true_theta_slope), colour='red',size=.5)+
	myTheme()+
	labs(x='Agents', y='theta_slope')#+ facet_grid(.~taskDifficulty)
soc_raw_plot_50 =ggplot(data=soc_raw_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_soc_raw_25, ymax=estimated_soc_raw_75), width=.05) +
	geom_point(aes(subject, estimated_soc_raw), colour='blue',size=.5)+
	#geom_point(aes(subject, bestfit_soc_raw), colour='orange',size=.5)+
	geom_point(aes(subject, true_soc_raw), colour='red',size=.5) +
	myTheme()+
	labs(x='Agents', y='soc_0')#+ facet_grid(.~taskDifficulty)
soc_reduction_plot_50 = ggplot(data=soc_reduction_recovery_data) +
	geom_errorbar(aes(x=subject, ymin=estimated_soc_reduction_25, ymax=estimated_soc_reduction_75, colour=estimated_soc_change), width=.05) +
	geom_point(aes(subject, estimated_soc_reduction, colour=estimated_soc_change),size=.5)+
	#geom_point(aes(subject, bestfit_soc_reduction), colour='orange',size=.5)+
	geom_point(aes(subject, true_soc_reduction),size=.5) +
	geom_hline(yintercept=0, colour='black') +
	myTheme()+
	scale_colour_manual(values=c("Unchanged"="#999999",'Decreasing'="blue", "Increasing"="red"), name="Strength of \nfrequency dependence") +
	labs(x='Agents', y='soc_change')#+ facet_grid(.~taskDifficulty)
param_recov_plot_2575 = plot_grid(alpha_plot_50, beta_plot_50, annealing_plot_50, theta_plot_50, theta_slope_plot_50, soc_raw_plot_50, soc_reduction_plot_50, labels = c("","","","","","",""), ncol = 2, align = 'v')
param_recov_plot_2575

ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs/param_recov_plot_2575.pdf", plot = param_recov_plot_2575, dpi = 600, width = 9, height = 6)


## Frequency dependence categorisation accuracy
theta_recovery_data$is_correct_95 = as.numeric(theta_recovery_data$true_mean_netTheta==theta_recovery_data$estimated_frequency_dependence)
theta_recovery_data$is_correct_50 = as.numeric(theta_recovery_data$true_mean_netTheta==theta_recovery_data$estimated_frequency_dependence2)
table(theta_recovery_data$is_correct_95)
table(theta_recovery_data$is_correct_50)

theta_recovery_data$true_mean_netTheta4 = theta_recovery_data$true_mean_netTheta
theta_recovery_data$true_mean_netTheta4[which(theta_recovery_data$true_mean_netTheta4!='neg-freq-dep')] = 'pos-freq-dep'

theta_recovery_data$estimated_frequency_dependence_point = "random-copying"
theta_recovery_data$estimated_frequency_dependence_point[which(theta_recovery_data$estimated_theta >= 0.2)] = "pos-freq-dep"
theta_recovery_data$estimated_frequency_dependence_point[which(theta_recovery_data$estimated_theta <= -0.2)] = "neg-freq-dep"


## Goodness of fit - global parameters
parameters_true_fitted = data.frame(
	params = rep(c('mu_alpha','mu_beta','mu_annealing','mu_soc','mu_soc_reduction','mu_theta','mu_theta_slope','s_alpha','s_beta','s_annealing','s_soc','s_soc_reduction','s_theta','s_theta_slope'), 3),
	type = rep(c(rep('mu',7),rep('variation',7)),3),
	taskDifficulty = c(rep('Low Uncertainty',14),rep('Moderate Uncertainty',14),rep('High Uncertainty',14)),
	true_values = c(
		0.97,1.97,3.17,-1.45,-1.35,1.31,1.09 ,1.86,1.49,1.49,0.77,0.82,1.56,1.39,
		0.95,2.08,2.00,-1.28,-2.34,0.39,8.24 ,1.58,0.74,1.77,1.62,2.50,2.36,3.86,
		0.62,1.80,2.16,-1.05,-3.21,0.97,4.03 ,2.58,1.01,1.99,1.35,3.86,2.24,6.81
		)
	#true_values = c(
	#	0.99,1.84,3.70,-1.55,-1.39,1.65,1.88,1.45,1.73,0.79,0.84,1.54,
	#	0.90,1.68,3.01,-2.37,-1.55,3.00,1.61,0.64,2.04,1.98,3.66,2.69,
	#	0.61,1.38,2.97,-2.16,-1.87,2.67,2.69,0.79,2.36,1.67,4.22,3.53
	#	)
	#true_values = c(
	#	0.64,0.5,0.5,-1.43,-1.34,1.47,0.27,0.1,0.1,0.39,0.29,0.95,
	#	0.60,0.4,0.4,-1.99,-1.37,2.65,0.29,0.1,0.1,1.39,2.50,1.47,
	#	0.57,0.3,0.3,-1.95,-1.46,2.33,0.34,0.1,0.1,1.15,3.30,2.25
	#	)
	)

parameters_true_fitted$fit_p2.5 = NA
parameters_true_fitted$fit_p25 = NA
parameters_true_fitted$fit_p50 = NA
parameters_true_fitted$fit_p75 = NA
parameters_true_fitted$fit_p97.5 = NA

for(p in names(table(parameters_true_fitted$params))){
	parameters_true_fitted$fit_p2.5[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_changingBST_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,1]
	parameters_true_fitted$fit_p25[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_changingBST_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,2]
	parameters_true_fitted$fit_p50[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_changingBST_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,3]
	parameters_true_fitted$fit_p75[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_changingBST_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,4]
	parameters_true_fitted$fit_p97.5[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_changingBST_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,5]
}

parameters_true_fitted$param_article_names = rep(c('mu_alpha','mu_beta','mu_epsilon','mu_sigma','mu_delta','mu_theta','mu_gamma','s_alpha','s_beta','s_epsilon','s_sigma','s_delta','s_theta','s_gamma'),3)
parameters_true_fitted$taskDifficulty = factor(parameters_true_fitted$taskDifficulty, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))


param_recov_global_plot = ggplot(parameters_true_fitted) +
	geom_errorbar(aes(x=param_article_names, ymin=fit_p2.5, ymax=fit_p97.5), colour='red',width=0.1)+
	geom_point(aes(param_article_names, fit_p50), colour='red', shape=2)+
	geom_point(aes(param_article_names, true_values))+
	myTheme()+
	theme(
		axis.text.x = element_text(angle=90),
		axis.title.y = element_text(angle=0))+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	labs(x='Parameters',y='',title='Parameter recovery test')+
	facet_grid(taskDifficulty~.)

ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs/param_recov_global_plot.pdf", plot = param_recov_global_plot, dpi = 600, width = 9, height = 9)


## ?% of parameters were recovered?
modelSimulation_changingBST$is_alpha_recov = 0
modelSimulation_changingBST$is_beta_recov = 0
modelSimulation_changingBST$is_annealing_recov = 0
modelSimulation_changingBST$is_soc_raw_recov = 0
modelSimulation_changingBST$is_soc_reduc_recov = 0
modelSimulation_changingBST$is_theta_recov = 0
modelSimulation_changingBST$is_theta_slope_recov = 0

modelSimulation_changingBST$is_alpha_recov[which(modelSimulation_changingBST$alpha>modelSimulation_changingBST$estimated_alpha_low & modelSimulation_changingBST$alpha<modelSimulation_changingBST$estimated_alpha_high)] <- 1
modelSimulation_changingBST$is_beta_recov[which(modelSimulation_changingBST$beta>modelSimulation_changingBST$estimated_beta_low & modelSimulation_changingBST$beta<modelSimulation_changingBST$estimated_beta_high)] <- 1
modelSimulation_changingBST$is_annealing_recov[which(modelSimulation_changingBST$annealing>modelSimulation_changingBST$estimated_annealing_low & modelSimulation_changingBST$annealing<modelSimulation_changingBST$estimated_annealing_high)] <- 1
modelSimulation_changingBST$is_soc_raw_recov[which(modelSimulation_changingBST$soc_raw>modelSimulation_changingBST$estimated_soc_raw_low & modelSimulation_changingBST$soc_raw<modelSimulation_changingBST$estimated_soc_raw_high)] <- 1
modelSimulation_changingBST$is_soc_reduc_recov[which(modelSimulation_changingBST$soc_reduc>modelSimulation_changingBST$estimated_soc_reduction_low & modelSimulation_changingBST$soc_reduc<modelSimulation_changingBST$estimated_soc_reduction_high)] <- 1
modelSimulation_changingBST$is_theta_recov[which(modelSimulation_changingBST$theta>modelSimulation_changingBST$estimated_theta_low & modelSimulation_changingBST$theta<modelSimulation_changingBST$estimated_theta_high)] <- 1
modelSimulation_changingBST$is_theta_slope_recov[which(modelSimulation_changingBST$theta_slope>modelSimulation_changingBST$estimated_theta_slope_low & modelSimulation_changingBST$theta_slope<modelSimulation_changingBST$estimated_theta_slope_high)] <- 1

length(which(subset(modelSimulation_changingBST, round==1)$is_alpha_recov==1))/572 * 100
length(which(subset(modelSimulation_changingBST, round==1)$is_beta_recov==1))/572 * 100
length(which(subset(modelSimulation_changingBST, round==1)$is_annealing_recov==1))/572 * 100
length(which(subset(modelSimulation_changingBST, round==1)$is_soc_raw_recov==1))/572 * 100
length(which(subset(modelSimulation_changingBST, round==1)$is_soc_reduc_recov==1))/572 * 100
length(which(subset(modelSimulation_changingBST, round==1)$is_theta_recov==1))/572 * 100
length(which(subset(modelSimulation_changingBST, round==1)$is_theta_slope_recov==1))/572 * 100
