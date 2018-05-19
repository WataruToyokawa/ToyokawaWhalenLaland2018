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
##                         UNC6_sReduc_annealing model
##
################################################################################
modelSimulation_UNC6_sReduc_annealing_25 <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/modelSimulation_UNC6_sReduc_annealing_25.csv", col_names = FALSE)
modelSimulation_UNC6_sReduc_annealing_50 <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/modelSimulation_UNC6_sReduc_annealing_50.csv", col_names = FALSE)
modelSimulation_UNC6_sReduc_annealing_78 <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/modelSimulation_UNC6_sReduc_annealing_78.csv", col_names = FALSE)

modelSimulation_UNC6_sReduc_annealing = rbind(rbind(modelSimulation_UNC6_sReduc_annealing_25,modelSimulation_UNC6_sReduc_annealing_50),modelSimulation_UNC6_sReduc_annealing_78)
colnames(modelSimulation_UNC6_sReduc_annealing) <- c('Rep','groupSize','alpha','beta','annealing','netBeta','soc_raw','soc_reduc','soc','theta','Y','outcome','performance','socialFreq0','socialFreq1','socialFreq2')
modelSimulation_UNC6_sReduc_annealing$taskDifficulty = c(rep('25%',nrow(modelSimulation_UNC6_sReduc_annealing_25)),rep('50%',nrow(modelSimulation_UNC6_sReduc_annealing_50)),rep('78%',nrow(modelSimulation_UNC6_sReduc_annealing_78)))
modelSimulation_UNC6_sReduc_annealing$round = rep(1:70, nrow(modelSimulation_UNC6_sReduc_annealing)/70)
modelSimulation_UNC6_sReduc_annealing$subject = NA
for(t in 1:70){
	modelSimulation_UNC6_sReduc_annealing$subject[which(modelSimulation_UNC6_sReduc_annealing$round==t)] = 1:(nrow(modelSimulation_UNC6_sReduc_annealing)/70)
}
modelSimulation_UNC6_sReduc_annealing$optimalChoice = modelSimulation_UNC6_sReduc_annealing$performance
modelSimulation_UNC6_sReduc_annealing$optimalChoice[which(modelSimulation_UNC6_sReduc_annealing$optimalChoice==1&modelSimulation_UNC6_sReduc_annealing$round>40)] <- 0
modelSimulation_UNC6_sReduc_annealing$optimalChoice[which(modelSimulation_UNC6_sReduc_annealing$optimalChoice==2)] <- 1
modelSimulation_UNC6_sReduc_annealing$con = c(rep(1,nrow(modelSimulation_UNC6_sReduc_annealing_25)),rep(2,nrow(modelSimulation_UNC6_sReduc_annealing_50)),rep(3,nrow(modelSimulation_UNC6_sReduc_annealing_78)))

ggplot(modelSimulation_UNC6_sReduc_annealing,aes(outcome,fill=as.factor(Y)))+geom_histogram()+facet_grid(.~taskDifficulty)
ggplot(modelSimulation_UNC6_sReduc_annealing,aes(round,performance,colour=groupSize)) +
	stat_summary(aes(group=groupSize), fun.y = mean, geom="line") + ylim(c(0,2))+
	facet_grid(.~taskDifficulty)
ggplot(modelSimulation_UNC6_sReduc_annealing,aes(round,optimalChoice,colour=groupSize)) +
	stat_summary(aes(group=groupSize), fun.y = mean, geom="line") + ylim(c(0,1))+
	facet_grid(.~taskDifficulty)

#ggplot(modelSimulation_UNC6_sReduc_annealing)+geom_line(aes(round,Y))+geom_point(aes(round,Y),alpha=1/4)+facet_grid(groupSize+Rep~taskDifficulty)

ReinforcementLearningStanData_simulation = list(
	All = nrow(modelSimulation_UNC6_sReduc_annealing),
	Nsub = length(table(modelSimulation_UNC6_sReduc_annealing$subject)),
	Ncue = 3, # number of options
	Ntrial = 70,
	Ncon = 3,
	sub = modelSimulation_UNC6_sReduc_annealing$subject,
	Y = modelSimulation_UNC6_sReduc_annealing$Y,
	trial = modelSimulation_UNC6_sReduc_annealing$round,
	outcome = modelSimulation_UNC6_sReduc_annealing$outcome,
	trueAlpha=modelSimulation_UNC6_sReduc_annealing$alpha,
	trueBeta = modelSimulation_UNC6_sReduc_annealing$beta,
	trueAnnealing = modelSimulation_UNC6_sReduc_annealing$annealing,
	trueNetBeta=modelSimulation_UNC6_sReduc_annealing$netBeta,
	trueSoc=modelSimulation_UNC6_sReduc_annealing$soc,
	trueSocRaw=modelSimulation_UNC6_sReduc_annealing$soc_raw,
	trueSocReduc=modelSimulation_UNC6_sReduc_annealing$soc_reduc,
	trueTheta=modelSimulation_UNC6_sReduc_annealing$theta
	)
con = numeric(length(table(modelSimulation_UNC6_sReduc_annealing$subject)))
for(i in 1:ReinforcementLearningStanData_simulation$Nsub){
	con[i] <- modelSimulation_UNC6_sReduc_annealing$con[which(modelSimulation_UNC6_sReduc_annealing$subject==i)][1]
}
ReinforcementLearningStanData_simulation$condition = con
startT = c()
MaxTrial = c()
for(i in 1:length(table(modelSimulation_UNC6_sReduc_annealing$subject))){
	startT = append(startT, min(subset(modelSimulation_UNC6_sReduc_annealing, subject==i)$round))
	MaxTrial = append(MaxTrial, length(which(subset(modelSimulation_UNC6_sReduc_annealing, subject==i)$Y>=1)))
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
F0_sim = array(rep(NA,nrow(modelSimulation_UNC6_sReduc_annealing)), c(ReinforcementLearningStanData_simulation$Nsub, ReinforcementLearningStanData_simulation$Ncue, ReinforcementLearningStanData_simulation$Ntrial))
F = F_calculation_simulation.compiled(F0_sim,ReinforcementLearningStanData_simulation$Nsub,ReinforcementLearningStanData_simulation$Ntrial,modelSimulation_UNC6_sReduc_annealing)
ReinforcementLearningStanData_simulation$F = round(F)


################################################################################
##
##                FITTING: UNC6_sReduc_annealing model -- STANDARD
##
################################################################################

	## FITTING
model_UNC6_sReduc_annealing = stan_model(file="~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model_UNC6_sReduc_annealing.stan")
fit_UNC6_sReduc_annealing_simulation = sampling(
	model_UNC6_sReduc_annealing, data=ReinforcementLearningStanData_simulation, seed=77,
	#control = list(adapt_delta = 0.9999, max_treedepth =20, stepsize = 0.001),
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
fit_UNC6_sReduc_annealing_simulation
saveRDS(fit_UNC6_sReduc_annealing_simulation, file='~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/mcmc_result/fit_UNC6_sReduc_annealing_simulation.rds')
ms_UNC6_sReduc_annealing_simulation <- rstan::extract(fit_UNC6_sReduc_annealing_simulation)
plot(fit_UNC6_sReduc_annealing_simulation, pars=c('mu_alpha','s_alpha','mu_beta','s_beta','mu_soc','s_soc','mu_theta','s_theta','mu_soc_reduction','s_soc_reduction','mu_annealing','s_annealing'))

	## 事後診断
traceplot(fit_UNC6_sReduc_annealing_simulation, pars=c('mu_alpha','mu_beta','mu_annealing'), inc_warmup=TRUE)
traceplot(fit_UNC6_sReduc_annealing_simulation, pars=c('mu_theta','mu_soc','mu_soc_reduction'), inc_warmup=TRUE)
traceplot(fit_UNC6_sReduc_annealing_simulation, pars=c('s_alpha','s_beta','s_theta','s_soc','s_soc_reduction'), inc_warmup=TRUE)

temp_UNC6_sReduc_annealing_simulation = data.frame(WAICi_UNC6_sReduc_annealing = WAIC_indv(fit_UNC6_sReduc_annealing_simulation)$waic)

temp_UNC6_sReduc_annealing_simulation <- cbind(temp_UNC6_sReduc_annealing_simulation, as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$alpha[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing_simulation <- cbind(temp_UNC6_sReduc_annealing_simulation, as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$beta[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing_simulation <- cbind(temp_UNC6_sReduc_annealing_simulation, as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$annealing[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing_simulation <- cbind(temp_UNC6_sReduc_annealing_simulation, as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$theta[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing_simulation <- cbind(temp_UNC6_sReduc_annealing_simulation, as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$soc_reduction[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
temp_UNC6_sReduc_annealing_simulation <- cbind(temp_UNC6_sReduc_annealing_simulation, as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$soc_raw[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(temp_UNC6_sReduc_annealing_simulation) <- c("UNC6_sReduc_annealing_simulation_WAICi",
							"UNC6_sReduc_annealing_simulation_alpha_p2.5", "UNC6_sReduc_annealing_simulation_alpha_p25", "UNC6_sReduc_annealing_simulation_alpha_p50", "UNC6_sReduc_annealing_simulation_alpha_p75", "UNC6_sReduc_annealing_simulation_alpha_p97.5",
							"UNC6_sReduc_annealing_simulation_beta_p2.5", "UNC6_sReduc_annealing_simulation_beta_p25", "UNC6_sReduc_annealing_simulation_beta_p50", "UNC6_sReduc_annealing_simulation_beta_p75", "UNC6_sReduc_annealing_simulation_beta_p97.5",
							"UNC6_sReduc_annealing_simulation_annealing_p2.5", "UNC6_sReduc_annealing_simulation_annealing_p25", "UNC6_sReduc_annealing_simulation_annealing_p50", "UNC6_sReduc_annealing_simulation_annealing_p75", "UNC6_sReduc_annealing_simulation_annealing_p97.5",
							"UNC6_sReduc_annealing_simulation_theta_p2.5", "UNC6_sReduc_annealing_simulation_theta_p25", "UNC6_sReduc_annealing_simulation_theta_p50", "UNC6_sReduc_annealing_simulation_theta_p75", "UNC6_sReduc_annealing_simulation_theta_p97.5",
							"UNC6_sReduc_annealing_simulation_soc_reduction_p2.5", "UNC6_sReduc_annealing_simulation_soc_reduction_p25", "UNC6_sReduc_annealing_simulation_soc_reduction_p50", "UNC6_sReduc_annealing_simulation_soc_reduction_p75", "UNC6_sReduc_annealing_simulation_soc_reduction_p97.5",
							"UNC6_sReduc_annealing_simulation_soc_raw_p2.5", "UNC6_sReduc_annealing_simulation_soc_raw_p25", "UNC6_sReduc_annealing_simulation_soc_raw_p50", "UNC6_sReduc_annealing_simulation_soc_raw_p75", "UNC6_sReduc_annealing_simulation_soc_raw_p97.5")
max_log_lik = apply(ms_UNC6_sReduc_annealing_simulation$log_lik[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)], 2, which.max)
max_UNC6_sReduc_annealing_simulation_alpha=rep(NA, ncol(ms_UNC6_sReduc_annealing_simulation$log_lik))
max_UNC6_sReduc_annealing_simulation_beta=rep(NA, ncol(ms_UNC6_sReduc_annealing_simulation$log_lik))
max_UNC6_sReduc_annealing_simulation_annealing=rep(NA, ncol(ms_UNC6_sReduc_annealing_simulation$log_lik))
max_UNC6_sReduc_annealing_simulation_theta=rep(NA, ncol(ms_UNC6_sReduc_annealing_simulation$log_lik))
max_UNC6_sReduc_annealing_simulation_soc_reduction=rep(NA, ncol(ms_UNC6_sReduc_annealing_simulation$log_lik))
max_UNC6_sReduc_annealing_simulation_soc_raw=rep(NA, ncol(ms_UNC6_sReduc_annealing_simulation$log_lik))
for(i in 1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)){
	max_UNC6_sReduc_annealing_simulation_alpha[i] <- ms_UNC6_sReduc_annealing_simulation$alpha[max_log_lik[i],i]
	max_UNC6_sReduc_annealing_simulation_beta[i] <- ms_UNC6_sReduc_annealing_simulation$beta[max_log_lik[i],i]
	max_UNC6_sReduc_annealing_simulation_annealing[i] <- ms_UNC6_sReduc_annealing_simulation$annealing[max_log_lik[i],i]
	max_UNC6_sReduc_annealing_simulation_theta[i] <- ms_UNC6_sReduc_annealing_simulation$theta[max_log_lik[i],i]
	max_UNC6_sReduc_annealing_simulation_soc_reduction[i] <- ms_UNC6_sReduc_annealing_simulation$soc_reduction[max_log_lik[i],i]
	max_UNC6_sReduc_annealing_simulation_soc_raw[i] <- ms_UNC6_sReduc_annealing_simulation$soc_raw[max_log_lik[i],i]
}
temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_alpha = max_UNC6_sReduc_annealing_simulation_alpha
temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_beta = max_UNC6_sReduc_annealing_simulation_beta
temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_annealing = max_UNC6_sReduc_annealing_simulation_annealing
temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_theta = max_UNC6_sReduc_annealing_simulation_theta
temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_soc_reduction = max_UNC6_sReduc_annealing_simulation_soc_reduction
temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_soc_raw = max_UNC6_sReduc_annealing_simulation_soc_raw
temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_log_lik = apply(ms_UNC6_sReduc_annealing_simulation$log_lik[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)], 2, max)
numParam = 6
temp_UNC6_sReduc_annealing_simulation$AIC_UNC6_sReduc_annealing = -temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_log_lik + numParam ## 正確には AIC/2
temp_UNC6_sReduc_annealing_simulation$BIC_UNC6_sReduc_annealing = -temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_log_lik + (numParam/2)*log(70) ## more precisely, BIC/2
temp_UNC6_sReduc_annealing_simulation$WAICi_UNC6_sReduc_annealing = WAIC_indv(fit_UNC6_sReduc_annealing_simulation)$waic


temp_UNC6_sReduc_annealing_simulation$log_lik_n_eff = NA
temp_UNC6_sReduc_annealing_simulation$log_lik_Rhat = NA
temp_UNC6_sReduc_annealing_simulation$alpha_n_eff = NA
temp_UNC6_sReduc_annealing_simulation$alpha_Rhat = NA
temp_UNC6_sReduc_annealing_simulation$beta_n_eff = NA
temp_UNC6_sReduc_annealing_simulation$beta_Rhat = NA
temp_UNC6_sReduc_annealing_simulation$soc_raw_n_eff = NA
temp_UNC6_sReduc_annealing_simulation$soc_raw_Rhat = NA
temp_UNC6_sReduc_annealing_simulation$soc_reduction_n_eff = NA
temp_UNC6_sReduc_annealing_simulation$soc_reduction_Rhat = NA
temp_UNC6_sReduc_annealing_simulation$theta_n_eff = NA
temp_UNC6_sReduc_annealing_simulation$theta_Rhat = NA
temp_UNC6_sReduc_annealing_simulation$annealing_n_eff = NA
temp_UNC6_sReduc_annealing_simulation$annealing_Rhat = NA


	## Estimated copying rate at each round
soc_table_UNC6_sReduc_annealing_simulation = data.frame(subject = 1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik), round = rep(1, ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)))
soc_table_UNC6_sReduc_annealing_simulation = cbind(soc_table_UNC6_sReduc_annealing_simulation, as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$soc[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(soc_table_UNC6_sReduc_annealing_simulation) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
	soc_table_UNC6_sReduc_annealing_simulation = rbind(soc_table_UNC6_sReduc_annealing_simulation, cbind(cbind(1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik), rep(t,ncol(ms_UNC6_sReduc_annealing_simulation$log_lik))), as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$soc[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(soc_table_UNC6_sReduc_annealing_simulation) <- c("subject", "round", "UNC6_sReduc_annealing_soc_p2.5", "UNC6_sReduc_annealing_soc_p25", "UNC6_sReduc_annealing_soc_p50", "UNC6_sReduc_annealing_soc_p75", "UNC6_sReduc_annealing_soc_p97.5")

write.csv(soc_table_UNC6_sReduc_annealing_simulation,
			"~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/soc_table_UNC6_sReduc_annealing_simulation.csv",
			row.names=FALSE)

	## Estimated netBeta at each round
netBeta_table_UNC6_sReduc_annealing_simulation = data.frame(subject = 1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik), round = rep(1, ncol(ms_UNC6_sReduc_annealing_simulation$log_lik)))
netBeta_table_UNC6_sReduc_annealing_simulation = cbind(netBeta_table_UNC6_sReduc_annealing_simulation, as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$netBeta[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik),1], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))))
colnames(netBeta_table_UNC6_sReduc_annealing_simulation) <- c("1","2","2.5%","25%","50%","75%","97.5%")
for(t in 2:70) {
	netBeta_table_UNC6_sReduc_annealing_simulation = rbind(netBeta_table_UNC6_sReduc_annealing_simulation, cbind(cbind(1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik), rep(t,ncol(ms_UNC6_sReduc_annealing_simulation$log_lik))), as.data.frame(t(apply(ms_UNC6_sReduc_annealing_simulation$netBeta[,1:ncol(ms_UNC6_sReduc_annealing_simulation$log_lik),t], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975))))))
}
colnames(netBeta_table_UNC6_sReduc_annealing_simulation) <- c("subject", "round", "UNC6_sReduc_annealing_netBeta_p2.5", "UNC6_sReduc_annealing_netBeta_p25", "UNC6_sReduc_annealing_netBeta_p50", "UNC6_sReduc_annealing_netBeta_p75", "UNC6_sReduc_annealing_netBeta_p97.5")

write.csv(netBeta_table_UNC6_sReduc_annealing_simulation,
			"~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/netBeta_table_UNC6_sReduc_annealing_simulation.csv", row.names=FALSE)


## Comparing with the TRUE values
## Were parameter values correctly estimated?
modelSimulation_UNC6_sReduc_annealing$estimated_alpha = NA
modelSimulation_UNC6_sReduc_annealing$estimated_alpha_high = NA
modelSimulation_UNC6_sReduc_annealing$estimated_alpha_low = NA
modelSimulation_UNC6_sReduc_annealing$estimated_alpha_25 = NA
modelSimulation_UNC6_sReduc_annealing$estimated_alpha_75 = NA

modelSimulation_UNC6_sReduc_annealing$estimated_beta = NA
modelSimulation_UNC6_sReduc_annealing$estimated_beta_high = NA
modelSimulation_UNC6_sReduc_annealing$estimated_beta_low = NA
modelSimulation_UNC6_sReduc_annealing$estimated_beta_25 = NA
modelSimulation_UNC6_sReduc_annealing$estimated_beta_75 = NA

modelSimulation_UNC6_sReduc_annealing$estimated_annealing = NA
modelSimulation_UNC6_sReduc_annealing$estimated_annealing_high = NA
modelSimulation_UNC6_sReduc_annealing$estimated_annealing_low = NA
modelSimulation_UNC6_sReduc_annealing$estimated_annealing_25 = NA
modelSimulation_UNC6_sReduc_annealing$estimated_annealing_75 = NA

modelSimulation_UNC6_sReduc_annealing$estimated_netBeta = NA
modelSimulation_UNC6_sReduc_annealing$estimated_netBeta_high = NA
modelSimulation_UNC6_sReduc_annealing$estimated_netBeta_low = NA
modelSimulation_UNC6_sReduc_annealing$estimated_netBeta_25 = NA
modelSimulation_UNC6_sReduc_annealing$estimated_netBeta_75 = NA

modelSimulation_UNC6_sReduc_annealing$estimated_theta = NA
modelSimulation_UNC6_sReduc_annealing$estimated_theta_high = NA
modelSimulation_UNC6_sReduc_annealing$estimated_theta_low = NA
modelSimulation_UNC6_sReduc_annealing$estimated_theta_25 = NA
modelSimulation_UNC6_sReduc_annealing$estimated_theta_75 = NA

modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_high = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_low = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_25 = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_75 = NA

modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_high = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_low = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_25 = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_75 = NA

modelSimulation_UNC6_sReduc_annealing$estimated_soc = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_high = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_low = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_25 = NA
modelSimulation_UNC6_sReduc_annealing$estimated_soc_75 = NA

modelSimulation_UNC6_sReduc_annealing$bestfit_alpha = NA
modelSimulation_UNC6_sReduc_annealing$bestfit_beta = NA
modelSimulation_UNC6_sReduc_annealing$bestfit_annealing = NA
modelSimulation_UNC6_sReduc_annealing$bestfit_netBeta = NA
modelSimulation_UNC6_sReduc_annealing$bestfit_theta = NA
modelSimulation_UNC6_sReduc_annealing$bestfit_soc_raw = NA
modelSimulation_UNC6_sReduc_annealing$bestfit_soc_reduction = NA
modelSimulation_UNC6_sReduc_annealing$bestfit_soc = NA

for (i in 1:nrow(temp_UNC6_sReduc_annealing_simulation)) {
        # alpha
    modelSimulation_UNC6_sReduc_annealing$bestfit_alpha[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_alpha[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_alpha[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_alpha_p50[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_alpha_low[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_alpha_p2.5[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_alpha_high[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_alpha_p97.5[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_alpha_25[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_alpha_p25[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_alpha_75[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_alpha_p75[i]
        # beta
    modelSimulation_UNC6_sReduc_annealing$estimated_beta[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_beta_p50[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_beta_low[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_beta_p2.5[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_beta_high[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_beta_p97.5[i]
    modelSimulation_UNC6_sReduc_annealing$bestfit_beta[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_beta[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_beta_25[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_beta_p25[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_beta_75[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_beta_p75[i]
        # annealing
    modelSimulation_UNC6_sReduc_annealing$estimated_annealing[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_annealing_p50[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_annealing_low[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_annealing_p2.5[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_annealing_high[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_annealing_p97.5[i]
    modelSimulation_UNC6_sReduc_annealing$bestfit_annealing[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_annealing[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_annealing_25[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_annealing_p25[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_annealing_75[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_annealing_p75[i]
        # theta
    modelSimulation_UNC6_sReduc_annealing$estimated_theta[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_theta_p50[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_theta_low[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_theta_p2.5[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_theta_high[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_theta_p97.5[i]
    modelSimulation_UNC6_sReduc_annealing$bestfit_theta[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_theta[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_theta_25[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_theta_p25[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_theta_75[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_theta_p75[i]
        # soc_raw
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_raw_p50[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_low[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_raw_p2.5[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_high[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_raw_p97.5[i]
    modelSimulation_UNC6_sReduc_annealing$bestfit_soc_raw[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_soc_raw[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_25[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_raw_p25[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_75[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_raw_p75[i]
        # soc_reduction
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_reduction_p50[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_low[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_reduction_p2.5[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_high[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_reduction_p97.5[i]
    modelSimulation_UNC6_sReduc_annealing$bestfit_soc_reduction[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$max_UNC6_sReduc_annealing_simulation_soc_reduction[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_25[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_reduction_p25[i]
    modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_75[which(modelSimulation_UNC6_sReduc_annealing$subject==i)] <- temp_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_simulation_soc_reduction_p75[i]
}

for (i in 1:nrow(temp_UNC6_sReduc_annealing_simulation)) {
	for (t in 1:70) {
		# soc
		modelSimulation_UNC6_sReduc_annealing$estimated_soc[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =soc_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_soc_p50[which(soc_table_UNC6_sReduc_annealing_simulation$subject==i&soc_table_UNC6_sReduc_annealing_simulation$round==t)]
		modelSimulation_UNC6_sReduc_annealing$estimated_soc_low[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =soc_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_soc_p2.5[which(soc_table_UNC6_sReduc_annealing_simulation$subject==i&soc_table_UNC6_sReduc_annealing_simulation$round==t)]
		modelSimulation_UNC6_sReduc_annealing$estimated_soc_high[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =soc_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_soc_p97.5[which(soc_table_UNC6_sReduc_annealing_simulation$subject==i&soc_table_UNC6_sReduc_annealing_simulation$round==t)]
		modelSimulation_UNC6_sReduc_annealing$estimated_soc_25[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =soc_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_soc_p25[which(soc_table_UNC6_sReduc_annealing_simulation$subject==i&soc_table_UNC6_sReduc_annealing_simulation$round==t)]
		modelSimulation_UNC6_sReduc_annealing$estimated_soc_75[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =soc_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_soc_p75[which(soc_table_UNC6_sReduc_annealing_simulation$subject==i&soc_table_UNC6_sReduc_annealing_simulation$round==t)]
		# netBeta
		modelSimulation_UNC6_sReduc_annealing$estimated_netBeta[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =netBeta_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_netBeta_p50[which(netBeta_table_UNC6_sReduc_annealing_simulation$subject==i&netBeta_table_UNC6_sReduc_annealing_simulation$round==t)]
		modelSimulation_UNC6_sReduc_annealing$estimated_netBeta_low[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =netBeta_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_netBeta_p2.5[which(netBeta_table_UNC6_sReduc_annealing_simulation$subject==i&netBeta_table_UNC6_sReduc_annealing_simulation$round==t)]
		modelSimulation_UNC6_sReduc_annealing$estimated_netBeta_high[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =netBeta_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_netBeta_p97.5[which(netBeta_table_UNC6_sReduc_annealing_simulation$subject==i&netBeta_table_UNC6_sReduc_annealing_simulation$round==t)]
		modelSimulation_UNC6_sReduc_annealing$estimated_netBeta_25[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =netBeta_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_netBeta_p25[which(netBeta_table_UNC6_sReduc_annealing_simulation$subject==i&netBeta_table_UNC6_sReduc_annealing_simulation$round==t)]
		modelSimulation_UNC6_sReduc_annealing$estimated_netBeta_75[which(modelSimulation_UNC6_sReduc_annealing$subject==i&modelSimulation_UNC6_sReduc_annealing$round==t)] =netBeta_table_UNC6_sReduc_annealing_simulation$UNC6_sReduc_annealing_netBeta_p75[which(netBeta_table_UNC6_sReduc_annealing_simulation$subject==i&netBeta_table_UNC6_sReduc_annealing_simulation$round==t)]
	}
}


###################################################################
##
## Parameter recovery test
##
###################################################################

## Alpha
alpha_recovery_data = data.frame(
								true_alpha = modelSimulation_UNC6_sReduc_annealing$alpha[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_alpha = modelSimulation_UNC6_sReduc_annealing$estimated_alpha[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								bestfit_alpha = modelSimulation_UNC6_sReduc_annealing$bestfit_alpha[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_alpha_low = modelSimulation_UNC6_sReduc_annealing$estimated_alpha_low[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								taskDifficulty = modelSimulation_UNC6_sReduc_annealing$taskDifficulty[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_alpha_25 = modelSimulation_UNC6_sReduc_annealing$estimated_alpha_25[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_alpha_75 = modelSimulation_UNC6_sReduc_annealing$estimated_alpha_75[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_alpha_high = modelSimulation_UNC6_sReduc_annealing$estimated_alpha_high[which(modelSimulation_UNC6_sReduc_annealing$round==1)]
								)
subject_order_alpha = order(modelSimulation_UNC6_sReduc_annealing$alpha[which(modelSimulation_UNC6_sReduc_annealing$round==1)])
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
								true_beta = modelSimulation_UNC6_sReduc_annealing$beta[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_beta = modelSimulation_UNC6_sReduc_annealing$estimated_beta[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								bestfit_beta = modelSimulation_UNC6_sReduc_annealing$bestfit_beta[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_beta_low = modelSimulation_UNC6_sReduc_annealing$estimated_beta_low[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								taskDifficulty = modelSimulation_UNC6_sReduc_annealing$taskDifficulty[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_beta_25 = modelSimulation_UNC6_sReduc_annealing$estimated_beta_25[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_beta_75 = modelSimulation_UNC6_sReduc_annealing$estimated_beta_75[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_beta_high = modelSimulation_UNC6_sReduc_annealing$estimated_beta_high[which(modelSimulation_UNC6_sReduc_annealing$round==1)]
								)
subject_order_beta = order(modelSimulation_UNC6_sReduc_annealing$beta[which(modelSimulation_UNC6_sReduc_annealing$round==1)])
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
								true_annealing = modelSimulation_UNC6_sReduc_annealing$annealing[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_annealing = modelSimulation_UNC6_sReduc_annealing$estimated_annealing[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								bestfit_annealing = modelSimulation_UNC6_sReduc_annealing$bestfit_annealing[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_annealing_low = modelSimulation_UNC6_sReduc_annealing$estimated_annealing_low[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								taskDifficulty = modelSimulation_UNC6_sReduc_annealing$taskDifficulty[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_annealing_25 = modelSimulation_UNC6_sReduc_annealing$estimated_annealing_25[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_annealing_75 = modelSimulation_UNC6_sReduc_annealing$estimated_annealing_75[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_annealing_high = modelSimulation_UNC6_sReduc_annealing$estimated_annealing_high[which(modelSimulation_UNC6_sReduc_annealing$round==1)]
								)
subject_order_annealing = order(modelSimulation_UNC6_sReduc_annealing$annealing[which(modelSimulation_UNC6_sReduc_annealing$round==1)])
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

## theta
theta_recovery_data = data.frame(
								true_theta = modelSimulation_UNC6_sReduc_annealing$theta[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_theta = modelSimulation_UNC6_sReduc_annealing$estimated_theta[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								bestfit_theta = modelSimulation_UNC6_sReduc_annealing$bestfit_theta[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_theta_low = modelSimulation_UNC6_sReduc_annealing$estimated_theta_low[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								taskDifficulty = modelSimulation_UNC6_sReduc_annealing$taskDifficulty[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_theta_25 = modelSimulation_UNC6_sReduc_annealing$estimated_theta_25[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_theta_75 = modelSimulation_UNC6_sReduc_annealing$estimated_theta_75[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_theta_high = modelSimulation_UNC6_sReduc_annealing$estimated_theta_high[which(modelSimulation_UNC6_sReduc_annealing$round==1)]
								)
theta_recovery_data$true_frequency_dependence = NA
theta_recovery_data$true_frequency_dependence[which(theta_recovery_data$true_theta<0)] = 'neg-freq-dep'
theta_recovery_data$true_frequency_dependence[which(theta_recovery_data$true_theta>=0&theta_recovery_data$true_theta<=1)] = 'weak-pos-freq-dep'
theta_recovery_data$true_frequency_dependence[which(theta_recovery_data$true_theta>1)] = 'conformist'
	# rand, neg, weak-pos, conf
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


subject_order_theta = order(modelSimulation_UNC6_sReduc_annealing$theta[which(modelSimulation_UNC6_sReduc_annealing$round==1)])
theta_recovery_data = theta_recovery_data[subject_order_theta,]
theta_recovery_data$subject = 1:nrow(theta_recovery_data)

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
								true_soc_raw = modelSimulation_UNC6_sReduc_annealing$soc_raw[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_raw = modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								bestfit_soc_raw = modelSimulation_UNC6_sReduc_annealing$bestfit_soc_raw[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								taskDifficulty = modelSimulation_UNC6_sReduc_annealing$taskDifficulty[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_raw_25 = modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_25[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_raw_75 = modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_75[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_raw_low = modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_low[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_raw_high = modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_high[which(modelSimulation_UNC6_sReduc_annealing$round==1)]
								)
subject_order_soc_raw = order(modelSimulation_UNC6_sReduc_annealing$soc_raw[which(modelSimulation_UNC6_sReduc_annealing$round==1)])
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
modelSimulation_UNC6_sReduc_annealing$soc_reduction = modelSimulation_UNC6_sReduc_annealing$soc_reduc
soc_reduction_recovery_data = data.frame(
								true_soc_reduction = modelSimulation_UNC6_sReduc_annealing$soc_reduction[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_reduction = modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								bestfit_soc_reduction = modelSimulation_UNC6_sReduc_annealing$bestfit_soc_reduction[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_reduction_25 = modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_25[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_reduction_75 = modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_75[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_reduction_low = modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_low[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								taskDifficulty = modelSimulation_UNC6_sReduc_annealing$taskDifficulty[which(modelSimulation_UNC6_sReduc_annealing$round==1)],
								estimated_soc_reduction_high = modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_high[which(modelSimulation_UNC6_sReduc_annealing$round==1)]
								)
subject_order_soc_reduction = order(modelSimulation_UNC6_sReduc_annealing$soc_reduction[which(modelSimulation_UNC6_sReduc_annealing$round==1)])
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
	labs(x='Agents', y=expression(paste('Individual ',delta,sep="")))


library(cowplot)
#plot_grid(alpha_plot, beta_plot, theta_plot, soc_raw_plot, soc_reduction_plot, soc_plot,  labels = c("","","","","",""), ncol = 2, align = 'v')
param_recov_plot = plot_grid(alpha_plot, beta_plot, annealing_plot, NULL, theta_plot_95, theta_plot_50, soc_raw_plot, soc_reduction_plot2, labels = c("a","b","c","","d","e","f","g"), ncol = 2, align = 'v')
#param_recov_plot

ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs_UNC6_sReduc_annealing_sim/param_recov_plot.pdf", plot = param_recov_plot, dpi = 600, width = 12, height = 9)

# soc (sigma)
modelSimulation_UNC6_sReduc_annealing$taskDifficulty2 = 'Low Uncertainty'
modelSimulation_UNC6_sReduc_annealing$taskDifficulty2[which(modelSimulation_UNC6_sReduc_annealing$taskDifficulty=='50%')] = 'Moderate Uncertainty'
modelSimulation_UNC6_sReduc_annealing$taskDifficulty2[which(modelSimulation_UNC6_sReduc_annealing$taskDifficulty=='78%')] = 'High Uncertainty'
modelSimulation_UNC6_sReduc_annealing$taskDifficulty2 = factor(modelSimulation_UNC6_sReduc_annealing$taskDifficulty2, levels=c('Low Uncertainty','Moderate Uncertainty','High Uncertainty'))
true_sigma = ggplot(modelSimulation_UNC6_sReduc_annealing,aes(round,soc,colour=theta))+
	geom_line(aes(group=subject)) + ylim(c(0,1))+
	facet_grid(.~taskDifficulty2)+
	labs(x='Round',y='True \nsocial learning weight')+
	scale_color_distiller(name='Theta',palette='RdYlBu',direction=-1)+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
fit_sigma = ggplot(modelSimulation_UNC6_sReduc_annealing,aes(round,estimated_soc,colour=theta))+
    geom_line(aes(group=subject)) + ylim(c(0,1))+ facet_grid(.~taskDifficulty2)+
    scale_color_distiller(name='Theta',palette='RdYlBu',direction=-1)+
    labs(x='Round',y='Fitted \nsocial learning weight')+
    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
sigma_plot = plot_grid(true_sigma, fit_sigma, labels = c("a","b"), ncol = 1, align = 'v')
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs_UNC6_sReduc_annealing_sim/sigma_plot.pdf", plot = sigma_plot, dpi = 600, width = 9, height = 6)

# netBeta
true_netBeta = ggplot(modelSimulation_UNC6_sReduc_annealing,aes(round,netBeta,colour=annealing))+
	geom_line(aes(group=subject)) + #ylim(c(0,6))+
	facet_grid(.~taskDifficulty2)+
	labs(x='Round',y=expression(paste('True inverse temperature ',beta,sep="")))+
	scale_color_distiller(name=expression(paste('',epsilon,sep="")),palette='RdYlBu',direction=-1)+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
fit_netBeta = ggplot(modelSimulation_UNC6_sReduc_annealing,aes(round,estimated_netBeta,colour=estimated_annealing))+
    geom_line(aes(group=subject)) + #ylim(c(0,6))+
    facet_grid(.~taskDifficulty2)+
    scale_color_distiller(name=expression(paste('Fitted ',epsilon,sep="")),palette='RdYlBu',direction=-1)+
    labs(x='Round',y=expression(paste('Fitted inverse temperature ',beta,sep="")))+
    panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE)+
	myTheme_legend()
netBeta_plot = plot_grid(true_netBeta, fit_netBeta, labels = c("a","b"), ncol = 1, align = 'v')
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs_UNC6_sReduc_annealing_sim/netBeta_plot.pdf", plot = netBeta_plot, dpi = 600, width = 9, height = 6)

## Correlation coefficients
cor.test(modelSimulation_UNC6_sReduc_annealing$alpha, modelSimulation_UNC6_sReduc_annealing$estimated_alpha)
cor.test(modelSimulation_UNC6_sReduc_annealing$beta, modelSimulation_UNC6_sReduc_annealing$estimated_beta)
cor.test(modelSimulation_UNC6_sReduc_annealing$annealing, modelSimulation_UNC6_sReduc_annealing$estimated_annealing)
cor.test(modelSimulation_UNC6_sReduc_annealing$theta, modelSimulation_UNC6_sReduc_annealing$estimated_theta)
cor.test(modelSimulation_UNC6_sReduc_annealing$soc_raw, modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw)
cor.test(modelSimulation_UNC6_sReduc_annealing$soc_reduction, modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction)

## alpha
alpha_scatter = ggplot(subset(modelSimulation_UNC6_sReduc_annealing, round==1),aes(alpha,estimated_alpha,colour=((alpha - estimated_alpha)^2)^(1/2)))+
    geom_point() + ylim(c(0,1))+ xlim(c(0,1))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    annotate("text", x=0.15, y=0.95, label='italic(r) == 0.80', parse = TRUE, size = 5)+
    labs(x=expression(paste('True ',alpha,sep="")),y=expression(paste('Fitted ',alpha,sep="")))+
    myTheme_legend()
## beta
beta_scatter = ggplot(subset(modelSimulation_UNC6_sReduc_annealing, round==1),aes(beta,estimated_beta,colour=((beta - estimated_beta)^2)^(1/2)))+
    geom_point() + #ylim(c(0,5)) + xlim(c(0,5))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',beta[0],sep="")),y=expression(paste('Fitted ',beta[0],sep="")))+
    annotate("text", x=-2, y=4, label='italic(r) == 0.69', parse = TRUE, size = 5)+
    myTheme_legend()
## annealing
annealing_scatter = ggplot(subset(modelSimulation_UNC6_sReduc_annealing, round==1),aes(annealing,estimated_annealing,colour=((annealing - estimated_annealing)^2)^(1/2)))+
    geom_point() + #ylim(c(0,5)) + xlim(c(0,5))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',epsilon,sep="")),y=expression(paste('Fitted ',epsilon,sep="")))+
    annotate("text", x=-3, y=6, label='italic(r) == 0.65', parse = TRUE, size = 5)+
    myTheme_legend()
## theta
theta_scatter = ggplot(subset(modelSimulation_UNC6_sReduc_annealing, round==1),aes(theta,estimated_theta,colour=((theta - estimated_theta)^2)^(1/2)))+
    geom_point() + #ylim(c(0,5)) + xlim(c(0,5))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',theta,sep="")),y=expression(paste('Fitted ',theta,sep="")))+
    annotate("text", x=-5, y=7.5, label='italic(r) == 0.54', parse = TRUE, size = 5)+
    myTheme_legend()
## soc_raw
soc_raw_scatter = ggplot(subset(modelSimulation_UNC6_sReduc_annealing, round==1),aes(soc_raw,estimated_soc_raw,colour=((soc_raw - estimated_soc_raw)^2)^(1/2)))+
    geom_point() + #ylim(c(-3,3)) + xlim(c(-3,3))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',sigma[0],sep="")),y=expression(paste('Fitted ',sigma[0],sep="")))+
    annotate("text", x=-5, y=2, label='italic(r) == 0.65', parse = TRUE, size = 5)+
    myTheme_legend()
## soc_reduction
soc_reduction_scatter = ggplot(subset(modelSimulation_UNC6_sReduc_annealing, round==1),aes(soc_reduc,estimated_soc_reduction,colour=((soc_reduc - estimated_soc_reduction)^2)^(1/2)))+
    geom_point() + #ylim(c(-6,0)) + xlim(c(-6,0))+
    scale_color_distiller(name='Difference',palette='RdYlBu',direction=-1)+
    labs(x=expression(paste('True ',delta,sep="")),y=expression(paste('Fitted ',delta,sep="")))+
    annotate("text", x=-10, y=5, label='italic(r) == 0.65', parse = TRUE, size = 5)+
    myTheme_legend()


param_recov_scatter_plot = plot_grid(alpha_scatter, beta_scatter, annealing_scatter, theta_scatter, soc_raw_scatter, soc_reduction_scatter, labels = c("a","b","c","d","e","f"), ncol = 2, align = 'v')
#param_recov_scatter_plot
ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs_UNC6_sReduc_annealing_sim/param_recov_scatter_plot.pdf", plot = param_recov_scatter_plot, dpi = 600, width = 9, height = 6)


### MCMC diagnostics
for(i in 1:nrow(temp_UNC6_sReduc_annealing_simulation)){
	temp_UNC6_sReduc_annealing_simulation$log_lik_n_eff[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('log_lik[',i,']',sep=''),'n_eff']
	temp_UNC6_sReduc_annealing_simulation$log_lik_Rhat[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('log_lik[',i,']',sep=''),'Rhat']
	temp_UNC6_sReduc_annealing_simulation$alpha_n_eff[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('alpha[',i,']',sep=''),'n_eff']
	temp_UNC6_sReduc_annealing_simulation$alpha_Rhat[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('alpha[',i,']',sep=''),'Rhat']
	temp_UNC6_sReduc_annealing_simulation$beta_n_eff[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('beta[',i,']',sep=''),'n_eff']
	temp_UNC6_sReduc_annealing_simulation$beta_Rhat[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('beta[',i,']',sep=''),'Rhat']
	temp_UNC6_sReduc_annealing_simulation$soc_raw_n_eff[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('soc_raw[',i,']',sep=''),'n_eff']
	temp_UNC6_sReduc_annealing_simulation$soc_raw_Rhat[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('soc_raw[',i,']',sep=''),'Rhat']
	temp_UNC6_sReduc_annealing_simulation$soc_reduction_n_eff[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('soc_reduction[',i,']',sep=''),'n_eff']
	temp_UNC6_sReduc_annealing_simulation$soc_reduction_Rhat[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('soc_reduction[',i,']',sep=''),'Rhat']
	temp_UNC6_sReduc_annealing_simulation$theta_n_eff[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('theta[',i,']',sep=''),'n_eff']
	temp_UNC6_sReduc_annealing_simulation$theta_Rhat[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('theta[',i,']',sep=''),'Rhat']
	temp_UNC6_sReduc_annealing_simulation$annealing_n_eff[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('annealing[',i,']',sep=''),'n_eff']
	temp_UNC6_sReduc_annealing_simulation$annealing_Rhat[i] <- summary(fit_UNC6_sReduc_annealing_simulation)$summary[paste('annealing[',i,']',sep=''),'Rhat']
}
### Save the file again
write.csv(temp_UNC6_sReduc_annealing_simulation,
			"~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/temp_UNC6_sReduc_annealing_simulation.csv",
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
	geom_errorbar(aes(x=subject, ymin=estimated_theta_25, ymax=estimated_theta_75, colour=estimated_frequency_dependence2), width=.05, alpha=1/1) +
	geom_point(aes(subject, estimated_theta, colour=estimated_frequency_dependence2),size=.5)+
	#geom_point(aes(subject, bestfit_theta), colour='orange',size=.5)+
	geom_point(aes(subject, true_theta),size=.5) +
	geom_hline(yintercept=0, colour='blue') +
	geom_hline(yintercept=1, colour='red') +
	myTheme()+
	scale_colour_manual(values=c("random-copying"="#999999",'neg-freq-dep'="blue", "weak-pos-freq-dep"="orange", "conformist"="red", "RANDOM"="#999999", "AL"="#999999"), name="Strength of \nfrequency dependence") +
	labs(x='Agents', y='theta')#+ facet_grid(.~taskDifficulty)
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
param_recov_plot_2575 = plot_grid(alpha_plot_50, beta_plot_50, annealing_plot_50, theta_plot_50, soc_raw_plot_50, soc_reduction_plot_50, labels = c("","","","","","",""), ncol = 2, align = 'v')
param_recov_plot_2575

ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs_UNC6_sReduc_annealing_sim/param_recov_plot_2575.pdf", plot = param_recov_plot_2575, dpi = 600, width = 9, height = 6)


## Frequency dependence categorisation accuracy
theta_recovery_data$is_correct_95 = as.numeric(theta_recovery_data$true_frequency_dependence==theta_recovery_data$estimated_frequency_dependence)
theta_recovery_data$is_correct_50 = as.numeric(theta_recovery_data$true_frequency_dependence==theta_recovery_data$estimated_frequency_dependence2)
table(theta_recovery_data$is_correct_95)
table(theta_recovery_data$is_correct_50)

theta_recovery_data$true_frequency_dependence4 = theta_recovery_data$true_frequency_dependence
theta_recovery_data$true_frequency_dependence4[which(theta_recovery_data$true_frequency_dependence4!='neg-freq-dep')] = 'pos-freq-dep'

theta_recovery_data$estimated_frequency_dependence_point = "random-copying"
theta_recovery_data$estimated_frequency_dependence_point[which(theta_recovery_data$estimated_theta >= 0.2)] = "pos-freq-dep"
theta_recovery_data$estimated_frequency_dependence_point[which(theta_recovery_data$estimated_theta <= -0.2)] = "neg-freq-dep"


## Goodness of fit - global parameters
parameters_true_fitted = data.frame(
	params = rep(c('mu_alpha','mu_beta','mu_annealing','mu_soc','mu_soc_reduction','mu_theta','s_alpha','s_beta','s_annealing','s_soc','s_soc_reduction','s_theta'),3),
	type = rep(c(rep('mu',6),rep('variation',6)),3),
	taskDifficulty = c(rep('Low',12),rep('Moderate',12),rep('High',12)),
	true_values = c(
		0.99,1.84,3.70,-1.55,-1.39,1.65,1.88,1.45,1.73,0.79,0.84,1.54,
		0.90,1.68,3.01,-2.37,-1.55,3.00,1.61,0.64,2.04,1.98,3.66,2.69,
		0.61,1.38,2.97,-2.16,-1.87,2.67,2.69,0.79,2.36,1.67,4.22,3.53
		)
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
	parameters_true_fitted$fit_p2.5[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_UNC6_sReduc_annealing_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,1]
	parameters_true_fitted$fit_p25[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_UNC6_sReduc_annealing_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,2]
	parameters_true_fitted$fit_p50[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_UNC6_sReduc_annealing_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,3]
	parameters_true_fitted$fit_p75[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_UNC6_sReduc_annealing_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,4]
	parameters_true_fitted$fit_p97.5[which(parameters_true_fitted$params==p)] = t(apply( eval(parse(text=paste("ms_UNC6_sReduc_annealing_simulation$",p,"[,1:3]",sep=""))), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)))[,5]
}

parameters_true_fitted$param_article_names = rep(c('mu_alpha','mu_beta','mu_epsilon','mu_sigma','mu_delta','mu_theta','s_alpha','s_beta','s_epsilon','s_sigma','s_delta','s_theta'),3)
parameters_true_fitted$taskDifficulty = factor(parameters_true_fitted$taskDifficulty, levels=c('Low','Moderate','High'))


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

ggsave(file = "~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/stan_model6/model6_simulation/figs_UNC6_sReduc_annealing_sim/param_recov_global_plot.pdf", plot = param_recov_global_plot, dpi = 600, width = 6, height = 6)


## ?% of parameters were recovered?
modelSimulation_UNC6_sReduc_annealing$is_alpha_recov = 0
modelSimulation_UNC6_sReduc_annealing$is_beta_recov = 0
modelSimulation_UNC6_sReduc_annealing$is_annealing_recov = 0
modelSimulation_UNC6_sReduc_annealing$is_soc_raw_recov = 0
modelSimulation_UNC6_sReduc_annealing$is_soc_reduc_recov = 0
modelSimulation_UNC6_sReduc_annealing$is_theta_recov = 0

modelSimulation_UNC6_sReduc_annealing$is_alpha_recov[which(modelSimulation_UNC6_sReduc_annealing$alpha>modelSimulation_UNC6_sReduc_annealing$estimated_alpha_low & modelSimulation_UNC6_sReduc_annealing$alpha<modelSimulation_UNC6_sReduc_annealing$estimated_alpha_high)] <- 1
modelSimulation_UNC6_sReduc_annealing$is_beta_recov[which(modelSimulation_UNC6_sReduc_annealing$beta>modelSimulation_UNC6_sReduc_annealing$estimated_beta_low & modelSimulation_UNC6_sReduc_annealing$beta<modelSimulation_UNC6_sReduc_annealing$estimated_beta_high)] <- 1
modelSimulation_UNC6_sReduc_annealing$is_annealing_recov[which(modelSimulation_UNC6_sReduc_annealing$annealing>modelSimulation_UNC6_sReduc_annealing$estimated_annealing_low & modelSimulation_UNC6_sReduc_annealing$annealing<modelSimulation_UNC6_sReduc_annealing$estimated_annealing_high)] <- 1
modelSimulation_UNC6_sReduc_annealing$is_soc_raw_recov[which(modelSimulation_UNC6_sReduc_annealing$soc_raw>modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_low & modelSimulation_UNC6_sReduc_annealing$soc_raw<modelSimulation_UNC6_sReduc_annealing$estimated_soc_raw_high)] <- 1
modelSimulation_UNC6_sReduc_annealing$is_soc_reduc_recov[which(modelSimulation_UNC6_sReduc_annealing$soc_reduc>modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_low & modelSimulation_UNC6_sReduc_annealing$soc_reduc<modelSimulation_UNC6_sReduc_annealing$estimated_soc_reduction_high)] <- 1
modelSimulation_UNC6_sReduc_annealing$is_theta_recov[which(modelSimulation_UNC6_sReduc_annealing$theta>modelSimulation_UNC6_sReduc_annealing$estimated_theta_low & modelSimulation_UNC6_sReduc_annealing$theta<modelSimulation_UNC6_sReduc_annealing$estimated_theta_high)] <- 1

length(which(subset(modelSimulation_UNC6_sReduc_annealing, round==1)$is_alpha_recov==1))/572 * 100
length(which(subset(modelSimulation_UNC6_sReduc_annealing, round==1)$is_beta_recov==1))/572 * 100
length(which(subset(modelSimulation_UNC6_sReduc_annealing, round==1)$is_annealing_recov==1))/572 * 100
length(which(subset(modelSimulation_UNC6_sReduc_annealing, round==1)$is_soc_raw_recov==1))/572 * 100
length(which(subset(modelSimulation_UNC6_sReduc_annealing, round==1)$is_soc_reduc_recov==1))/572 * 100
length(which(subset(modelSimulation_UNC6_sReduc_annealing, round==1)$is_theta_recov==1))/572 * 100
