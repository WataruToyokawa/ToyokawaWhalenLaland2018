# Individual Based Simulation using fitted parameter sets
# Considering: AL6_annealing and UNC6_sReduc_ann
# Requirement: Because this simulation is based on the model fit, I have to have fitted parameter sets listed in 'performanceSummary4' data.frame
# see also learningModelAnalysis.R
#
# 11 August 2018
# Wataru Toyokawa

library(readr)
library(pforeach) # If you don't have one -> install.packages("devtools"); devtools::install_github("hoxo-m/pforeach")
performanceSummary4 <- read_csv("~/Dropbox/wataru/St Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/performanceSummary4.csv")
performanceSummary4$smallestWaicModel[which(performanceSummary4$frequency_dependence_4=='--')] = NA
performanceSummary4$smallestWaicModel[which(is.na(performanceSummary4$alpha))] = NA
# Setting for all conditions
diff_25 = 1.264
diff_50 = 0.742
diff_78 = 0.3
repetition = 10000 #5000
lifetime = 90 #70
changingPoint = 41
slotmachines = 3
basePayoffMean = 3.1
payoffSD = 0.55
boxDistribution1 = c(0, 0, 1)
boxDistribution2 = c(2, 0, 1)
payoffGenerate = function(m) {
	payoff = rnorm(length(m), mean=m, sd=payoffSD)
	payoff[which(payoff<0)] = 0
	return (payoff)
}
softmax = function(Q, beta) {
	return( exp(beta*Q)/sum(exp(beta*Q)) )
}
test = function(behav_values) {
	return( exp(beta*Q)/sum(exp(beta*Q)) )
}
multiplyBeta = function(Q, beta) {
	return( beta*Q )
}
divideVector = function(numerator, denominator) {
		return( numerator/denominator )
}
multiplyTheta = function(f, theta) {
	return( f^theta )
}
frequencyDependentCopy = function(F, lastChoice, theta) {
	totalN=length(lastChoice)
	f = rep(F, totalN)
	f[lastChoice+((1:totalN)-1)*3] = f[lastChoice+((1:totalN)-1)*3] - 1
	f_matrix = matrix(f, ncol=totalN)
	ftheta = apply(f_matrix, 1, multiplyTheta, theta=theta) %>% t()
	return ( apply(ftheta, 1, divideVector, denominator=apply(ftheta,2,sum)) %>% t() )
}

system.time({
	#(* -- settings _25% condition -- *)
	qualityDiff = diff_25
	boxMeans1 = rep(basePayoffMean, slotmachines) + (boxDistribution1 - 2)*qualityDiff
	boxMeans2 = rep(basePayoffMean, slotmachines) + (boxDistribution2 - 2)*qualityDiff
	# Excluding the cases of groupSize == 1 (i.e. single condition)
	groupSizes_25 = performanceSummary4 %>% dplyr::filter(taskDifficulty=="25%") %>% dplyr::select(groupSize) %>% dplyr::count(groupSize) %>% spread(groupSize,n) %>% names() %>% as.numeric()

	#(* -- Initial Settings -- *)
	posthoc_sim_result_25 = list()

	#(* -- Simulation -- *)
	for(g in 1:length(groupSizes_25)) {
		size = groupSizes_25[g]
		if(size==1) size =2 # make arrays' dimension three (if getting one subject only, "apply()" function doesn't work well)
		thisSizeAmazonIDs = performanceSummary4 %>% dplyr::filter(groupSize==groupSizes_25[g]&taskDifficulty=='25%'&environment=='e1+e2'&as.numeric(as.factor(smallestWaicModel))>0) %>% dplyr::pull(amazonID) %>% table() %>% names()
		posthoc_sim_result_25[[paste("n=",groupSizes_25[g])]] = pforeach(rep = 1:repetition,.cores=8,.export=ls(envir=parent.frame()),.combine=rbind)({
			# this group's amazonID
			amazonIDs = sample(thisSizeAmazonIDs, size, replace = TRUE)
			models = rep(NA, size)
			alpha = rep(NA, size)
			beta = rep(NA, size)
			theta = rep(NA, size)
			soc_raw = rep(NA, size)
			soc_change = rep(NA, size)
			annealing = rep(NA, size)
			counter = 1
			for(i in amazonIDs) {
				# learning strategies
				this_subject = performanceSummary4 %>% dplyr::filter(subjectID == i & environment=='e1+e2')
				models[counter] = this_subject %>% dplyr::pull(smallestWaicModel)
				# parameters
				if(models[counter]=='UNC6_sReduc_ann') {
					alpha[counter] = this_subject %>% dplyr::pull(alpha)
					beta[counter] = this_subject %>% dplyr::pull(beta)
					theta[counter] = this_subject %>% dplyr::pull(theta)
					soc_raw[counter] = this_subject %>% dplyr::pull(soc_raw)
					soc_change[counter] = this_subject %>% dplyr::pull(soc_change)
					annealing[counter] = this_subject %>% dplyr::pull(annealing)
				} else {
					alpha[counter] = this_subject %>% dplyr::pull(alpha)
					beta[counter] = this_subject %>% dplyr::pull(beta)
					theta[counter] = 0
					soc_raw[counter] = -100
					soc_change[counter] = 0
					annealing[counter] = this_subject %>% dplyr::pull(annealing)
				}
				counter = counter + 1
			}
			# Setting initial values
			choices = matrix(nrow=size, ncol=lifetime)
			payoffs = matrix(nrow=size, ncol=lifetime)
			performance = matrix(nrow=size, ncol=lifetime)
			isThisBestOption = matrix(nrow=size, ncol=lifetime)
			optimalChoiceProb = matrix(nrow=size, ncol=lifetime)
			expectedPerformance = matrix(nrow=size, ncol=lifetime)
			softmaxChoiceProb = array(dim = c(slotmachines, lifetime, size))
			copyingChoiceProb = array(dim = c(slotmachines, lifetime, size))
			Q = array(dim = c(slotmachines, lifetime, size))
			Q[,1,] = 1e-10
			netChoiceProb = array(dim = c(slotmachines, lifetime, size))
			netChoiceProb[,1,] = 1/3
			socialFrequency = matrix(nrow=slotmachines, ncol=lifetime)
			socialFrequency[,] = 1e-1
			for(t in 1:lifetime) {
				# choices are done by each agent
				choices[,t] = mapply(function(p1,p2,p3){ sample(1:slotmachines, 1, prob=c(p1,p2,p3), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,], netChoiceProb[3,t,] )
				# this subject earns some money and update the performance
				if(t < changingPoint) {
					payoffs[,t] = payoffGenerate(boxMeans1[choices[,t]])
					performance[,t] = boxDistribution1[choices[,t]]
					isThisBestOption[,t] = 0
					isThisBestOption[,t][which(performance[,t]==1)] = 1
				}else{
					payoffs[,t] = payoffGenerate(boxMeans2[choices[,t]])
					performance[,t] = boxDistribution2[choices[,t]]
					isThisBestOption[,t] = 0
					isThisBestOption[,t][which(performance[,t]==2)] = 1
				}
				if(t < lifetime) {
					# update Q value based on Rescorla-Wagner model
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(slotmachines*size, lifetime)
					QQ[(choices[,t] + slotmachines*(1:size-1)),t+1] = QQ[(choices[,t] + slotmachines*(1:size-1)),t] + alpha * (payoffs[,t] - QQ[(choices[,t] + slotmachines*(1:size-1)),t])
					dim(QQ) = c(slotmachines, size, lifetime)
					Q = aperm(QQ, c(1,3,2))
					# update socialFrequency
					if(length(which(names(table(choices[,t]))==1))>0) {
						socialFrequency[1,t+1] = socialFrequency[1,t+1] + table(choices[,t])[which(names(table(choices[,t]))==1)][1]
					}
					if(length(which(names(table(choices[,t]))==2))>0) {
						socialFrequency[2,t+1] = socialFrequency[2,t+1] + table(choices[,t])[which(names(table(choices[,t]))==2)][1]
					}
					if(length(which(names(table(choices[,t]))==3))>0) {
						socialFrequency[3,t+1] = socialFrequency[3,t+1] + table(choices[,t])[which(names(table(choices[,t]))==3)][1]
					}
					# Calculating the next choice probability
					# It depends on what strategy each agent deploys
					# basic asocial softmax choice function
					Q_exp = apply(Q[,t+1,], 1, multiplyBeta, beta=(beta+annealing*(t+1)/lifetime)) %>% t() %>% apply(2,exp)
					softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
					freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], theta)
					# Calculating copying rate at this round
					soc = 1/(1+exp(-(soc_raw + soc_change*(t+1)/lifetime)))
					if(length(which(models!="UNC6_sReduc_ann"))>0){
						# Asocial learners have zero copying rate
						soc[which(models=="AL_ann")] = 0
					}
					# (1-soc) * softmax + soc * frequency_dependent_copying
					netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-soc)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=soc) %>% t()
					netChoiceProbAperm = aperm(netChoiceProb, c(1,3,2))
					dim(netChoiceProbAperm) = c(slotmachines*size, lifetime)
					dim(netMatrix) = c(slotmachines*size, 1)
					netChoiceProbAperm[,t+1] = netMatrix
					dim(netChoiceProbAperm) = c(slotmachines, size, lifetime)
					netChoiceProb = aperm(netChoiceProbAperm, c(1,3,2))
				}
			}
			# Calculating optimal choice probability for each agent
			for(i in 1:size){
				optimalChoiceProb[i,] = c(netChoiceProb[3,1:(changingPoint-1),i], netChoiceProb[1,changingPoint:lifetime,i])
				expectedPerformance[i,] = c(
					(netChoiceProb[1,1:(changingPoint-1),i]+netChoiceProb[2,1:(changingPoint-1),i])+2*netChoiceProb[3,1:(changingPoint-1),i],
					3*netChoiceProb[1,changingPoint:lifetime,i]+ netChoiceProb[2,changingPoint:lifetime,i]+ 2*netChoiceProb[3,changingPoint:lifetime,i]
				)
			}
			# submit this repetition's output
			print(
				c(apply(optimalChoiceProb,2,mean)
					,apply(expectedPerformance,2,mean)
					,apply(payoffs,2,mean)
					,apply(performance,2,mean)
					,apply(isThisBestOption,2,mean)
				)
			)
		})
		gc();gc() # rubbish collection
	}
})


posthoc_sim_summary_25 =
	apply(posthoc_sim_result_25[[paste("n=",1)]], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame() %>%
	cbind(groupSize = rep(1,lifetime*5)) %>%
	cbind(round = rep(1:lifetime, 5)) %>%
	cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime)))

for(g in groupSizes_25) {
	if(g != 1) {
		posthoc_sim_summary_25 = posthoc_sim_summary_25 %>%
		rbind(
			apply(posthoc_sim_result_25[[paste("n=",g)]], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame()%>%
			cbind(groupSize = rep(g,lifetime*5)) %>%
			cbind(round = rep(1:lifetime, 5)) %>%
			cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime)))
		)
	}
}

posthoc_sim_summary_25$sizeCategory = posthoc_sim_summary_25$groupSize
posthoc_sim_summary_25$indiv_small_large = 'Small group'
posthoc_sim_summary_25$indiv_small_large[which(posthoc_sim_summary_25$groupSize>9)] = 'Large group'
posthoc_sim_summary_25$indiv_small_large[which(posthoc_sim_summary_25$groupSize==1)] = 'Individual'

write.csv(posthoc_sim_summary_25,
			"~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/posthoc_sim_summary_25.csv",
			row.names=FALSE)



system.time({
	#(* -- settings _50% condition -- *)
	qualityDiff = diff_50
	boxMeans1 = rep(basePayoffMean, slotmachines) + (boxDistribution1 - 2)*qualityDiff
	boxMeans2 = rep(basePayoffMean, slotmachines) + (boxDistribution2 - 2)*qualityDiff
	# Excluding the cases of groupSize == 1 (i.e. single condition)
	groupSizes_50 = performanceSummary4 %>% dplyr::filter(taskDifficulty=="50%") %>% dplyr::select(groupSize) %>% dplyr::count(groupSize) %>% spread(groupSize,n) %>% names() %>% as.numeric()

	#(* -- Initial Settings -- *)
	posthoc_sim_result_50 = list()

	#(* -- Simulation -- *)
	for(g in 1:length(groupSizes_50)) {
		size = groupSizes_50[g]
		if(size==1) size =2 # make arrays' dimension three (if getting one subject only, "apply()" function doesn't work well)
		thisSizeAmazonIDs = performanceSummary4 %>% dplyr::filter(groupSize==groupSizes_50[g]&taskDifficulty=='50%'&environment=='e1+e2'&as.numeric(as.factor(smallestWaicModel))>0) %>% dplyr::pull(amazonID) %>% table() %>% names()
		posthoc_sim_result_50[[paste("n=",groupSizes_50[g])]] = pforeach(rep = 1:repetition,.cores=8,.export=ls(envir=parent.frame()),.combine=rbind)({
			# this group's amazonID
			amazonIDs = sample(thisSizeAmazonIDs, size, replace = TRUE)
			models = rep(NA, size)
			alpha = rep(NA, size)
			beta = rep(NA, size)
			theta = rep(NA, size)
			soc_raw = rep(NA, size)
			soc_change = rep(NA, size)
			annealing = rep(NA, size)
			counter = 1
			for(i in amazonIDs) {
				# learning strategies
				this_subject = performanceSummary4 %>% dplyr::filter(subjectID == i & environment=='e1+e2')
				models[counter] = this_subject %>% dplyr::pull(smallestWaicModel)
				# parameters
				if(models[counter]=='UNC6_sReduc_ann') {
					alpha[counter] = this_subject %>% dplyr::pull(alpha)
					beta[counter] = this_subject %>% dplyr::pull(beta)
					theta[counter] = this_subject %>% dplyr::pull(theta)
					soc_raw[counter] = this_subject %>% dplyr::pull(soc_raw)
					soc_change[counter] = this_subject %>% dplyr::pull(soc_change)
					annealing[counter] = this_subject %>% dplyr::pull(annealing)
				} else {
					alpha[counter] = this_subject %>% dplyr::pull(alpha)
					beta[counter] = this_subject %>% dplyr::pull(beta)
					theta[counter] = 0
					soc_raw[counter] = -100
					soc_change[counter] = 0
					annealing[counter] = this_subject %>% dplyr::pull(annealing)
				}
				counter = counter + 1
			}
			# Setting initial values
			choices = matrix(nrow=size, ncol=lifetime)
			payoffs = matrix(nrow=size, ncol=lifetime)
			performance = matrix(nrow=size, ncol=lifetime)
			isThisBestOption = matrix(nrow=size, ncol=lifetime)
			optimalChoiceProb = matrix(nrow=size, ncol=lifetime)
			expectedPerformance = matrix(nrow=size, ncol=lifetime)
			softmaxChoiceProb = array(dim = c(slotmachines, lifetime, size))
			copyingChoiceProb = array(dim = c(slotmachines, lifetime, size))
			Q = array(dim = c(slotmachines, lifetime, size))
			Q[,1,] = 1e-10
			netChoiceProb = array(dim = c(slotmachines, lifetime, size))
			netChoiceProb[,1,] = 1/3
			socialFrequency = matrix(nrow=slotmachines, ncol=lifetime)
			socialFrequency[,] = 1e-1
			for(t in 1:lifetime) {
				# choices are done by each agent
				choices[,t] = mapply(function(p1,p2,p3){ sample(1:slotmachines, 1, prob=c(p1,p2,p3), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,], netChoiceProb[3,t,] )
				# this subject earns some money and update the performance
				if(t < changingPoint) {
					payoffs[,t] = payoffGenerate(boxMeans1[choices[,t]])
					performance[,t] = boxDistribution1[choices[,t]]
					isThisBestOption[,t] = 0
					isThisBestOption[,t][which(performance[,t]==1)] = 1
				}else{
					payoffs[,t] = payoffGenerate(boxMeans2[choices[,t]])
					performance[,t] = boxDistribution2[choices[,t]]
					isThisBestOption[,t] = 0
					isThisBestOption[,t][which(performance[,t]==2)] = 1
				}
				if(t < lifetime) {
					# update Q value based on Rescorla-Wagner model
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(slotmachines*size, lifetime)
					QQ[(choices[,t] + slotmachines*(1:size-1)),t+1] = QQ[(choices[,t] + slotmachines*(1:size-1)),t] + alpha * (payoffs[,t] - QQ[(choices[,t] + slotmachines*(1:size-1)),t])
					dim(QQ) = c(slotmachines, size, lifetime)
					Q = aperm(QQ, c(1,3,2))
					# update socialFrequency
					if(length(which(names(table(choices[,t]))==1))>0) {
						socialFrequency[1,t+1] = socialFrequency[1,t+1] + table(choices[,t])[which(names(table(choices[,t]))==1)][1]
					}
					if(length(which(names(table(choices[,t]))==2))>0) {
						socialFrequency[2,t+1] = socialFrequency[2,t+1] + table(choices[,t])[which(names(table(choices[,t]))==2)][1]
					}
					if(length(which(names(table(choices[,t]))==3))>0) {
						socialFrequency[3,t+1] = socialFrequency[3,t+1] + table(choices[,t])[which(names(table(choices[,t]))==3)][1]
					}
					# Calculating the next choice probability
					# It depends on what strategy each agent deploys
					# basic asocial softmax choice function
					Q_exp = apply(Q[,t+1,], 1, multiplyBeta, beta=(beta+annealing*(t+1)/lifetime)) %>% t() %>% apply(2,exp)
					softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
					freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], theta)
					# Calculating copying rate at this round
					soc = 1/(1+exp(-(soc_raw + soc_change*(t+1)/lifetime)))
					if(length(which(models!="UNC6_sReduc_ann"))>0){
						# Asocial learners have zero copying rate
						soc[which(models=="AL_ann")] = 0
					}
					# (1-soc) * softmax + soc * frequency_dependent_copying
					netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-soc)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=soc) %>% t()
					netChoiceProbAperm = aperm(netChoiceProb, c(1,3,2))
					dim(netChoiceProbAperm) = c(slotmachines*size, lifetime)
					dim(netMatrix) = c(slotmachines*size, 1)
					netChoiceProbAperm[,t+1] = netMatrix
					dim(netChoiceProbAperm) = c(slotmachines, size, lifetime)
					netChoiceProb = aperm(netChoiceProbAperm, c(1,3,2))
				}
			}
			# Calculating optimal choice probability for each agent
			for(i in 1:size){
				optimalChoiceProb[i,] = c(netChoiceProb[3,1:(changingPoint-1),i], netChoiceProb[1,changingPoint:lifetime,i])
				expectedPerformance[i,] = c(
					(netChoiceProb[1,1:(changingPoint-1),i]+netChoiceProb[2,1:(changingPoint-1),i])+2*netChoiceProb[3,1:(changingPoint-1),i],
					3*netChoiceProb[1,changingPoint:lifetime,i]+ netChoiceProb[2,changingPoint:lifetime,i]+ 2*netChoiceProb[3,changingPoint:lifetime,i]
				)
			}
			# submit this repetition's output
			print(
				c(apply(optimalChoiceProb,2,mean)
					,apply(expectedPerformance,2,mean)
					,apply(payoffs,2,mean)
					,apply(performance,2,mean)
					,apply(isThisBestOption,2,mean)
				)
			)
		})
		gc();gc() # rubbish collection
	}
})


posthoc_sim_summary_50 =
	apply(posthoc_sim_result_50[[paste("n=",1)]], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame() %>%
	cbind(groupSize = rep(1,lifetime*5)) %>%
	cbind(round = rep(1:lifetime, 5)) %>%
	cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime)))

for(g in groupSizes_50) {
	if(g != 1) {
		posthoc_sim_summary_50 = posthoc_sim_summary_50 %>%
		rbind(
			apply(posthoc_sim_result_50[[paste("n=",g)]], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame()%>%
			cbind(groupSize = rep(g,lifetime*5)) %>%
			cbind(round = rep(1:lifetime, 5)) %>%
			cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime)))
		)
	}
}

posthoc_sim_summary_50$sizeCategory = posthoc_sim_summary_50$groupSize
posthoc_sim_summary_50$indiv_small_large = 'Small group'
posthoc_sim_summary_50$indiv_small_large[which(posthoc_sim_summary_50$groupSize>9)] = 'Large group'
posthoc_sim_summary_50$indiv_small_large[which(posthoc_sim_summary_50$groupSize==1)] = 'Individual'

write.csv(posthoc_sim_summary_50,
			"~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/posthoc_sim_summary_50.csv",
			row.names=FALSE)



system.time({
	#(* -- settings _78% condition -- *)
	qualityDiff = diff_78
	boxMeans1 = rep(basePayoffMean, slotmachines) + (boxDistribution1 - 2)*qualityDiff
	boxMeans2 = rep(basePayoffMean, slotmachines) + (boxDistribution2 - 2)*qualityDiff
	# Excluding the cases of groupSize == 1 (i.e. single condition)
	groupSizes_78 = performanceSummary4 %>% dplyr::filter(taskDifficulty=="78%") %>% dplyr::select(groupSize) %>% dplyr::count(groupSize) %>% spread(groupSize,n) %>% names() %>% as.numeric()

	#(* -- Initial Settings -- *)
	posthoc_sim_result_78 = list()

	#(* -- Simulation -- *)
	for(g in 1:length(groupSizes_78)) {
		size = groupSizes_78[g]
		if(size==1) size =2 # make arrays' dimension three (if getting one subject only, "apply()" function doesn't work well)
		thisSizeAmazonIDs = performanceSummary4 %>% dplyr::filter(groupSize==groupSizes_78[g]&taskDifficulty=='78%'&environment=='e1+e2'&as.numeric(as.factor(smallestWaicModel))>0) %>% dplyr::pull(amazonID) %>% table() %>% names()
		posthoc_sim_result_78[[paste("n=",groupSizes_78[g])]] = pforeach(rep = 1:repetition,.cores=8,.export=ls(envir=parent.frame()),.combine=rbind)({
			# this group's amazonID
			amazonIDs = sample(thisSizeAmazonIDs, size, replace = TRUE)
			models = rep(NA, size)
			alpha = rep(NA, size)
			beta = rep(NA, size)
			theta = rep(NA, size)
			soc_raw = rep(NA, size)
			soc_change = rep(NA, size)
			annealing = rep(NA, size)
			counter = 1
			for(i in amazonIDs) {
				# learning strategies
				this_subject = performanceSummary4 %>% dplyr::filter(subjectID == i & environment=='e1+e2')
				models[counter] = this_subject %>% dplyr::pull(smallestWaicModel)
				# parameters
				if(models[counter]=='UNC6_sReduc_ann') {
					alpha[counter] = this_subject %>% dplyr::pull(alpha)
					beta[counter] = this_subject %>% dplyr::pull(beta)
					theta[counter] = this_subject %>% dplyr::pull(theta)
					soc_raw[counter] = this_subject %>% dplyr::pull(soc_raw)
					soc_change[counter] = this_subject %>% dplyr::pull(soc_change)
					annealing[counter] = this_subject %>% dplyr::pull(annealing)
				} else {
					alpha[counter] = this_subject %>% dplyr::pull(alpha)
					beta[counter] = this_subject %>% dplyr::pull(beta)
					theta[counter] = 0
					soc_raw[counter] = -100
					soc_change[counter] = 0
					annealing[counter] = this_subject %>% dplyr::pull(annealing)
				}
				counter = counter + 1
			}
			# Setting initial values
			choices = matrix(nrow=size, ncol=lifetime)
			payoffs = matrix(nrow=size, ncol=lifetime)
			performance = matrix(nrow=size, ncol=lifetime)
			isThisBestOption = matrix(nrow=size, ncol=lifetime)
			optimalChoiceProb = matrix(nrow=size, ncol=lifetime)
			expectedPerformance = matrix(nrow=size, ncol=lifetime)
			softmaxChoiceProb = array(dim = c(slotmachines, lifetime, size))
			copyingChoiceProb = array(dim = c(slotmachines, lifetime, size))
			Q = array(dim = c(slotmachines, lifetime, size))
			Q[,1,] = 1e-10
			netChoiceProb = array(dim = c(slotmachines, lifetime, size))
			netChoiceProb[,1,] = 1/3
			socialFrequency = matrix(nrow=slotmachines, ncol=lifetime)
			socialFrequency[,] = 1e-1
			for(t in 1:lifetime) {
				# choices are done by each agent
				choices[,t] = mapply(function(p1,p2,p3){ sample(1:slotmachines, 1, prob=c(p1,p2,p3), replace=FALSE) }, netChoiceProb[1,t,], netChoiceProb[2,t,], netChoiceProb[3,t,] )
				# this subject earns some money and update the performance
				if(t < changingPoint) {
					payoffs[,t] = payoffGenerate(boxMeans1[choices[,t]])
					performance[,t] = boxDistribution1[choices[,t]]
					isThisBestOption[,t] = 0
					isThisBestOption[,t][which(performance[,t]==1)] = 1
				}else{
					payoffs[,t] = payoffGenerate(boxMeans2[choices[,t]])
					performance[,t] = boxDistribution2[choices[,t]]
					isThisBestOption[,t] = 0
					isThisBestOption[,t][which(performance[,t]==2)] = 1
				}
				if(t < lifetime) {
					# update Q value based on Rescorla-Wagner model
					Q[,t+1,] = Q[,t,]
					QQ = aperm(Q, c(1,3,2))
					dim(QQ) = c(slotmachines*size, lifetime)
					QQ[(choices[,t] + slotmachines*(1:size-1)),t+1] = QQ[(choices[,t] + slotmachines*(1:size-1)),t] + alpha * (payoffs[,t] - QQ[(choices[,t] + slotmachines*(1:size-1)),t])
					dim(QQ) = c(slotmachines, size, lifetime)
					Q = aperm(QQ, c(1,3,2))
					# update socialFrequency
					if(length(which(names(table(choices[,t]))==1))>0) {
						socialFrequency[1,t+1] = socialFrequency[1,t+1] + table(choices[,t])[which(names(table(choices[,t]))==1)][1]
					}
					if(length(which(names(table(choices[,t]))==2))>0) {
						socialFrequency[2,t+1] = socialFrequency[2,t+1] + table(choices[,t])[which(names(table(choices[,t]))==2)][1]
					}
					if(length(which(names(table(choices[,t]))==3))>0) {
						socialFrequency[3,t+1] = socialFrequency[3,t+1] + table(choices[,t])[which(names(table(choices[,t]))==3)][1]
					}
					# Calculating the next choice probability
					# It depends on what strategy each agent deploys
					# basic asocial softmax choice function
					Q_exp = apply(Q[,t+1,], 1, multiplyBeta, beta=(beta+annealing*(t+1)/lifetime)) %>% t() %>% apply(2,exp)
					softmaxMatrix = Q_exp %>% apply(1, divideVector, denominator = apply(Q_exp,2,sum)) %>% t()
					freqDepenMatrix = frequencyDependentCopy(socialFrequency[,t+1], choices[,t], theta)
					# Calculating copying rate at this round
					soc = 1/(1+exp(-(soc_raw + soc_change*(t+1)/lifetime)))
					if(length(which(models!="UNC6_sReduc_ann"))>0){
						# Asocial learners have zero copying rate
						soc[which(models=="AL_ann")] = 0
					}
					# (1-soc) * softmax + soc * frequency_dependent_copying
					netMatrix = apply(softmaxMatrix, 1, multiplyBeta, beta=(1-soc)) %>% t() + apply(freqDepenMatrix, 1, multiplyBeta, beta=soc) %>% t()
					netChoiceProbAperm = aperm(netChoiceProb, c(1,3,2))
					dim(netChoiceProbAperm) = c(slotmachines*size, lifetime)
					dim(netMatrix) = c(slotmachines*size, 1)
					netChoiceProbAperm[,t+1] = netMatrix
					dim(netChoiceProbAperm) = c(slotmachines, size, lifetime)
					netChoiceProb = aperm(netChoiceProbAperm, c(1,3,2))
				}
			}
			# Calculating optimal choice probability for each agent
			for(i in 1:size){
				optimalChoiceProb[i,] = c(netChoiceProb[3,1:(changingPoint-1),i], netChoiceProb[1,changingPoint:lifetime,i])
				expectedPerformance[i,] = c(
					(netChoiceProb[1,1:(changingPoint-1),i]+netChoiceProb[2,1:(changingPoint-1),i])+2*netChoiceProb[3,1:(changingPoint-1),i],
					3*netChoiceProb[1,changingPoint:lifetime,i]+ netChoiceProb[2,changingPoint:lifetime,i]+ 2*netChoiceProb[3,changingPoint:lifetime,i]
				)
			}
			# submit this repetition's output
			print(
				c(apply(optimalChoiceProb,2,mean)
					,apply(expectedPerformance,2,mean)
					,apply(payoffs,2,mean)
					,apply(performance,2,mean)
					,apply(isThisBestOption,2,mean)
				)
			)
		})
		gc();gc() # rubbish collection
	}
})


posthoc_sim_summary_78 =
	apply(posthoc_sim_result_78[[paste("n=",1)]], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame() %>%
	cbind(groupSize = rep(1,lifetime*5)) %>%
	cbind(round = rep(1:lifetime, 5)) %>%
	cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime)))

for(g in groupSizes_78) {
	if(g != 1) {
		posthoc_sim_summary_78 = posthoc_sim_summary_78 %>%
		rbind(
			apply(posthoc_sim_result_78[[paste("n=",g)]], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame()%>%
			cbind(groupSize = rep(g,lifetime*5)) %>%
			cbind(round = rep(1:lifetime, 5)) %>%
			cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime)))
		)
	}
}

posthoc_sim_summary_78$sizeCategory = posthoc_sim_summary_78$groupSize
posthoc_sim_summary_78$indiv_small_large = 'Small group'
posthoc_sim_summary_78$indiv_small_large[which(posthoc_sim_summary_78$groupSize>9)] = 'Large group'
posthoc_sim_summary_78$indiv_small_large[which(posthoc_sim_summary_78$groupSize==1)] = 'Individual'

write.csv(posthoc_sim_summary_78,
			"~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/posthoc_sim_summary_78.csv",
			row.names=FALSE)




posthoc_sim_summary_25$taskDifficulty3 = 'Low Uncertainty'
posthoc_sim_summary_50$taskDifficulty3 = 'Moderate Uncertainty'
posthoc_sim_summary_78$taskDifficulty3 = 'High Uncertainty'
posthoc_sim_summary_all = rbind(posthoc_sim_summary_25, posthoc_sim_summary_50) %>% rbind(posthoc_sim_summary_78)
posthoc_sim_summary_all$taskDifficulty3 = factor(posthoc_sim_summary_all$taskDifficulty3, levels=c('Low Uncertainty', 'Moderate Uncertainty', 'High Uncertainty'))
posthoc_sim_summary_all$sizeCategory_string = apply(rbind(rep('n = ', nrow(posthoc_sim_summary_all)),as.numeric(as.character(posthoc_sim_summary_all$groupSize))), 2, str_c, collapse="")
posthoc_sim_summary_all$sizeCategory_string = factor(posthoc_sim_summary_all$sizeCategory_string, levels=c("n = 1","n = 2","n = 3","n = 4","n = 5","n = 6","n = 7","n = 8","n = 9","n = 10","n = 11","n = 12","n = 13","n = 14","n = 15","n = 16","n = 17","n = 18","n = 19","n = 20","n = 21","n = 22","n = 23","n = 24","n = 25","n = 26","n = 27"))
allBehaviouralData2$sizeCategory_string = apply(rbind(rep('n = ', nrow(allBehaviouralData2)),allBehaviouralData2$sizeCategory), 2, str_c, collapse="")
allBehaviouralData2$sizeCategory_string = factor(allBehaviouralData2$sizeCategory_string, levels=c("n = 1","n = 2","n = 3","n = 4","n = 5","n = 6","n = 7","n = 8","n = 9","n = 10","n = 11","n = 12","n = 13","n = 14","n = 15","n = 16","n = 17","n = 18","n = 19","n = 20","n = 21","n = 22","n = 23","n = 24","n = 25","n = 26","n = 27"))

write.csv(posthoc_sim_summary_all,
			"~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/posthoc_sim_summary_all.csv",
			row.names=FALSE)


## Figure: Time evolution of the proportion of choosing the best option for each group
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
	axis.text.y = element_text(size=11, family="Times" ,colour='black'),
	strip.text.y = element_text(angle = 0))+
	NULL

ggsave(file = "~/Dropbox/wataru/papers/ConformityOnlineExp/data_and_analysis/FigureS_eachGroupsPerformance.pdf", plot = FigureS_eachGroupsPerformance, dpi = 600, width = 9, height = 9)


## Figure: A relation ship between group size and average proportion of choosing the optimal option.
### _25%
performanceSummary_posthocSimulation_25 = apply(posthoc_sim_result_25[[paste("n=",1)]][,1:(changingPoint-1)], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>%
	rbind(apply(posthoc_sim_result_25[[paste("n=",1)]][,changingPoint:70], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
	rbind(apply(posthoc_sim_result_25[[paste("n=",1)]][,71:lifetime], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
	data.frame() %>%
	cbind(environment = c('environment 1','environment 2','environment 3')) %>%
	cbind(groupSize = rep(1, 3)) %>%
	cbind(value = rep('optimalChoiceProb', 3))

for(g in groupSizes_25) {
	if(g != 1) {
		performanceSummary_posthocSimulation_25 = performanceSummary_posthocSimulation_25 %>%
		rbind(
			apply(posthoc_sim_result_25[[paste("n=",g)]][,1:(changingPoint-1)], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>%
			rbind(apply(posthoc_sim_result_25[[paste("n=",g)]][,changingPoint:70], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
			rbind(apply(posthoc_sim_result_25[[paste("n=",g)]][,71:lifetime], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
			data.frame() %>%
			cbind(environment = c('environment 1','environment 2','environment 3')) %>%
			cbind(groupSize = rep(g, 3)) %>%
			cbind(value = rep('optimalChoiceProb', 3))
		)
	}
}

performanceSummary_posthocSimulation_25$taskDifficulty3 = 'Low Uncertainty'

### _50%
performanceSummary_posthocSimulation_50 = apply(posthoc_sim_result_50[[paste("n=",1)]][,1:(changingPoint-1)], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>%
	rbind(apply(posthoc_sim_result_50[[paste("n=",1)]][,changingPoint:70], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
	rbind(apply(posthoc_sim_result_50[[paste("n=",1)]][,71:lifetime], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
	data.frame() %>%
	cbind(environment = c('environment 1','environment 2','environment 3')) %>%
	cbind(groupSize = rep(1, 3)) %>%
	cbind(value = rep('optimalChoiceProb', 3))

for(g in groupSizes_50) {
	if(g != 1) {
		performanceSummary_posthocSimulation_50 = performanceSummary_posthocSimulation_50 %>%
		rbind(
			apply(posthoc_sim_result_50[[paste("n=",g)]][,1:(changingPoint-1)], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>%
			rbind(apply(posthoc_sim_result_50[[paste("n=",g)]][,changingPoint:70], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
			rbind(apply(posthoc_sim_result_50[[paste("n=",g)]][,71:lifetime], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
			data.frame() %>%
			cbind(environment = c('environment 1','environment 2','environment 3')) %>%
			cbind(groupSize = rep(g, 3)) %>%
			cbind(value = rep('optimalChoiceProb', 3))
		)
	}
}

performanceSummary_posthocSimulation_50$taskDifficulty3 = 'Moderate Uncertainty'

### _78%
performanceSummary_posthocSimulation_78 = apply(posthoc_sim_result_78[[paste("n=",1)]][,1:(changingPoint-1)], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>%
	rbind(apply(posthoc_sim_result_78[[paste("n=",1)]][,changingPoint:70], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
	rbind(apply(posthoc_sim_result_78[[paste("n=",1)]][,71:lifetime], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
	data.frame() %>%
	cbind(environment = c('environment 1','environment 2','environment 3')) %>%
	cbind(groupSize = rep(1, 3)) %>%
	cbind(value = rep('optimalChoiceProb', 3))

for(g in groupSizes_78) {
	if(g != 1) {
		performanceSummary_posthocSimulation_78 = performanceSummary_posthocSimulation_78 %>%
		rbind(
			apply(posthoc_sim_result_78[[paste("n=",g)]][,1:(changingPoint-1)], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>%
			rbind(apply(posthoc_sim_result_78[[paste("n=",g)]][,changingPoint:70], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
			rbind(apply(posthoc_sim_result_78[[paste("n=",g)]][,71:lifetime], 1, mean) %>% quantile(probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t()) %>%
			data.frame() %>%
			cbind(environment = c('environment 1','environment 2','environment 3')) %>%
			cbind(groupSize = rep(g, 3)) %>%
			cbind(value = rep('optimalChoiceProb', 3))
		)
	}
}

performanceSummary_posthocSimulation_78$taskDifficulty3 = 'High Uncertainty'

performanceSummary_posthocSimulation = rbind(performanceSummary_posthocSimulation_25, performanceSummary_posthocSimulation_50) %>% rbind(performanceSummary_posthocSimulation_78) %>% data.frame()
performanceSummary_posthocSimulation$taskDifficulty3 = factor(performanceSummary_posthocSimulation$taskDifficulty3, levels=c('Low Uncertainty', 'Moderate Uncertainty', 'High Uncertainty'))
names(performanceSummary_posthocSimulation) = c('p2.5','p25','p50','p75','p97.5','environment','groupSize','value','taskDifficulty3')

## Figure: A relation ship between group size and average proportion of choosing the optimal option.
ggplot() +
	geom_point(data = performanceSummary_posthocSimulation %>% dplyr::filter(environment!='environment 3'), mapping=aes(x=groupSize+0.3,y=p50), colour='red')+
	geom_errorbar(data = performanceSummary_posthocSimulation %>% dplyr::filter(environment!='environment 3'), mapping=aes(x=groupSize+0.3,ymin=p2.5,ymax=p97.5),colour='grey80')+
	geom_errorbar(data = performanceSummary_posthocSimulation %>% dplyr::filter(environment!='environment 3'), mapping=aes(x=groupSize+0.3,ymin=p25,ymax=p75),colour='grey50')+
	#stat_summary(data=subset(performanceSummary4,environment!='e1+e2'), mapping=aes(groupSize, optimalChoiceRate, group=groupSize),fun.y = mean,geom="point",size=3) +
	facet_grid(environment ~ taskDifficulty3)+
	scale_y_continuous(limits = c(0,1), breaks=c(0.0,0.33,0.5,0.75,1.0)) +
	labs(x='Group size', y='Average proportion of choosing\n the best option')+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
	geom_hline(yintercept=1/3, linetype="dashed")+
	myTheme()+
	theme(legend.position='none',
		axis.text.y = element_text(size=13, family="Times" ,colour='black'),
		strip.text.y = element_text(angle = 0))+
	NULL


## individual-small-large comparison plot
posthoc_sim_summary_indivSmallLarge_25 =
	apply(posthoc_sim_result_25[[paste("n=",1)]], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame() %>%
	cbind(mean = apply(posthoc_sim_result_25[[paste("n=",1)]], 2, mean)) %>%
	cbind(groupSize = rep(1,lifetime*5)) %>%
	cbind(round = rep(1:lifetime, 5)) %>%
	cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime))) %>%
	cbind(indiv_small_large = rep('Individual', lifetime*5))

posthoc_sim_summary_indivSmallLarge_25 = posthoc_sim_summary_indivSmallLarge_25 %>%
	rbind(
		apply(rbind(posthoc_sim_result_25[[paste("n=",3)]],posthoc_sim_result_25[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",6)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",9)]]), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame()%>%
		cbind(mean = apply(rbind(posthoc_sim_result_25[[paste("n=",3)]],posthoc_sim_result_25[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",6)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",9)]]), 2, mean)) %>%
		cbind(groupSize = rep(g,lifetime*5)) %>%
		cbind(round = rep(1:lifetime, 5)) %>%
		cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime))) %>%
		cbind(indiv_small_large = rep('Small', lifetime*5))
	)
posthoc_sim_summary_indivSmallLarge_25 = posthoc_sim_summary_indivSmallLarge_25 %>%
	rbind(
		apply(rbind(posthoc_sim_result_25[[paste("n=",12)]],posthoc_sim_result_25[[paste("n=",15)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",16)]]), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame() %>%
		cbind(mean=apply(rbind(posthoc_sim_result_25[[paste("n=",12)]],posthoc_sim_result_25[[paste("n=",15)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",16)]]), 2, mean)) %>%
		cbind(groupSize = rep(g,lifetime*5)) %>%
		cbind(round = rep(1:lifetime, 5)) %>%
		cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime))) %>%
		cbind(indiv_small_large = rep('Large', lifetime*5))
	)

posthoc_sim_summary_indivSmallLarge_50 =
	apply(posthoc_sim_result_50[[paste("n=",1)]], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame() %>%
	cbind(mean = apply(posthoc_sim_result_50[[paste("n=",1)]], 2, mean)) %>%
	cbind(groupSize = rep(1,lifetime*5)) %>%
	cbind(round = rep(1:lifetime, 5)) %>%
	cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime))) %>%
	cbind(indiv_small_large = rep('Individual', lifetime*5))

posthoc_sim_summary_indivSmallLarge_50 = posthoc_sim_summary_indivSmallLarge_50 %>%
	rbind(
		apply(rbind(posthoc_sim_result_50[[paste("n=",2)]],posthoc_sim_result_50[[paste("n=",3)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",4)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",6)]]), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame()%>%
		cbind(mean = apply(rbind(posthoc_sim_result_50[[paste("n=",2)]],posthoc_sim_result_50[[paste("n=",3)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",4)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",6)]]), 2, mean)) %>%
		cbind(groupSize = rep(g,lifetime*5)) %>%
		cbind(round = rep(1:lifetime, 5)) %>%
		cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime))) %>%
		cbind(indiv_small_large = rep('Small', lifetime*5))
	)
posthoc_sim_summary_indivSmallLarge_50 = posthoc_sim_summary_indivSmallLarge_50 %>%
	rbind(
		apply(rbind(posthoc_sim_result_50[[paste("n=",8)]],posthoc_sim_result_50[[paste("n=",9)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",17)]]), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame() %>%
		cbind(mean = apply(rbind(posthoc_sim_result_50[[paste("n=",8)]],posthoc_sim_result_50[[paste("n=",9)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",17)]]), 2, mean)) %>%
		cbind(groupSize = rep(g,lifetime*5)) %>%
		cbind(round = rep(1:lifetime, 5)) %>%
		cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime))) %>%
		cbind(indiv_small_large = rep('Large', lifetime*5))
	)

posthoc_sim_summary_indivSmallLarge_78 =
	apply(posthoc_sim_result_78[[paste("n=",1)]], 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame() %>%
	cbind(mean = apply(posthoc_sim_result_78[[paste("n=",1)]], 2, mean)) %>%
	cbind(groupSize = rep(1,lifetime*5)) %>%
	cbind(round = rep(1:lifetime, 5)) %>%
	cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime))) %>%
	cbind(indiv_small_large = rep('Individual', lifetime*5))

posthoc_sim_summary_indivSmallLarge_78 = posthoc_sim_summary_indivSmallLarge_78 %>%
	rbind(
		apply(rbind(posthoc_sim_result_78[[paste("n=",2)]],posthoc_sim_result_78[[paste("n=",3)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",6)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",7)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",8)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",9)]]), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame()%>%
		cbind(mean = apply(rbind(posthoc_sim_result_78[[paste("n=",2)]],posthoc_sim_result_78[[paste("n=",3)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",6)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",7)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",8)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",9)]]), 2, mean)) %>%
		cbind(groupSize = rep(g,lifetime*5)) %>%
		cbind(round = rep(1:lifetime, 5)) %>%
		cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime))) %>%
		cbind(indiv_small_large = rep('Small', lifetime*5))
	)
posthoc_sim_summary_indivSmallLarge_78 = posthoc_sim_summary_indivSmallLarge_78 %>%
	rbind(
		apply(rbind(posthoc_sim_result_78[[paste("n=",12)]],posthoc_sim_result_78[[paste("n=",14)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",15)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",18)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",20)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",23)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",24)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",27)]]), 2, quantile, probs = c(0.025,0.25,0.5,0.75,0.975)) %>% t() %>% data.frame() %>%
		cbind(mean = apply(rbind(posthoc_sim_result_78[[paste("n=",12)]],posthoc_sim_result_78[[paste("n=",14)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",15)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",18)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",20)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",23)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",24)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",27)]]), 2, mean)) %>%
		cbind(groupSize = rep(g,lifetime*5)) %>%
		cbind(round = rep(1:lifetime, 5)) %>%
		cbind(value = c(rep('optimalChoiceProb', lifetime),rep('expectedPerformance', lifetime),rep('payoff', lifetime),rep('performance', lifetime),rep('optimalChoiceFreq', lifetime))) %>%
		cbind(indiv_small_large = rep('Large', lifetime*5))
	)


posthoc_sim_summary_indivSmallLarge_25$taskDifficulty3 = 'Low Uncertainty'
posthoc_sim_summary_indivSmallLarge_50$taskDifficulty3 = 'Moderate Uncertainty'
posthoc_sim_summary_indivSmallLarge_78$taskDifficulty3 = 'High Uncertainty'

posthoc_sim_summary_indivSmallLarge = rbind(posthoc_sim_summary_indivSmallLarge_25, posthoc_sim_summary_indivSmallLarge_50) %>% rbind(posthoc_sim_summary_indivSmallLarge_78)

posthoc_sim_summary_indivSmallLarge$taskDifficulty3 = factor(posthoc_sim_summary_indivSmallLarge$taskDifficulty3, levels=c('Low Uncertainty', 'Moderate Uncertainty', 'High Uncertainty'))


write.csv(posthoc_sim_summary_indivSmallLarge,
			"~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/posthoc_sim_summary_indivSmallLarge.csv",
			row.names=FALSE)


ggplot() +
	geom_ribbon(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Individual'), mapping=aes(round, ymin=X25., ymax=X75.), fill='black', alpha = 1/6) +
	geom_line(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Individual'), mapping=aes(round, y=mean), colour='black', size=1, linetype='dashed') +
	geom_ribbon(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Small'), mapping=aes(round, ymin=X25., ymax=X75.), fill='orange', alpha = 1/6) +
	geom_line(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Small'), mapping=aes(round, y=mean), colour='orange', size=1, linetype='twodash') +
	geom_ribbon(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Large'), mapping=aes(round, ymin=X25., ymax=X75.), fill='red', alpha = 1/6) +
	geom_line(data=posthoc_sim_summary_indivSmallLarge %>% dplyr::filter(value=="optimalChoiceProb"&indiv_small_large=='Large'), mapping=aes(round, y=mean), colour='red', size=1, linetype='solid') +
	facet_grid(. ~ taskDifficulty3) +
	scale_y_continuous(limits = c(-0.05,1.05), breaks=c(0.0,0.25,0.75,1.0)) +
	geom_vline(xintercept=70)+
	labs(x='Rounds', y='Proportion of choosing\n the best option')+
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
	geom_hline(yintercept=1/3, linetype="dashed")+
	myTheme()+
	theme(legend.position='none',
	axis.text.y = element_text(size=11, family="Times" ,colour='black'),
	strip.text.y = element_text(angle = 0))+
	NULL


## Each group's average performance
## Does the performance distribution shape unimordal or bimodal?
## Environment 1
## _25%
posthoc_sim_eachGroup_aveOptim_env1_indiv_25 =  apply(posthoc_sim_result_25[[paste("n=",1)]][,1:(changingPoint-1)], 2, mean) %>% data.frame(taskDifficulty='25%',indiv_small_large='Single',environment='1st - 40th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env1_small_25 =  apply(( rbind(posthoc_sim_result_25[[paste("n=",3)]],posthoc_sim_result_25[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",6)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",9)]]) )[,1:(changingPoint-1)], 2, mean) %>% data.frame(taskDifficulty='25%',indiv_small_large='Small',environment='1st - 40th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env1_large_25 =  apply(( rbind(posthoc_sim_result_25[[paste("n=",12)]],posthoc_sim_result_25[[paste("n=",15)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",16)]]) )[,1:(changingPoint-1)], 2, mean) %>% data.frame(taskDifficulty='25%',indiv_small_large='Large',environment='1st - 40th') %>% dplyr::rename(meanOptimProb = ".")


## _50%
posthoc_sim_eachGroup_aveOptim_env1_indiv_50 =  apply(posthoc_sim_result_50[[paste("n=",1)]][,1:(changingPoint-1)], 2, mean) %>% data.frame(taskDifficulty='50%',indiv_small_large='Single',environment='1st - 40th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env1_small_50 =  apply(( rbind(posthoc_sim_result_50[[paste("n=",2)]],posthoc_sim_result_50[[paste("n=",3)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",4)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",6)]]) )[,1:(changingPoint-1)], 2, mean) %>% data.frame(taskDifficulty='50%',indiv_small_large='Small',environment='1st - 40th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env1_large_50 =  apply(( rbind(posthoc_sim_result_50[[paste("n=",8)]],posthoc_sim_result_50[[paste("n=",9)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",17)]]) )[,1:(changingPoint-1)], 2, mean) %>% data.frame(taskDifficulty='50%',indiv_small_large='Large',environment='1st - 40th') %>% dplyr::rename(meanOptimProb = ".")


## _78%
posthoc_sim_eachGroup_aveOptim_env1_indiv_78 =  apply(posthoc_sim_result_78[[paste("n=",1)]][,1:(changingPoint-1)], 2, mean) %>% data.frame(taskDifficulty='78%',indiv_small_large='Single',environment='1st - 40th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env1_small_78 =  apply(( rbind(posthoc_sim_result_78[[paste("n=",2)]],posthoc_sim_result_78[[paste("n=",3)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",6)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",7)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",8)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",9)]]) )[,1:(changingPoint-1)], 2, mean) %>% data.frame(taskDifficulty='78%',indiv_small_large='Small',environment='1st - 40th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env1_large_78 =  apply(( rbind(posthoc_sim_result_78[[paste("n=",12)]],posthoc_sim_result_78[[paste("n=",14)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",15)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",18)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",20)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",23)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",24)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",27)]]) )[,1:(changingPoint-1)], 2, mean) %>% data.frame(taskDifficulty='78%',indiv_small_large='Large',environment='1st - 40th') %>% dplyr::rename(meanOptimProb = ".")


## Environment 2
## _25%
posthoc_sim_eachGroup_aveOptim_env2_indiv_25 =  apply(posthoc_sim_result_25[[paste("n=",1)]][,changingPoint:70], 2, mean) %>% data.frame(taskDifficulty='25%',indiv_small_large='Single',environment='41st - 70th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env2_small_25 =  apply(( rbind(posthoc_sim_result_25[[paste("n=",3)]],posthoc_sim_result_25[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",6)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",9)]]) )[,changingPoint:70], 2, mean) %>% data.frame(taskDifficulty='25%',indiv_small_large='Small',environment='41st - 70th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env2_large_25 =  apply(( rbind(posthoc_sim_result_25[[paste("n=",12)]],posthoc_sim_result_25[[paste("n=",15)]]) %>% rbind(posthoc_sim_result_25[[paste("n=",16)]]) )[,changingPoint:70], 2, mean) %>% data.frame(taskDifficulty='25%',indiv_small_large='Large',environment='41st - 70th') %>% dplyr::rename(meanOptimProb = ".")

## _50%
posthoc_sim_eachGroup_aveOptim_env2_indiv_50 =  apply(posthoc_sim_result_50[[paste("n=",1)]][,changingPoint:70], 2, mean) %>% data.frame(taskDifficulty='50%',indiv_small_large='Single',environment='41st - 70th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env2_small_50 =  apply(( rbind(posthoc_sim_result_50[[paste("n=",2)]],posthoc_sim_result_50[[paste("n=",3)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",4)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",6)]]) )[,changingPoint:70], 2, mean) %>% data.frame(taskDifficulty='50%',indiv_small_large='Small',environment='41st - 70th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env2_large_50 =  apply(( rbind(posthoc_sim_result_50[[paste("n=",8)]],posthoc_sim_result_50[[paste("n=",9)]]) %>% rbind(posthoc_sim_result_50[[paste("n=",17)]]) )[,changingPoint:70], 2, mean) %>% data.frame(taskDifficulty='50%',indiv_small_large='Large',environment='41st - 70th') %>% dplyr::rename(meanOptimProb = ".")

## _78%
posthoc_sim_eachGroup_aveOptim_env2_indiv_78 =  apply(posthoc_sim_result_78[[paste("n=",1)]][,changingPoint:70], 2, mean) %>% data.frame(taskDifficulty='78%',indiv_small_large='Single',environment='41st - 70th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env2_small_78 =  apply(( rbind(posthoc_sim_result_78[[paste("n=",2)]],posthoc_sim_result_78[[paste("n=",3)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",5)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",6)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",7)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",8)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",9)]]) )[,changingPoint:70], 2, mean) %>% data.frame(taskDifficulty='78%',indiv_small_large='Small',environment='41st - 70th') %>% dplyr::rename(meanOptimProb = ".")
posthoc_sim_eachGroup_aveOptim_env2_large_78 =  apply(( rbind(posthoc_sim_result_78[[paste("n=",12)]],posthoc_sim_result_78[[paste("n=",14)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",15)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",18)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",20)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",23)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",24)]]) %>% rbind(posthoc_sim_result_78[[paste("n=",27)]]) )[,changingPoint:70], 2, mean) %>% data.frame(taskDifficulty='78%',indiv_small_large='Large',environment='41st - 70th') %>% dplyr::rename(meanOptimProb = ".")

posthoc_sim_eachGroup_aveOptim =
	rbind(posthoc_sim_eachGroup_aveOptim_env1_indiv_25,posthoc_sim_eachGroup_aveOptim_env1_small_25) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env1_large_25) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env1_indiv_50) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env1_small_50) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env1_large_50) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env1_indiv_78) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env1_small_78) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env1_large_78) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env2_indiv_25) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env2_small_25) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env2_large_25) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env2_indiv_50) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env2_small_50) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env2_large_50) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env2_indiv_78) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env2_small_78) %>%
	rbind(posthoc_sim_eachGroup_aveOptim_env2_large_78)
posthoc_sim_eachGroup_aveOptim$taskDifficulty3 = 'Low Uncertainty'
posthoc_sim_eachGroup_aveOptim$taskDifficulty3[which(posthoc_sim_eachGroup_aveOptim$taskDifficulty=='50%')] = 'Moderate Uncertainty'
posthoc_sim_eachGroup_aveOptim$taskDifficulty3[which(posthoc_sim_eachGroup_aveOptim$taskDifficulty=='78%')] = 'High Uncertainty'
posthoc_sim_eachGroup_aveOptim$taskDifficulty3 = factor(posthoc_sim_eachGroup_aveOptim$taskDifficulty3, levels=c('Low Uncertainty', 'Moderate Uncertainty', 'High Uncertainty'))

write.csv(posthoc_sim_eachGroup_aveOptim,
			"~/Dropbox/wataru/St\ Andrews/onlineExperiment/st-andrews-webexp2/analysis/stan_mdelFitting/posthoc_sim_eachGroup_aveOptim.csv",
			row.names=FALSE)

ggplot() +
	geom_density(data = posthoc_sim_eachGroup_aveOptim, mapping=aes(x=meanOptimProb, y=..scaled.., colour=indiv_small_large,fill=indiv_small_large,linetype=indiv_small_large),alpha=1/3) +
	labs(y='Proportion',x='Groups\' average proportion of \nchoosing the best option')+
	scale_colour_manual(values=c("Single"='black','Small'='orange', 'Large'='red'))+
	scale_fill_manual(values=c("Single"='black','Small'='orange', 'Large'='red'))+
	scale_linetype_manual(values=c("Single"='dashed','Small'='twodash', 'Large'='solid'))+
	scale_y_continuous(limits = c(-0.05,1.05), breaks=c(0.0,0.5,1.0)) +
	scale_x_continuous(limits = c(-0.05,1.05), breaks=c(0.0,0.5,1.0)) +
	panel_border(colour = "black", size = 0.5, linetype = 1,remove = FALSE) +
	myTheme_small()+
	theme(axis.text.x = element_text(size=9, family="Times" ,colour='black'),
		axis.text.y = element_text(size=9, family="Times" ,colour='black'))+
	facet_grid(.~taskDifficulty3 + environment)+
	NULL
























































