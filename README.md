# ToyokawaWhalenLaland2018
Individual-based model and statistical code used in an article entitled "Social learning strategies regulate the wisdom and madness of interactive crowds"

The data used are basically included in this repository, but some are missing due to their size. You could generate them by running simulaiton and analyses in your side. WT will upload all the data to somewhere open, but in the mean time, you can contact him to get them. 

Once you have all the data set (`"./data/allBehaviouralData2.csv"`), you can draw all figures shown in the main text of the manuscript using `ToyokawaWhalenLaland2018_figuresScript_revision.R`, by running the script from the top to the end.

## Agent-based model simulation
- The simulation was written in a Mathematica notebook (`IBM_simulation.nb`)
 
## Computational models
- All code used to fit the main model can be found in `model_fitting_for_exp.R`
- For the alternative model (where both social learning parameters are time dependent) can be found in `model_fitting_fullmodel_for_supp.R`
- Full details of the model are written in stan files (e.g. `model_UNC6_sReduc_annealing.stan`).

## Parameter recovery test
- Synthetic data are generated in `parameterRecovery/parameterRecoverySimulation.nb`
- Then, you can fit the computational model using `parameterRecovery/model_UNC6_sReduc_annealing_fitting_for_sim.R`.
