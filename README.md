# ToyokawaWhalenLaland2018
Individual-based model and statistical code used in an article entitled "Social learning strategies regulate the wisdom and madness of interactive online crowds"

All the data used are also included in this repository. You can draw all figures shown in the main text of the manuscript using `ToyokawaWhalenLaland2018_figuresScript.R`, by running the script from the top to the end.

## Agent-based model simulation
- The simulation was written in a Mathematica notebook (`IBM_simulation.nb`)
 
## Computational models
- All code used to fit the main model can be found in `model_fitting_for_exp.R`
- For the alternative model (where both social learning parameters are time dependent) can be found in `model_fitting_fullmodel_for_supp.R`
- Full details of the model are written in stan (`model_UNC6_sReduc_annealing.stan`).

## Parameter recovery test
- Synthetic data are generated in `parameterRecovery/parameterRecoverySimulation.nb`
- Then, you can fit the computational model using `parameterRecovery/model_UNC6_sReduc_annealing_fitting_for_sim.R`.
