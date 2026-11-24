# np-errorprone-aux

The R scripts execute Monte Carlo simulation studies comparing the performance of several methods for making inference using non-probability samples when there is error in some of the auxiliary variables of the reference probability sample.

## KernelSmoothing_for_ME_nonprob.R

The script KernelSmoothing_for_ME_nonprob.R examines 4 different covariate measurement error scenarios.  The user will need to comment/uncomment the appropriate scenarios in the following lines:
* Line 637: Classical measurement error scenario
* Line 640: Downwards-biased measurement error
* Line 643: Classical measurement error scenario with larger variance than Scenario 1
* Line 646: Non-linear measurement error

Additional parameters that can be changed include:
* Line 626: Random seed
* Line 627: Boolean variable to determine if bootstrap variance estimates are produced
* Line 628: Number of bootstrap replicates
* Line 630: Number of Monte Carlo simulations

Comparative Relative Bias and Relative Root Mean Square Error results are produced in the output1 dataset.  Additionally, the coverage rate and average confidence interval length for the bootstrap variance estimator are contained in the *c_rate* and *al* variables (lines 1104 and 1107).

## KernelSmoothing_SNAR_Scenario.R

The script KernelSmoothing_SNAR_Scenario.R examines the extension of the measurement error scenario where one of the required covariates in the model for participation probability is the data item of interest, and it is missing.

Simulation parameters that can be changed include:
* Line 611: Random seed
* Line 612: Boolean variable to determine if bootstrap variance estimates are produced
* Line 613: Number of bootstrap replicates
* Line 615: Number of Monte Carlo simulations
* Lines 620 and 621: Standard deviation of the x covariates

Comparative Relative Bias and Relative Root Mean Square Error results are produced in the output1 dataset.  Additionally, the coverage rate and average confidence interval length for the bootstrap variance estimator are contained in the *c_rate* and *al* variables (lines 1061 and 1064).
 

## Reference

The results of the simulation studies are contained in the paper:

Ang, L., R. Clark, B. Loong, and A. Holmberg. 2026. "Inference for Non-Probability Samples Using Error-Prone Auxiliary Information". In submission. Pre-print available at http://arxiv.org/abs/2611.10975
