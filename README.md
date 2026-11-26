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
* Line 669: Correlation between the variable of interest $Y$ and $X\beta$
* Line 670: Sample size of the non-probability sample
* Line 671: Sample size of the probability sample

To get Tables 1 to 4 in the paper referenced below, the user will need to:
* Run the script by un-commenting the appropriate scenario:
  * Line 637 for Table 1
  * Line 640 for Table 2
  * Line 643 for Table 3
  * Line 646 for Table 4
* Set the sample size of the non-probability sample to 2,000 or 8,000 as required
* Set the sample size of the probability sample to 1,000
* Set the correlation to 0.3
* The number of Monte Carlo simulations is set to 2,000
* The number of Bootstrap replicates is set to 50

Comparative Relative Bias and Relative Root Mean Square Error results are produced in the output1 dataset.  Additionally, the coverage rate and average confidence interval length for the bootstrap variance estimator are contained in the *c_rate* and *al* variables (lines 1104 and 1107).


## KernelSmoothing_SNAR_Scenario.R

The script KernelSmoothing_SNAR_Scenario.R examines the extension of the measurement error scenario where one of the required covariates in the model for participation probability is the data item of interest, and it is missing.

Simulation parameters that can be changed include:
* Line 611: Random seed
* Line 612: Boolean variable to determine if bootstrap variance estimates are produced
* Line 613: Number of bootstrap replicates
* Line 615: Number of Monte Carlo simulations
* Lines 620 and 621: Standard deviation of the x covariates

The simulation settings we used included:
* Set population size to 5,000
* Set probability sample size to 1,000 and non-probability sample size to 2,500
* Set standard deviation of the $x_1$ covariate to 0.3 or 1
* Nunmber of simulations was 1,000
* Number of bootstrap replicates was 50

Comparative Relative Bias and Relative Root Mean Square Error results are produced in the output1 dataset.  Additionally, the coverage rate and average confidence interval length for the bootstrap variance estimator are contained in the *c_rate* and *al* variables (lines 1061 and 1064).
 

## Reference

The results of the simulation studies from the first R script are contained in the paper:

Ang, L., R. Clark, B. Loong, and A. Holmberg. 2026. "Inference for Non-Probability Samples Using Error-Prone Auxiliary Information". In submission. Pre-print available at http://arxiv.org/abs/2611.10975
