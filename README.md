# np-errorprone-aux

This R code executes a Monte Carlo simulation study comparing the performance of several methods for making inference using non-probability samples when there is error in some of the auxiliary variables of the reference probability sample.

The code will produce results for 4 different measurement error scenarios.  The user will need to comment/uncomment the appropriate scenarios in the following lines:
* Line 637: Classical measurement error scenario
* Line 640: Downwards-biased measurement error
* Line 643: Classical measurement error scenario with larger variance than Scenario 1
* Line 646: Non-linear measurement error

The code KernelSmoothing_SNAR_Scenario.R examines the extension of the measurement error scenario where the required covariate (in this case the variable of interest) is completely missing.

The results of the simulation study in KernelSmoothing_SNAR_Scenario.R are contained in the paper:

Ang, L., R. Clark, B. Loong, and A. Holmberg. 2026. "Inference for Non-Probability Samples Using Error-Prone Auxiliary Information". In submission.
