########################################################################
### This code executes a simulation comparing a number of methods for
### producing propensities of selection for a non-probability dataset when
### some of the auxiliary variables in the reference probability sample
### suffer from measurement error.
###
### This code proposes an approach which takes advantage of the ability
### to link the reference probability sample and the non-probability sample
### to account for the measurement error and estimate participation probabilities.
###
### The simulation setup draws inspiration from a paper by Chen, Li and Wu (2020)
### exploring doubly robust estimation for non-probability samples.  I acknowledge
### the authors of that paper for sharing their original code.
### Any errors in this program are my own.
###
### 
### Author: Lyndon Ang
###
### Date: November 2025
###
######################################################################

### Note: In this code, the dataset A is the non-probability sample, while
###       the dataset B represents the probability sample.

rm(list=ls())

library("xtable")
library("VIM")
library("sampling")
#library("ks")
library("np")


### Libraries to facilitate parallel computing
library(foreach)
library(doParallel)
library(doRNG)

starttime <- Sys.time()

### Set the number of core processors to use
registerDoParallel(cores=4)



#################################################
### Functions used in the simulation
#################################################


#randomized systematic sampling

syspps=function(N_f,x){ 
  ##Population is first randomized
  U=sample(N_f,N_f)
  xx=x[U]
  cx=cumsum(xx)
  r=runif(1)
  s=numeric()
  for(i in 1:N_f){
    if(cx[i]>=r){
      s=c(s,U[i])
      r=r+1
    }
  }
  return(s)
}


#Estimate theta with Chen et al method

hat_theta_f=function(xsa,xsb,wsb){
  #xsa <- X_SA_res
  #xsb <- X_SB_res
  #wsb <- w_SB
  
  col_xsa=colSums(xsa)
  theta_0=rep(0,length(col_xsa))
  theta_1=solve(t(xsb)%*%(wsb*xsb),col_xsa-t(xsb)%*%wsb)
  while (abs(max(theta_1-theta_0))>10^(-8)){
    theta_0=theta_1
    ps=1/(1+exp(-xsb%*%theta_0))
    theta_1=theta_0+solve(t(xsb)%*%(c(wsb*ps*(1-ps))*xsb),col_xsa-t(xsb)%*%c(wsb*ps))
  }
  
  return(theta_1)
}


### Estimate theta using subset of B excluding A

hat_theta2_f <- function(xsa, xsb, wsb) {
  #xsa <- X_SA_res
  #xsb <- X_SBnA_res
  #wsb <- w_SBnA
  
  col_xsa <- colSums(xsa) ### Calculate column sums of x variables in A
  theta_0 <- rep(0.1, length(col_xsa)) ### Initial values for theta
  theta_1 <- theta_0 - 0.1
  
  while(abs(max(theta_1 - theta_0)) > 10^(-8)) {
    theta_0 <- theta_1
    ps_A <- 1/(1 + exp(-xsa %*% theta_0)) ### Calculate propensity for units in A
    ps_BnA <- 1/(1 + exp(-xsb %*% theta_0)) ### Calculate propensity for units in B not A
    
    ### Now find the next iteration of theta
    theta_1 <- theta_0 + solve(t(xsa) %*% (c(ps_A*(1 - ps_A)) * xsa) + t(xsb)%*%(c(wsb*ps_BnA*(1-ps_BnA))*xsb),
                               col_xsa - t(xsa) %*% ps_A - t(xsb)%*%c(wsb*ps_BnA))
    
  }
  
  return(theta_1)
  
}


### Estimate theta using A and B linked together (Kim and Wang method)

hat_theta_KW_f <- function(xsbpa, xsbna, wsbpa, wsbna) {
  
  col_xsbpa <- t(xsbpa) %*% wsbpa  ### Calculate weighted column sums of x variables in B
  theta_0 <- rep(0.1, length(col_xsbpa)) ### Initial values for theta
  theta_1 <- theta_0 - 0.1
  
  while(abs(max(theta_1 - theta_0)) > 10^(-8)) {
    theta_0 <- theta_1
    ps_BpA <- 1/(1 + exp(-xsbpa %*% theta_0)) ### Calculate propensity for units in B linked to A
    ps_BnA <- 1/(1 + exp(-xsbna %*% theta_0)) ### Calculate propensity for units in B not A
    
    ### Now find the next iteration of theta
    theta_1 <- theta_0 + solve(t(xsbpa) %*% (c(wsbpa*ps_BpA*(1 - ps_BpA)) * xsbpa) + t(xsbna)%*%(c(wsbna*ps_BnA*(1-ps_BnA))*xsbna),
                               col_xsbpa - t(xsbpa) %*% c(wsbpa*ps_BpA) - t(xsbna)%*%c(wsbna*ps_BnA))
    
  }
  
  return(theta_1)
  
}


### Multivariate normal kernel function
mvn_kernel_r <- function(xi, x, H) {
  diff <- x - xi
  kernel_value <- exp(-0.5 * rowSums(diff %*% solve(H) * diff)) / sqrt((2 * pi)^ncol(x) * det(H))
  return(kernel_value)
}

### Logistic function
logistic_r <- function(X, theta) {
  return (1 / (1 + exp(-X %*% theta)))
}

### Function to execute the E step
e_step_R <- function(ZX, theta_em, krn, rowAB, wgtAB) {
  # nrowA <- rowBnA
  # wgtAB <- wijwgt
  
  # Calculate eta
  eta <- ZX %*% theta_em
  
  # Calculate the Odds O(x,y)
  O_ZX <- exp(-eta)
  
  # Calculate the Expected Odds
  EOk <- O_ZX * krn * wgtAB
  
  # Reshape the odds vector to a matrix, each column corresponds to a different set of z
  prod_matrix <- matrix(EOk, nrow = rowAB, byrow = FALSE)
  
  # Calculate the column sums
  col_sums <- colSums(prod_matrix)
  
  # Divide the cells in the row by column sums
  wij <- sweep(prod_matrix, 2, col_sums, "/")
  
  # Return the vectorized wij
  return(as.vector(wij))
}

### Function to execute the M step

m_step_loop_R <- function(xsa, ZX, theta_init, col_xsa, wsb, wij_vec, tol = 1e-6, max_iter = 100) {
  #theta_init <- hat_theta_me_init
  #wij_vec <- wij.vec
  
  theta_0 <- theta_init - 0.1
  theta_1 <- theta_init
  
  for (iter in 1:max_iter) {
    theta_0 <- theta_1
    
    # Calculate propensities
    ps_A <- 1 / (1 + exp(-xsa %*% theta_0))
    ps_BnA <- 1 / (1 + exp(-ZX %*% theta_0))
    
    # Calculate weights
    w1 <- ps_BnA * wsb * wij_vec
    w2 <- ps_BnA * (1 - ps_BnA) * wsb * wij_vec
    
    
    # Calculate gradient
    grad <- col_xsa - t(xsa) %*% ps_A - t(ZX) %*% w1
    
    # Calculate Hessians
    hess_A <- t(xsa) %*% (c(ps_A * (1 - ps_A))*xsa)
    hess_B <- t(ZX) %*% (c(w2)*ZX)
    hess <- hess_A + hess_B
    
    # Update theta
    delta <- solve(hess, grad)
    theta_1 <- theta_0 + delta
    
    # Check for convergence
    if (sqrt(sum((theta_1 - theta_0)^2)) < tol) {
      break
    }
  }
  
  return(theta_1)
}



### Function to estimate theta assuming SNAR using EM algorithm - assume can link
hat_theta_me_fn <- function(hat_theta_init, ZX, krn, xsa, rowAB, ysa, wsb, wijwgt, printout=FALSE) {
  
  ### Set up initial values for theta
  
  mod.coef.old <- hat_theta_init
  mod.coef.new <- hat_theta_init + 2
  theta_em <- hat_theta_init
  
  i <- 0
  allest <- numeric()
  col_xsa <- colSums(xsa)
  
  ### EM iterations
  while (abs(max(mod.coef.new-mod.coef.old)) > 10^(-3)) {
    mod.coef.old <- theta_em
    i<-i+1
    
    ### Run the E-step function (port out to C++)
    
    #wij.vec <- e_step(ZX, theta_em, krn, rowAB, wijwgt)
    wij.vec <- e_step_R(ZX, theta_em, krn, rowAB, wijwgt)
    
    ### Now the M step - undertakes the Newton-Raphson method to find updated theta
    
    #theta_1 <- m_step_loop(xsa, ZX, theta_em, colSums(xsa), wsb, wij.vec)
    theta_1<- m_step_loop_R(xsa, ZX, theta_em, colSums(xsa), wsb, wij.vec)
    
    
    ### If printout is TRUE, print the current theta estimate and iteration number
    if(printout) {
      print(paste("Iteration:", i, "Theta:", paste(round(theta_1, 4), collapse=", ")))
    }
    
    theta_em <- theta_1
    
    hat_score <- logistic_r(xsa, theta_em)
    mu_new <- sum(ysa/hat_score)/sum(1/hat_score)
    
    ### Record the estimates for each EM iteration for inspection if needed
    allest[i] <- mu_new
    
    mod.coef.new <- theta_em
    
    ### To stop infinite loops
    if(i>100) break
    
  }
  
  results <- list("theta" = theta_em, "allest" = allest, "EM.iter" = i)
  
  return(results)
  
}


hat_theta_me_fn_old <- function(hat_theta_init, ZX, krn, xsa, ysa, wsb, printout=FALSE) {
  
  ### Set up initial values for theta
  
  mod.coef.old <- hat_theta_init
  mod.coef.new <- hat_theta_init + 2
  theta_em <- hat_theta_init
  
  i <- 0
  allest <- numeric()
  
  ### EM iterations
  while (abs(max(mod.coef.new-mod.coef.old)) > 10^(-3)) {
    mod.coef.old <- theta_em
    i<-i+1
    
    ### E-step
    
    wij.vec <- e_step(ZX, theta_em, krn, nrow(xsa))
    
    #wsb <- dwgts
    
    ### Now the M step
    
    col_xsa <- colSums(xsa) ### Calculate column sums of x variables in A
    theta_0 <- theta_em - 0.1 ### Initial values for theta
    theta_1 <- theta_em
    
    j <- 0
    
    while(abs(max(theta_1 - theta_0)) > 10^(-6)) {
      j <- j + 1
      theta_0 <- theta_1
      
      ps_A <- logistic_r (xsa, theta_0) ### Calculate propensity for units in A
      ps_BnA <- logistic_r(ZX, theta_0)
      
      xpw <- t(ZX) %*% c((ps_BnA)*wsb*wij.vec)
      
      xpw2 <- t(ZX) %*% (c(ps_BnA*(1 - ps_BnA)*wsb*wij.vec)*ZX)
      
      ### Now find the next iteration of theta
      theta_1 <- theta_0 + solve(t(xsa) %*% (c(ps_A*(1 - ps_A)) * xsa) + xpw2,
                                 col_xsa - t(xsa) %*% ps_A - xpw)
      
      ### If printout is TRUE, print the current theta estimate and iteration number
      if(printout) {
        print(paste("Iteration:", i, "Theta:", paste(round(theta_1, 4), collapse=", ")))
      }
      #if(j>100) break
    }
    
    
    theta_em <- theta_1
    
    hat_score <- 1/(1 + exp(-xsa%*%theta_em))
    mu_new <- sum(ysa/hat_score)/sum(1/hat_score)
    
    ### Record the estimates for each EM iteration for inspection if needed
    allest[i] <- mu_new
    
    mod.coef.new <- theta_em
    
    ### To stop infinite loops
    if(i>100) break
    
  }
  
  ### Put results on to list
  results <- list("theta" = theta_em, "allest" = allest, "EM.iter" = i)
  
  return(results)
  
}

### Function which we will use as input to the optimise function, to find optimal bandwidth
### This function is similar to the calc_theta_pwr function below, except that this function
### also includes the calculation of the error norm

build_objective_fn <- function(
    ZXs_AB, wgtk, ZXs_BA.matrix, ZXs_BnA.repmatrix,
    w_SBnA, samp_AB, ZX, xsa, Y_SA,
    ZXs_AB.dupmatrix, ZXs_AB.repmatrix,
    ZX_AB_matrix, ZX_AB
) {
  
  function(pwr_val) {
    
    ### Calculate the bandwidth using the given power pwr_val
    h <- apply(ZXs_AB[,-1], 2, sd) * nrow(ZXs_AB)^(-pwr_val)
    ### Calculate kernel value for the ij combinations
    #krn <- wgtk * mvn_kernel_r(ZXs_BA.matrix, ZXs_BnA.repmatrix, diag(h))
    krn <- mvn_kernel_r(ZXs_BA.matrix, ZXs_BnA.repmatrix, diag(h))
    ### Weights for B not A - each one duplicated nrow(AB) times
    dwgts <- rep(w_SBnA, each = nrow(samp_AB))
    
    ### Initialise hat_theta
    hat_theta_me_init <- rep(0, ncol(ZX))
    
    ### Get number or rows in BnA
    rowAB <- nrow(samp_AB)
    
    ### Now calculate theta_hat
    me_em <- hat_theta_me_fn(hat_theta_me_init, ZX, krn, xsa, rowAB, Y_SA, dwgts, wgtk, printout = FALSE)
    
    ### Find the distribution of zx in A - this is our f1(x)
    h1 <- apply(ZXs_AB[,-1], 2, sd) * nrow(ZXs_AB)^(-0.4)
    f1x_comp <- mvn_kernel_r(ZXs_AB.dupmatrix, ZXs_AB.repmatrix, diag(h1))
    f1x <- colSums(matrix(f1x_comp, nrow = rowAB, byrow = FALSE))
    
    ### Kernel values K(A,A)
    krnf1x <- mvn_kernel_r(ZXs_AB.dupmatrix, ZXs_AB.repmatrix, diag(h))
    
    ### Now calculate the error using the chosen bandwidth h
    
    # Odds using A
    O_ZXtrue <- exp(-ZX_AB %*% me_em$theta)
    O_ZX <- exp(-ZX_AB_matrix %*% me_em$theta)
    EOk <- O_ZX * as.vector(krnf1x)
    ### Place into matrix now where each set of z values is down a column (and rows are the number of rows in A)
    prod.matrix <- matrix(EOk, nrow = rowAB, byrow = FALSE)
    D_hat <- colSums(prod.matrix)
    
    ### Difference for D
    
    diff_D <- (D_hat / f1x - O_ZXtrue)^2
    
    ### Difference for C
    
    ### Calculate propensities
    ps_B <- logistic_r(ZX_AB_matrix, me_em$theta)
    psB.EOk <- matrix(ps_B, nrow = rowAB, byrow = FALSE) * prod.matrix
    
    ps_Btrue <- logistic_r(ZX_AB, me_em$theta)
    
    ### Calculate C_hat for each x covariate (column) separately
    C_hat <- sapply(1:ncol(ZX_AB_matrix), function(j) {
      Zj <- matrix(ZX_AB_matrix[, j], nrow = rowAB, byrow = FALSE)
      colSums(Zj * psB.EOk) / f1x
    })
    
    ### Calculate true C value for each x covariate (column) separately
    C_true <- ZX_AB * c(ps_Btrue * O_ZXtrue)
    diff_C <- sapply(1:ncol(C_hat), function(j) (C_hat[,j] - C_true[, j])^2)
    
    ### Calculate the norm
    normK <- sum((sqrt(diff_D + rowSums(diff_C)))^2)
    return(normK)
  }
}

### Function to calculate theta using the optimal power
calc_theta_pwr <- function(
    pwr_val,
    ZXs_AB, wgtk, ZXs_BA.matrix, ZXs_BnA.repmatrix,
    w_SBnA, samp_AB, ZX, xsa, Y_SA
) {
  ### Calculate the bandwidth using the given power pwr_val
  h <- apply(ZXs_AB[,-1], 2, sd) * nrow(ZXs_AB)^(-pwr_val)
  ### Calculate kernel value for the ij combinations
  krn <- mvn_kernel_r(ZXs_BA.matrix, ZXs_BnA.repmatrix, diag(h))
  ### Weights for B not A - each one duplicated nrow(B) times
  dwgts <- rep(w_SBnA, each = nrow(samp_AB))
  
  ### Initialise hat_theta
  hat_theta_me_init <- rep(0, ncol(xsa))
  rowAB <- nrow(samp_AB)
  
  ### Now calculate theta_hat
  me_em <- hat_theta_me_fn(hat_theta_me_init, ZX, krn, xsa, rowAB, Y_SA, dwgts, wgtk, printout = FALSE)
  
  return(me_em)
}


### Function to calculate bootstrap estimates

boot_var_me <- function(
    samp_A, X_SA_reg, Y_SA, hat_score_ME_SA, 
    N, mean_n_SB, B, pwr)
{
  
  # Create the pseudo population using Booth et al (1994) approach
  # Number of duplicates based on inverse probability of being in the NP dataset
  np.wgts <- 1/hat_score_ME_SA
  pp.prob <- (np.wgts/sum(np.wgts))*N
  pp.prob[pp.prob<1] <- 1
  duptimes_floor <- floor(pp.prob)
  
  # Create rowindex which will duplicate items in NP sample based on inverse probability
  rowindex <- rep(1:nrow(samp_A), duptimes_floor)
  
  # Create first segment of the pseudo population based on the duplicates
  Uf <- samp_A[rowindex,]
  Uf.probselA <- hat_score_ME_SA[rowindex]
  Uf.X_reg <- X_SA_reg[rowindex,]
  Uf.Y <- Y_SA[rowindex]
  
  # Take a PPS sample with probability based on remainder from pp.prob from samp_A
  pps_remainder <- pp.prob - duptimes_floor
  remsamp <- rbinom(nrow(samp_A), 1, pps_remainder)
  Ue_index <- order(-remsamp)[1:sum(remsamp)]
  
  Ue <- samp_A[Ue_index,]
  Ue.probselA <- hat_score_ME_SA[Ue_index]
  Ue.X_reg <- X_SA_reg[Ue_index,]
  Ue.Y <- Y_SA[Ue_index]
  
  # Create final pseudo population
  pseudopopn <- rbind(Uf, Ue)
  probselA <- c(Uf.probselA, Ue.probselA)
  ppX_reg <- rbind(Uf.X_reg, Ue.X_reg)
  ppY <- c(Uf.Y, Ue.Y)
  
  pseudopopn$id <- 1:nrow(pseudopopn)
  pseudopopn$probselA <- probselA
  
  # Initialize bootstrap estimates vector
  byest <- numeric(B)
  
  # Run bootstrap iterations in parallel if possible
  for (bi in 1:B) {
    # 1 - Select the NP sample A
    bsam <- rbinom(nrow(pseudopopn), 1, probselA)
    bnSA <- sum(bsam)
    bsa_index <- order(-bsam)[1:bnSA]
    bxsa_reg <- ppX_reg[bsa_index,]
    bysa <- ppY[bsa_index]
    bsa <- pseudopopn[bsa_index,]
    
    # 2 - Select the P sample B using Pareto sampling
    bz1 <- (1/pseudopopn$w)/sum(1/pseudopopn$w)
    bpik <- sampling::inclusionprobabilities(bz1, mean_n_SB)
    omega <- runif(nrow(pseudopopn))
    lambda <- bpik
    bsb_index <- order(omega * (1 - lambda)/((1 - omega) * lambda))[1:mean_n_SB]
    bsb <- pseudopopn[bsb_index,]
    
    # 3 - Get datasets needed for estimation
    bx_ab <- as.matrix(bsb[bsb$id %in% bsa$id, c("x4")])
    bzxsab_all <- as.matrix(bsb[bsb$id %in% bsa$id, c("V1","x1","x2","x3","x4me_boot")])
    bzbna <- as.matrix(bsb[!(bsb$id %in% bsa$id), c("V1","x1","x2","x3")])
    bzxsbna_all <- as.matrix(bsb[!(bsb$id %in% bsa$id), c("V1","x1","x2","x3","x4me_boot")])
    bz_ab <- as.matrix(bsb[bsb$id %in% bsa$id, c("V1","x1","x2","x3")])
    bzxsab_all <- as.matrix(bsb[bsb$id %in% bsa$id, c("V1","x1","x2","x3","x4me_boot")])
    bzxab_all <- as.matrix(bsb[bsb$id %in% bsa$id, c("V1","x1","x2","x3","x4")])
    
    bsamp_bna <- bsb[!(bsb$id %in% bsa$id),]
    bsamp_ab <- bsb[bsb$id %in% bsa$id,]
    
    ### Weights for segments of sample B
    bw_sbna <- bsamp_bna$w
    bw_sba <- bsamp_ab$w
    
    # Create matrices for ij combinations
    bz.repmatrix <- apply(t(bzbna[,-c(1)]), 1, function(row) rep(row, each = nrow(bsamp_ab)))
    bx.repmatrix <- rep(bx_ab, nrow(bsamp_bna))
    ones <- rep(1, nrow(bsamp_bna))
    ones.repmatrix <- apply(t(ones), 1, function(row) rep(row, each=nrow(bsamp_ab)))
    bZX <- cbind(ones.repmatrix, bz.repmatrix, bx.repmatrix)
    
    bzA.repmatrix <- apply(t(bz_ab[,-c(1)]), 1, function(row) rep(row, each = nrow(bsamp_ab)))
    bxA.repmatrix <- rep(bx_ab, nrow(bsamp_ab))
    onesA.repmatrix <- rep(1, nrow(bsamp_ab))
    bZX_A <- cbind(onesA.repmatrix, bzA.repmatrix, bxA.repmatrix)
    
    # Replicate matrices for kernel calculations
    bzxba.matrix <- do.call(rbind, replicate(nrow(bsamp_bna), bzxsab_all[,-c(1)], simplify=FALSE))
    bzxbna.repmatrix <- apply(t(bzxsbna_all[,-c(1)]), 1, function(row) rep(row, each = nrow(bsamp_ab)))
    bzxa.dupmatrix <- apply(t(bzxsab_all[,-c(1)]), 1, function(row) rep(row, each = nrow(bsamp_ab)))
    bzxa.repmatrix <- do.call(rbind, replicate(nrow(bsamp_ab), bzxsab_all[,-c(1)], simplify=FALSE))
    
    # Calculate bandwidth
    bwgtk <- rep(bw_sba, nrow(bsamp_bna))
    
    ### If doCalcPower is TRUE then run the code below to calculate the optimal power
    if (is.null(pwr)) {
      # If pwr is NULL, we will calculate the optimal power level
      # Create objective function for optimization
      obj_fn <- build_objective_fn(
        bzxsab_all, bwgtk, bzxba.matrix, bzxbna.repmatrix,
        bw_sbna, bsamp_ab, bZX, bxsa_reg, bysa,
        bzxa.repmatrix, bzxa.dupmatrix,
        bZX_A, bzxab_all
      )
      
      # Find optimal power level for bandwidth
      b.result <- optimize(obj_fn, interval = c(0.1, 0.5), tol = 1e-1)
      
      # Calculate theta using optimal power
      b.opt_theta <- calc_theta_pwr(
        b.result$minimum,
        bzxsab_all, bwgtk, bzxba.matrix, bzxbna.repmatrix,
        bw_sbna, bsamp_ab, bZX, bxsa_reg, bysa
      )
    } else {
      # If pwr is not null, use the power value
      b.opt_theta <- calc_theta_pwr(
        pwr,
        bzxsab_all, bwgtk, bzxba.matrix, bzxbna.repmatrix,
        bw_sbna, bsamp_ab, bZX, bxsa_reg, bysa
      )
    }
    
    # 4 - Obtain bootstrap estimate of theta
    b.hat_theta_me <- b.opt_theta$theta
    
    # 5 - Calculate bootstrap estimate of the score
    bhat_score_me <- logistic_r(bxsa_reg, b.hat_theta_me)
    
    # 6 - Calculate bootstrap estimate of Y-hat
    byest[bi] <- sum(bysa/bhat_score_me)/sum(1/bhat_score_me)
    
    if (bi %% 10 == 0) {
      print(paste("Bootstrap iteration", bi))
    }
  }
  
  # Calculate variance of bootstrap estimates
  est_var_me <- var(byest)  # Taking simple variance
  est_var_me_unb <- 1/(B - 1) * sum((byest - mean(pseudopopn$Y))^2)  # Variance estimator based on pseudo-popn mean
  
  
  return(list(
    varof_mean = est_var_me,
    variance = est_var_me_unb,
    estimates = byest
  ))
}



#####################################################

### Set up parameters for the simulation


set.seed(202506)
doBoot <- TRUE #whether to do the bootstrap variance
B=10                            #bootstrap sample size
N=20000                             #population total
iteration=50                     #simulation iteration                        
beta=c(2,1.2,1,1,1)                 #True beta for "y=x beta + e"
theta=c(0.1,0.2,0.3,0.4)          #True theta for PS, no intercept
#theta=c(0.1,0.2,0.3,0)
#Covariate
x1=rbinom(N,1,0.5) 
x2=runif(N,min=0,max=2)+0.3*x1                                                   
x3=rexp(N,1)+0.2*(x1+x2)                                                    
x4=rchisq(N,4)+0.1*(x1+x2+x3)
#x4=rchisq(N,7)+0.1*(x1)

### A measurement error version of x4 - 4 scenarios
### Un-comment the scenario of interest

x4me <-  x4 + rnorm(length(x4), 0, 1)
#x4me <- 0.8*x4 + rnorm(length(x4), 0, 1)
#x4me <- x4 + rnorm(length(x4), 0 , 2)
#x4me <- 10*(x4^(1/4)) + rnorm(length(x4), 0, 0.7)
#plot(x4,x4me)

#Design matrix
X=cbind(rep(1,N),x1,x2,x3,x4)
#working matrices
X_reg=cbind(rep(1,N),x1,x2,x3,x4)
X_res=cbind(rep(1,N),x1,x2,x3,x4)
e=rnorm(N,mean=0,sd=1)

### Measurement error versions of the X matrix
X_reg_me=cbind(rep(1,N),x1,x2,x3,x4me)
X_res_me=cbind(rep(1,N),x1,x2,x3,x4me)



rho=0.3             #cor(Y,XB)
mean_n_SA=2000   #sample size of S_A 
mean_n_SB=1000  #sample size of S_B 


#find sigma
sigma=sqrt((1/rho^2-1)*var(c(X%*%beta)))

#Generate y
Y=X%*%beta+e*sigma
#Check cor(Y,XB)
cor(Y,X%*%beta)
mu=mean(Y)


### We will now put Y and X into a data frame
popn <- data.frame(Y,X,x4me)

### Add a column with a unit identifier
popn$id <- 1:nrow(popn)

###############################################

#find intercept

etap=as.vector(X[,-1]%*%theta)
dif=1
L=-10
R=-1
while(dif>0){
  M=(L+R)/2
  ps=exp(M+etap)/(1+exp(M+etap))
  if(sum(ps)<=mean_n_SA) L=M
  if(sum(ps)>=mean_n_SA) R=M
  if(abs(sum(ps)-mean_n_SA)<=0.5) dif=0
}


#true score
score=exp(M+etap)/(1+exp(M+etap))   
(min(score))
(max(score))
(sum(score))

#set max(z1)/min(z1)=50
#find c

z=x3
c=(max(z)-20*min(z))/19
s=sum(z+c)
w=s/(z+c)/mean_n_SB
z1=(z+c)/s


### Now calculate probability of selection into S_B (probability sample)
cost.i <- 1

### Form final population level dataset including the design weights for probability sample B
popn <- cbind(popn,w)


#################################################
### Start simulation iterations
#################################################


Bias=0
MSE=0
est_var=0
c_rate=0


simresults <- foreach(g = 1:iteration, .combine=rbind) %dorng% {
  
  #draw the non-probability sample 
  sam=rbinom(N,1,score)
  n_SA=sum(sam)
  index_SA=order(-sam)[1:n_SA] 
  Y_SA=Y[index_SA]
  X_SA_reg=X_reg[index_SA,]
  X_SA_res=X_res[index_SA,]
  w_SA=w[index_SA]
  
  samp_A <- popn[index_SA,]
  
  
  ### Draw probability sample using Pareto sampling
  pik <- sampling::inclusionprobabilities(z1, mean_n_SB)
  omega <- runif(N)
  lambda <- pik
  index_SB = order(omega * (1 - lambda)/((1 - omega) * lambda))[1:mean_n_SB]
  
  n_SB=mean_n_SB
  X_SB_res=X_res[index_SB,]
  X_SB_reg=X_reg[index_SB,]
  w_SB=w[index_SB]
  Y_SB <- Y[index_SB]
  
  ### Measurement Error versions
  X_SB_res_me=X_res_me[index_SB,]
  X_SB_reg_me=X_reg_me[index_SB,]
  
  samp_B <- popn[index_SB,]
  
  
  ############################################################
  
  
  ### Obtain the subsets which we will use to correct selection bias
  
  ### We assume that the x values are available from a sample
  
  samp_B_not_A <- samp_B[!(samp_B$id %in% samp_A$id),] # Sample B \ A
  popn_not_A <- popn[!(popn$id %in% samp_A$id),] # Segment C in the population
  samp_AB <- samp_B[samp_B$id %in% samp_A$id,] # Sample B and A
  
  
  ### Covariate datasets assuming a MAR/SAR setup
  X_SBnA_MAR_res <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id), c("V1","x1","x2","x3","x4me")])
  X_SBnA_MAR_reg <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id), c("V1","x1","x2","x3","x4me")])
  w_SBnA_MAR <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id),c("w")])
  
  ### Covariate datasets assuming true values available for all x variables
  X_SBnA_res <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  X_SBnA_reg <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id), c("V1","x1","x2","x3","x4")])
  w_SBnA_nome <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id),c("w")])
  
  ### Now put together the various datasets we need for B and A
  X_SBpA_res <- as.matrix(samp_B[samp_B$id %in% samp_A$id, c("V1","x1","x2","x3","x4")])
  X_SBpA_reg <- as.matrix(samp_B[samp_B$id %in% samp_A$id, c("V1","x1","x2","x3","x4")])
  w_SBpA <- as.matrix(samp_B[samp_B$id %in% samp_A$id,c("w")])
  
  ### The datasets assuming the sample B is available
  ### Use X here to represent the measurement error prone variables
  ### Use Z here to represent the variables that do not suffer from measurement error
  
  X_AB <- as.matrix(samp_B[samp_B$id %in% samp_A$id, c("x4")])
  ZXs_AB <- as.matrix(samp_B[samp_B$id %in% samp_A$id, c("V1","x1","x2","x3","x4me")])
  
  Z_BnA <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id), c("V1","x1","x2","x3")])
  ZXs_BnA <- as.matrix(samp_B[!(samp_B$id %in% samp_A$id), c("V1","x1","x2","x3","x4me")])
  
  Z_AB <- as.matrix(samp_B[samp_B$id %in% samp_A$id, c("V1","x1","x2","x3")])
  ZXs_AB <- as.matrix(samp_B[samp_B$id %in% samp_A$id, c("V1","x1","x2","x3","x4me")])
  ZX_AB <- as.matrix(samp_B[samp_B$id %in% samp_A$id, c("V1","x1","x2","x3","x4")])
  
  w_SBnA <- samp_B_not_A$w
  w_SBA <- samp_AB$w
  
  
  ###############################################
  
  ### Create the datasets we need, which require the ij combinations for the data
  
  ### First - assuming population x available
  
  ### This code results in a matrix to feed into expected propensities ie over the
  ### double sum A\B and A and B. ZX contains Z variables over A\B, for each row of Z
  ### there is a full copy of the X variables from A and B
  z.repmatrix <- apply(t(Z_BnA[,-c(1)]), 1, function(row) rep(row, each = nrow(samp_AB)))
  x.repmatrix <- rep(X_AB, nrow(samp_B_not_A))
  ones <- rep(1, nrow(samp_B_not_A))
  ones.repmatrix <- apply(t(ones),1,function(row) rep(row, each=nrow(samp_AB)))
  ZX <- cbind(ones.repmatrix,z.repmatrix,x.repmatrix)
  
  
  ### Also get ZX for the units in AB
  zAB.repmatrix <- apply(t(Z_AB[,-c(1)]), 1, function(row) rep(row, each = nrow(samp_AB)))
  xAB.repmatrix <- rep(X_AB, nrow(samp_AB))
  onesA.repmatrix <- rep(1, nrow(samp_AB))
  ZX_AB_matrix <- cbind(onesA.repmatrix, zAB.repmatrix, xAB.repmatrix)
  
  
  ### Create matrices that will feed into the kernel estimation
  ### These matrices contain Z and X-star (measurement error versions of X)
  ZXs_BA.matrix <- do.call(rbind, replicate(nrow(samp_B_not_A), ZXs_AB[,-c(1)], simplify=FALSE))
  ZXs_BnA.repmatrix <- apply(t(ZXs_BnA[,-c(1)]), 1, function(row) rep(row, each = nrow(samp_AB)))
  
  ### These matrices will feed into kernel estimation
  ### They contain Z and X-star over A and B
  ZXs_AB.repmatrix <- apply(t(ZXs_AB[,-c(1)]), 1, function(row) rep(row, each = nrow(samp_AB)))
  ZXs_AB.dupmatrix <- do.call(rbind, replicate(nrow(samp_AB), ZXs_AB[,-c(1)], simplify=FALSE))
  
  
  #################################################################
  
  ### Select a bandwidth and use it to calculate theta
  
  wgtk <- rep(w_SBA, nrow(samp_B_not_A))
  xsa <- X_SA_reg
  
  ### Set up the function to calculate the objective function
  obj_fn <- build_objective_fn(
    ZXs_AB, wgtk, ZXs_BA.matrix, ZXs_BnA.repmatrix,
    w_SBnA, samp_AB, ZX, xsa, Y_SA,
    ZXs_AB.dupmatrix, ZXs_AB.repmatrix,
    ZX_AB_matrix, ZX_AB
  )
  
  ### Find the optimal power level for the bandwidth
  result <- optimize(obj_fn, interval = c(0.1, 0.5), tol = 1e-1)
  
  ### Calculate theta using the optimal power
  opt_theta <- calc_theta_pwr(result$minimum,
                              ZXs_AB, wgtk, ZXs_BA.matrix, ZXs_BnA.repmatrix,
                              w_SBnA, samp_AB, ZX, xsa, Y_SA
  ) 
  
  
  ### Extract the theta estimate
  hat_theta_me <- opt_theta$theta
  
  #########################################
  
  ### Create a model to estimate x4 using A and B sample
  ### This is a regression model for the measurement error correction
  ### Use samp_AB to create a model for x4me ~ x4
  
  me.model <- glm(x4me ~ x4, data = samp_AB, family = gaussian(link = "identity"))
  
  ### Now apply the model to estimate x4 in samp_B_not_A
  samp_B_not_A$pred_x4 <- (samp_B_not_A$x4me - me.model$coefficients[1])/me.model$coefficients[2]
  X_SBnA_est <- as.matrix(samp_B_not_A[, c("V1","x1","x2","x3","pred_x4")])
  
  
  #########################################
  
  ### Produce estimates for measurement error correction approach
  hat_score_ME_SA <- logistic_r(X_SA_reg, hat_theta_me)
  
  ### Produce estimates using the model to correct
  hat_model_SA <- hat_theta2_f(X_SA_res, X_SBnA_est, w_SBnA_MAR)
  hat_score_model_SA <- 1/(1 + exp(-X_SA_res %*% hat_model_SA))
  
  
  #estimated theta and scores for CLW
  hat_theta=hat_theta_f(X_SA_res,X_SB_res,w_SB)
  hat_score_SA=1/(1+exp(-X_SA_res%*%hat_theta))
  hat_score_SB=1/(1+exp(-X_SB_res_me%*%hat_theta))
  
  ### Estimated theta for KW method
  hat_theta3 <- hat_theta_KW_f(X_SBpA_res, X_SBnA_MAR_res, w_SBpA, w_SBnA_MAR)
  hat_score3_SA <- 1/(1+exp(-X_SA_res%*%hat_theta3))
  hat_score3_SB <- 1/(1+exp(-X_SB_res%*%hat_theta3))
  
  ### Estimated theta for KW method - with model to correct
  hat_theta3_model <- hat_theta_KW_f(X_SBpA_res, X_SBnA_est, w_SBpA, w_SBnA_MAR)
  hat_score3_model_SA <- 1/(1+exp(-X_SA_res%*%hat_theta3_model))
  
  
  ### Estimated theta for the alternative method
  hat_theta2 <- hat_theta2_f(X_SA_res, X_SBnA_MAR_res, w_SBnA_MAR)
  hat_score2_SA=1/(1+exp(-X_SA_res%*%hat_theta2))
  hat_score2_SB=1/(1+exp(-X_SB_res%*%hat_theta2))
  
  hat_theta2_nome <- hat_theta2_f(X_SA_res,X_SBnA_res,w_SBnA_nome)
  hat_score2_nome_SA=1/(1+exp(-X_SA_res%*%hat_theta2_nome))
  hat_score2_nome_SB=1/(1+exp(-X_SB_res%*%hat_theta2_nome))
  
  ### Produce estimates of population size
  
  #check sum
  N_B=sum(w_SB)
  
  N_A=sum(1/hat_score_SA)
  N_A2=sum(1/hat_score2_SA)
  N_A_me <- sum(1/hat_score_ME_SA)
  N_A_model <- sum(1/hat_score_model_SA)
  N_A_nome <- sum(1/hat_score2_nome_SA)
  
  N_A_KW <- sum(1/hat_score3_SA)
  NA_KW_model <- sum(1/hat_score3_model_SA)
  
  #################################################
  ### Point estimators
  #################################################
  
  
  mu_naive=mean(Y_SA)
  
  mu_ipw2=sum(Y_SA/hat_score_SA)/N_A
  
  mu_kw <- sum(Y_SA/hat_score3_SA)/N_A_KW
  
  mu_kw_model <- sum(Y_SA/hat_score3_model_SA)/NA_KW_model
  
  mu_new <- sum(Y_SA/hat_score2_SA)/N_A2
  
  mu_new_nome <- sum(Y_SA/hat_score2_nome_SA)/N_A_nome
  
  mu_ME <- sum(Y_SA/hat_score_ME_SA)/N_A_me
  
  mu_model <- sum(Y_SA/hat_score_model_SA)/N_A_model
  
  ### Create estimates of mean
  muhat=c(mu_naive, mu_kw, mu_new, mu_kw_model, mu_model, mu_ME, mu_new_nome)
  
  
  #################################################
  ### Variance calculations
  #################################################
  
  a=t(X_SB_res)%*%(c(w_SB*hat_score_SB*(1-hat_score_SB))*X_SB_res)
  anew <- t(X_SB_res)%*%(c(w_SB*hat_score2_SB*(1-hat_score2_SB))*X_SB_res)
  anew2 <- t(X_SA_res)%*%(c(1-hat_score2_SA)*X_SA_res)
  
  hat_b1=solve(a,t(X_SA_res)%*%c((1/hat_score_SA-1)*Y_SA))
  hat_b2=solve(a,t(X_SA_res)%*%c((1/hat_score_SA-1)*(Y_SA-mu_ipw2)))

  hat_bnew <- solve(anew2,t(X_SA_res)%*%c((1/hat_score2_SA-1)*(Y_SA-mu_new)))
  
  ### Variance component for IPW2 in CLW
  tt=c(X_SB_res%*%hat_b2)*c(hat_score_SB)
  st=sum(w_SB*tt)
  est_vd2=1/N^2*sum((w_SB*tt-st/n_SB)^2)
  
  ### variance component for new method
  tt=c(X_SB_res%*%hat_bnew)*(c(hat_score2_SB)*(1-c(hat_score2_SB)))
  st=sum(w_SB*tt)
  est_vdnew=1/N^2*sum((w_SB*tt-st/n_SB)^2)
  
  est_vr1=(1/N^2)*sum((1-hat_score_SA)*(Y_SA/hat_score_SA-X_SA_res%*%hat_b1)^2)
  est_vr2=(1/N^2)*sum((1-hat_score_SA)*((Y_SA-mu_ipw2)/hat_score_SA-X_SA_res%*%hat_b2)^2)

  ### estimated variance for new method
  est_vr2new <- (1/N^2)*sum((1-hat_score2_SA)*((Y_SA-mu_new)/hat_score2_SA-X_SA_res%*%hat_bnew)^2)
  
  est_vrnew <- (1/N^2) * sum((hat_score2_SA)^2*(1 - hat_score2_SA)*(w_SA - 1)*(X_SA_res%*%hat_bnew)^2)
  est_var_new <- est_vr2new + est_vdnew + est_vrnew
  
  
  
  ### estimated variances for IPW2 in CLW
  est_var_ipw2=est_vd2+est_vr2
  
  ### See if we want to do bootstrap variances
  
  if (doBoot) {
    ### Create a model for x4me using samp_AB
    mod.me <- glm(x4me ~ x4, data = samp_AB, family = gaussian(link = "identity"))
    ### What is the variance of the error term
    summary(mod.me)$dispersion
    ### Now apply the model to samp_A to get predicted x4me values
    samp_A$pred_x4me <- mod.me$coefficients[1] + mod.me$coefficients[2]*samp_A$x4 + rnorm(nrow(samp_A), 0, sqrt(summary(mod.me)$dispersion))
    
    ### Now create variable x4me.boot which takes value x4 if in samp_AB and predicted value if in samp_A only
    samp_A$x4me_boot <- ifelse(samp_A$id %in% samp_AB$id,
                               samp_A$x4me,
                               samp_A$pred_x4me)
    
    ### estimated variance using bootstrap for the non-parametric method
    boot_results <- boot_var_me(samp_A, X_SA_reg, Y_SA, hat_score_ME_SA, 
                                N, mean_n_SB, B=B,pwr=result$minimum)
    
    ### Extract the bootstrap variance
    est_var_me <- boot_results$variance
  } else {
    est_var_me <- 1
  }
  
  #estimated variance
  #est_var=est_var+c(est_var_ipw1,est_var_ipw2,est_var_new,est_var_dr1,est_var_dr2,boot_var_dr2,est_var_hk)/iteration
  
  est_var=c(est_var_ipw2,est_var_new, est_var_me)
  
  
  #sd
  sd_ipw2=sqrt(est_var_ipw2)
  sd_new <- sqrt(est_var_new)
  sd_me <- sqrt(est_var_me)
  

  
  #CI
  CI_ipw2=c(mu_ipw2-1.96*sd_ipw2,mu_ipw2+1.96*sd_ipw2)
  CI_new <- c(mu_new - 1.96*sd_new, mu_new + 1.96*sd_new)
  CI_me <- c(mu_ME - 1.96*sd_me, mu_ME + 1.96*sd_me)
  
  ### Average length of the confidence intervals
  al_ipw2=CI_ipw2[2]-CI_ipw2[1]
  al_new=CI_new[2]-CI_new[1]
  al_me=CI_me[2]-CI_me[1]
  
  avelength <- c(al_ipw2, al_new, al_me)
  
  
  #cp
  
  CI_ipw2_c=(mu>=CI_ipw2[1])*(mu<=CI_ipw2[2])
  CI_new_c <- (mu >= CI_new[1])*(mu <= CI_new[2])
  CI_me_c <- (mu >= CI_me[1])*(mu <= CI_me[2])

  
  c_rate=c(CI_ipw2_c, CI_new_c, CI_me_c)
  
  
  if(floor(g/10)==g/10) print(g)
  
  ### Create list of results to output
  list(ests=muhat,est_var=est_var,c_rate=c_rate,avelength=avelength)
  
}

endtime <- Sys.time()
endtime-starttime


##########################################

### Collate results and produce monte carlo estimates of
### relative bias, MSE, coverage rate and average length
### of confidence intervals

sim.ests <- matrix(unlist(simresults[,1]),byrow=T,nrow=iteration)
sim.vars <- matrix(unlist(simresults[,2]), byrow=T,nrow = iteration)
sim.c <- matrix(unlist(simresults[,3]), byrow=T, nrow=iteration)
sim.al <- matrix(unlist(simresults[,4]), byrow=T, nrow=iteration)

### Produce the column means of sim.ests and then bias
mean.ests <- colMeans(sim.ests)
Bias <- mean.ests - mu
RB <- 100*Bias/mu

### Produce estimated MSE using variance methods
est.var <- colMeans(sim.vars)
est.RSE <- 100*sqrt(est.var)/mu

### Produce MC MSE
MC.MSE <- apply(sim.ests,2,var) + Bias^2
MC.RRMSE <- 100*sqrt(MC.MSE)/mu
MC.variance <- MC.MSE - Bias^2
MC.RSE <- 100*sqrt(MC.variance)/mu


### Produce coverage rate
c_rate <- colMeans(sim.c)

### Average length of CIs
al <- colMeans(sim.al)

################################################

### Output results to a data frame and then create latex table

library(xtable)
output1=data.frame("Estimator"=c("mu_naive","mu_kw","mu_new","mu_kw_model","mu_model", "mu_ME","mu_new_nome"),
                   "RB"=RB,"RRMSE"=MC.RRMSE)
rownames(output1) <- NULL

table_title=paste0("na=",mean_n_SA,",","nb=",mean_n_SB,", ","cor(Y,XB)=",rho*100)
print(xtable(output1, caption = table_title,digits=3),
      caption.placement = 'top',include.colnames =T, include.rownames = F)


