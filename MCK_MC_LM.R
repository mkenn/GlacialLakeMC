# ###########
# Monte Carlo analysis of predicted lake volume
# For Shugar et al. Rapid worldwide growth of glacial lakes since 1990. Nature Climate Change
# script adapted from https://www.r-bloggers.com/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/
# web site last accessed May 27, 2020
# The goal of this script is to obtain nsim plausible volumes for each lake area
# The sum of each of these across all lakes in a given time frame will then be 
# a plausible total lake volume. The distribution of these totals will be 
# a plausible distribution of total lake volume
# This is in the spirit of Monte Carlo prediction intervals for model outputs
# across a reasonable characterization of parameter uncertainty
############
# Function for linear model
predictLM.sum <- function(
  object, 
  newdata,
  nsim,
  x.var="area.log",
  ...
)
{
  require(MASS, quietly = TRUE)
  
  ## get right-hand side of formula
  RHS <- as.list(object$call$formula)[[3]]
  EXPR <- as.expression(RHS)
  
  ## all variables in model
  VARS <- all.vars(EXPR)
  
  ## coefficients
  COEF <- coef(object)
  
  ## extract predictor variable    
  predNAME <- setdiff(VARS, names(COEF))  
  
  ## take fitted values, if 'newdata' is missing
  if (missing(newdata)) {
    newdata <- eval(object$data)[predNAME]
    colnames(newdata) <- predNAME
  }
  
  ## get variance-covariance matrix
  VCOV <- vcov(object)
  
  ## augment variance-covariance matrix for 'mvrnorm' 
  ## by adding a column/row for 'error in x'
  # Here we use 0 for all
  NCOL <- ncol(VCOV)
  ADD1 <- c(rep(0, NCOL))
  ADD1 <- matrix(ADD1, ncol = 1)
  colnames(ADD1) <- predNAME
  VCOV <- cbind(VCOV, ADD1)
  ADD2 <- c(rep(0, NCOL + 1))
  ADD2 <- matrix(ADD2, nrow = 1)
  rownames(ADD2) <- predNAME
  VCOV <- rbind(VCOV, ADD2) 
  
  ## iterate over all entries in 'newdata' as in usual 'predict.' functions
  NR <- nrow(newdata)
  respVEC <- numeric(NR)
  seVEC <- numeric(NR)
  varPLACE <- ncol(VCOV)   
  
  ## define counter function
  counter <- function (i) 
  {
    if (i%%10 == 0) 
      cat(i)
    else cat(".")
    if (i%%50 == 0) 
      cat("\n")
    flush.console()
  }
  
  outMAT<-matrix(NA,ncol=nsim,nrow=nrow(newdata))
  outPred.df <- data.frame(Year_Start=newdata$Year_Start,LakeArea_m2=newdata$LakeArea_m2)
  outPred.df<-cbind(outPred.df,outMAT)
  
  
  for (i in 1:NR) {
    counter(i)
    
    ## get predictor values and optional errors
    predVAL <- newdata[i, x.var]
    predERROR<-0
    # if (ncol(newdata) == 2) predERROR <- newdata[i, 2] else predERROR <- 0
    names(predVAL) <- predNAME  
    names(predERROR) <- predNAME  
    
    ## create mean vector for 'mvrnorm'
    #   K2<-rnorm(nsim,COEF[2],VCOV[2,2])
    
    
    MU <- c(COEF, predVAL)
    # 
    # ## create variance-covariance matrix for 'mvrnorm'
    # ## by putting error^2 in lower-right position of VCOV
    newVCOV <- VCOV
    newVCOV[varPLACE, varPLACE] <- predERROR^2
    # 
    ## create MC simulation matrix
    simMAT <- rmvnorm(n = nsim, mean = MU, sigma = newVCOV)#, empirical = TRUE)
    
# Predict volume for log-log linear model
    outPred.df[i,3:ncol(outPred.df)]<-exp(simMAT[,1]+simMAT[,2]*simMAT[,3])
  }
  return(list(predictions=outPred.df,OneSim=simMAT))  
}
