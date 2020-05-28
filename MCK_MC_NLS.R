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
# Function for nls
predictNLS.sum <- function(
  object, 
  newdata,
  nsim,
  x.var="LakeArea_m2",
  ...
)
{
#  require(MASS, quietly = TRUE)
  
  ## get right-hand side of formula
  RHS <- as.list(object$call$formula)[[3]]
  EXPR <- as.expression(RHS)
  
  ## all variables in model
  VARS <- all.vars(EXPR)
  
  ## get parameter coefficients
  COEF <- coef(object)
  
  # extract predictor variable    
  predNAME <- setdiff(VARS, names(COEF))  
  
  ## take fitted values, if 'newdata' is missing
  if (missing(newdata)) {
    newdata <- eval(object$data)[predNAME]
    colnames(newdata) <- predNAME
  }
  
  
  ## get variance-covariance matrix
  VCOV <- vcov(object)
  
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
  
  # create dataframe for output, specific to glacial lake volume
  outMAT<-matrix(NA,ncol=nsim,nrow=nrow(newdata))
  outPred.df <- data.frame(Year_Start=newdata$Year_Start,LakeArea_m2=newdata$LakeArea_m2)
  outPred.df<-cbind(outPred.df,outMAT)
  
  
  for (i in 1:NR) {
    counter(i)
    
    ## get predictor values and optional errors
    predVAL <- newdata[i, x.var]
    names(predVAL) <- predNAME  

    ## create mean vector for 'mvrnorm'
   K2<-rnorm(nsim,COEF[2],sqrt(VCOV[2,2]))
    
    
    ## create MC simulation matrix
    #    simMAT <- mvrnorm(n = nsim, mu = MU, Sigma = newVCOV, empirical = TRUE)
    simMAT<-cbind(k1=rep(COEF[1],nsim),k2=K2,LakeArea_m2=rep(predVAL,nsim))

    ## evaluate expression on rows of simMAT
    EVAL <- try(eval(EXPR, envir = as.data.frame(simMAT)), silent = TRUE)
    if (inherits(EVAL, "try-error")) stop("There was an error evaluating the simulations!")
    
    outPred.df[i,3:ncol(outPred.df)]<-EVAL
    ## collect statistics
  }

  cat("\n")
  
  return(list(predictions=outPred.df,OneSim=simMAT))  
}
