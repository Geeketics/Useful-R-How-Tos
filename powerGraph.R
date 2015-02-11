# Tanya Flynn
# August 2013
# University of Otago
#
powerOR=function(caseN=NULL,controlN=NULL,controlA1=NULL,oddsratio=NULL,siglevel=NULL){
  
  #warnings
  if(sum(sapply(list(caseN,controlN,controlA1,oddsratio),is.null))>0)
    stop("caseN, controlN, controlA1, and oddsratio require values")
  if(!is.null(oddsratio) && oddsratio < 0)
    stop("oddsratio must be positive: suggest using 1.2 (weak), 1.5 (moderate) or 2.0 (strong)")
  if((!is.null(caseN) && caseN <1) || (!is.null(controlN) && controlN < 1))
    stop("need at least 1 case and 1 control")
  if(!is.null(controlA1) && any(controlA1 < 0 | controlA1 >1))
    stop("controlA1 must be between 0 and 1")
  if(!is.null(siglevel) && any(siglevel <0 | siglevel >1))
    stop("siglevel must be between 0 and 1")
  
  #making required values
  if(is.null(siglevel))
    siglevel=0.05
  else siglevel=siglevel
  mfactor <- controlN/caseN
  controlA2 <- 1-controlA1
  caseA1 <- (controlA1*oddsratio)/(1+(controlA1*oddsratio)-controlA1)
  caseA2 <- 1-caseA1
  caseChr <- caseN*2
  logodds <- log(oddsratio)
  
  if(siglevel==0.05)
    criticalvalue <- 1.96
  else if(siglevel==0.01)
    criticalvalue <- 2.58
  else if(siglevel==0.001)
    criticalvalue <- 3.29
  else if(siglevel==0.0001)
    criticalvalue <- 3.89
  else if(siglevel==0.00001)
    criticalvalue <- 4.27
  else if(siglevel==0.000001)
    criticalvalue <- 4.75
  else if(siglevel==0.0000001)
    criticalvalue <- 5.33
  else if(siglevel==0.00000001)
    criticalvalue <- 5.73
  else if(siglevel==0.00000005)
    criticalvalue <- 5.45
  else stop("significance level unknown")
  
  #power calculation
  numerator <- criticalvalue*(sqrt((1/(caseChr*controlA1))+(1/(mfactor*caseChr*controlA1))+(1/(caseChr*controlA2))+(1/(mfactor*caseChr*controlA2))))-logodds
  denominator <- sqrt((1/(caseChr*caseA1))+(1/(mfactor*caseChr*controlA1))+(1/(caseChr*caseA2))+(1/(mfactor*caseChr*controlA2)))
  zscore <- numerator/denominator
  power <- 1-pnorm(zscore)
  
  #output file
  note <- "based on the calculations published by Johnson et al. (2001) Haplotype tagging for the identification of common disease genes. Nature Genetics vol29 pg233-237"
  structure(list(controlN=controlN,caseN=caseN,controlA1=controlA1,controlA2=controlA2,caseA1=caseA1,caseA2=caseA2,oddsratio=oddsratio,siglevel=siglevel,criticalvalue=criticalvalue,power=power,note=note),class="list")
}

powerGraph=function(caseN=NULL,controlN=NULL,controlA1=NULL,oddsratio=NULL,siglevel=NULL,title=NULL){
  if(is.null(oddsratio))
    oddsratio <- c(1.2,1.5,2.0)
  else oddsratio=oddsratio  
  if(is.null(controlA1))
    controlA1 <- seq(0.01,0.5,0.005)
  else controlA1=controlA1
  
  #specify how many variable to use  
  noddsratio <- length(oddsratio)
  ncontrolA1 <- length(controlA1)
  #make file of graphable variables
  result <- data.frame(matrix(nrow=(ncontrolA1),ncol=(noddsratio+1)))
  names(result)=c("allelefreq",(oddsratio))
  result$allelefreq <- controlA1
  
  #calculate detection power
  for(i in 1:ncontrolA1){
    for(j in 1:noddsratio){
      output <- powerOR(caseN=caseN, controlN=controlN, controlA1=controlA1[i],oddsratio=oddsratio[j],siglevel=siglevel)
      result[i,(j+1)] <- output[10]
    }
  }
  
  # set up graph
  xrange <- range(controlA1)
  yrange <- 0:1
  plot(xrange, yrange, type="n",
       xlab="Allele Frequency",
       ylab="Detection Power" )
  
  # add power curves
  for (i in 1:noddsratio){
    lines(controlA1, result[,(i+1)], type="l", lwd=2, lty=i)
  }
  
  # add legend
  title(as.character(title))
  legend("topright", title="Odds Ratio", as.character(oddsratio), lty=1:noddsratio)
  
}