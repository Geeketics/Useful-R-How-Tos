########
# Tanya Flynn
# Nov 2013
# University of Otago
#

#if using multiple covariates input as eg: Cov=list(Cov1=ARIC$SEX,Cov2=ARIC$AGE,Cov3=ARIC$BMI)
Formal_GxE=function(SNP=NULL,Env=NULL,DepVar=NULL,Cov=NULL){
  #make two new variables
  SNP_2 <- ifelse(SNP==2,1,0)
  SNP_3 <- ifelse(SNP==3,1,0)
  
  #make two interaction terms
  INT_2 <- SNP_2*Env
  INT_3 <- SNP_3*Env
  
  #make new variables matrix
  data <- data.frame(DepVar=DepVar,SNP=as.factor(SNP),Env=Env,INT_2=INT_2,INT_3=INT_3)
  
  #merge data and Cov
  data <- merge(data,Cov,by="row.names")
  data$Row.names <- NULL

  #I think this is necessary
  str(data)
  
  #linear regression including new terms
  fit <- lm(DepVar~.,data=data)
  
  #show results
  print(summary(fit))
}