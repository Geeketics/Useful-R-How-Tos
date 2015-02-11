# Tanya Flynn
# October 2013
# University of Otago
#

ES_Compare=function(EFF1=NULL, SE1=NULL, EFF2=NULL, SE2=NULL, type="beta"){
  if(type=="OR"){
    EFF1=log(EFF1)
    EFF2=log(EFF2)
  }
  ##calculating difference, standard error, p-value and confidence intervals
  diff=EFF1-EFF2
  diff.se=sqrt((SE1^2)+(SE2^2))
  diff.lci=diff-(1.96*diff.se)
  diff.uci=diff+(1.96*diff.se)
  diff.z=diff/diff.se
  diff.p=pnorm(abs(diff.z),lower.tail=FALSE)*2
  ##calculating values for Q
  #weights
  EFF1.w=1/(SE1^2)
  EFF2.w=1/(SE2^2)
  EFF.w.sum=EFF1.w+EFF2.w
  #weight*effect size
  EFF1.w.es=EFF1.w*EFF1
  EFF2.w.es=EFF2.w*EFF2
  EFF.w.es.sum=EFF1.w.es+EFF2.w.es
  #weight*effect size squared
  EFF1.w.es.sqr=EFF1.w*(EFF1^2)
  EFF2.w.es.sqr=EFF2.w*(EFF2^2)
  EFF.w.es.sqr.sum=EFF1.w.es.sqr+EFF2.w.es.sqr
  #calculating Q and het-p (confirms difference)
  Q=EFF.w.es.sqr.sum-((EFF.w.es.sum^2)/EFF.w.sum)
  k=2
  het.p=1-pchisq(Q,(k-1))
  #output file
  output <- c(EFF1=EFF1,EFF2=EFF2,SE1=SE1,SE2=SE2,diff=diff,diff.se=diff.se,diff.lci=diff.lci,diff.uci=diff.uci,diff.z=diff.z,diff.p=diff.p,Q=Q,df=k-1,het.p=het.p)
}

show=function(list){
  print("Input Values")
  print(c(round(list[1],4),round(list[3],4)))
  print(c(round(list[2],4),round(list[4],4)))
  print("Difference")
  print(c(round(list[5],4),round(list[6],4),round(list[7],4),round(list[8],4)))
  print("Significance")
  print(c((round(list[9],4)),list[10]))
  print("Heterogeneity")
  print(c(round(list[11],4),round(list[12],1),list[13]))
}