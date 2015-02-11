basic_info=function(x){
  print("BMI_con")
  print(summary(x$BMI[x$AFFSTAT==1]))
  print("BMI_case")
  print(summary(x$BMI[x$AFFSTAT==2]))
  print("AGECOL_con")
  print(summary(x$AGECOL[x$AFFSTAT==1]))
  print("AGEATK_case")
  print(summary(x$AGE1ATK[x$AFFSTAT==2]))
  print("SEXFULL")
  print(table(x$SEXFULL,x$AFFSTAT),exclude=NULL)
  print("SEXCOVAR")
  print(table(x$SEXCOVAR,x$AFFSTAT),exclude=NULL)
  print("RS475414")
  print(table(x$RS475414,x$AFFSTAT,exclude=NULL))
  print("RS17299124")
  print(table(x$RS17299124,x$AFFSTAT,exclude=NULL))
  print("RS693591")
  print(table(x$RS693591,x$AFFSTAT,exclude=NULL))
  print("RS17300741")
  print(table(x$RS17300741,x$AFFSTAT,exclude=NULL))
  print("RS2078267")
  print(table(x$RS2078267,x$AFFSTAT,exclude=NULL))
  print("RS3825018")
  print(table(x$RS3825108,x$AFFSTAT,exclude=NULL))
  print("RS475688")
  print(table(x$RS475688,x$AFFSTAT,exclude=NULL))
  print("RS7932775")
  print(table(x$RS7932775,x$AFFSTAT,exclude=NULL))
  print("RS476037")
  print(table(x$RS476037,x$AFFSTAT,exclude=NULL))
  print("RS478607")
  print(table(x$RS478607,x$AFFSTAT,exclude=NULL))
  print("RS12289836")
  print(table(x$RS12289836,x$AFFSTAT,exclude=NULL))
  print("RS642803")
  print(table(x$RS642803,x$AFFSTAT,exclude=NULL))
}

for i in 1060{
  if(CAU[i,i]=="."){
    CAU[i,i] <- NA
  }
}

basic_info2=function(x){
print("BMI_con")
print(summary(x$BMI[x$AFFSTAT==1]))
print("BMI_case")
print(summary(x$BMI[x$AFFSTAT==2]))
print("AGECOL_con")
print(summary(x$AGECOL[x$AFFSTAT==1]))
print("AGEATK_case")
print(summary(x$AGE1ATK[x$AFFSTAT==2]))
}


basic_info3=function(x){
print("BMI_con")
print(sd(x$BMI[x$AFFSTAT==1],na.rm=T))
print("BMI_case")
print(sd(x$BMI[x$AFFSTAT==2],na.rm=T))
print("AGECOL_con")
print(sd(x$AGECOL.x[x$AFFSTAT==1],na.rm=T))
print("AGEATK_case")
print(sd(x$AGE1ATK[x$AFFSTAT==2],na.rm=T))
}
