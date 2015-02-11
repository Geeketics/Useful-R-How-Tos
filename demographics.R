demographics=function(x){
  for(i in 1:length(x[1,])){
    print(names(x[i]))
    print(c("mean",mean(x[,i],na.rm=T)))
    print(c("sd",sd(x[,i],na.rm=T)))
  }
}