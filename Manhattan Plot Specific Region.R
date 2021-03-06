###########
# Tanya Flynn
# Nov 2013
# Manhattan Plot for Specific Region

METAL_Plot=function(BP=NULL,P=NULL,SNP=NULL,Gene=NULL,Genestart=NULL,Geneend=NULL,Title=NULL){

#convert P values to log P values
log.P <- -log10℗
#multiple testing threshold
threshold <- -log10(0.05/length(P))
#plot P values into manhattan plot
plot(BP,log.P,main=as.character(Title),xlim=c(min(BP-10),max(BP+10)), ylim=c(-0.5,max(c(threshold+0.5,max(log.P)),
     xlab="Base Position",ylab=expression(-log[10](italic(p))),las=1,pch=20,col="grey50")
#put x-axis back in
abline(a=0,b=0,lwd=1.2)
#add multiple testing threshold
abline(a=threshold,b=0,col="red")
#add suggestive threshold
abline(a=(-log10(0.05)),b=0,col="blue")
#annotate top hits
text(x=(BP[log.P>(-log10(0.05/length(P)))]),y=(log.P[log.P>(-log10(0.05/length(P)))]+0.5),labels=as.character(SNP[log.P>(-log10(0.05/length(P)))]),srt=90,cex=0.7)
#annotate gene positions
lines(x=c(genestart,geneend),y=c(-0.2,-0.2),col="black",lwd=2)
text(x=(genestart+geneend)/2,y=-0.4,labels=as.character(Gene))
}