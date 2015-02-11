library("MASS")
library("matrixcalc")
library("sem")
library("aod")

output_OS <- data.frame(matrix(nrow=18, ncol=11))
names(output_OS) <- c("group","F","R2","ols_B","ols_se","ols_p","tsls_B","tsls_se","tsls_p","DH_F","DH_p")
output_OS[,1] <- c("All_Crude_Total", "All_Adj_Total","All_Crude_Neck","All_Adj_Neck","All_Crude_Spine","All_Adj_Spine","Male_Crude_Total", "Male_Adj_Total","Male_Crude_Neck","Male_Adj_Neck","Male_Crude_Spine","Male_Adj_Spine","Female_Crude_Total", "Female_Adj_Total","Female_Crude_Neck","Female_Adj_Neck","Female_Crude_Spine","Female_Adj_Spine")

MendRand=function(Fstat_lm, OLS_lm, TSLS_tsls, Reg2_lm){
	MRoutput = matrix(nrow=1, ncol=10)
	#FStat
	MRoutput[1] <- summary(Fstat_lm)[[10]][1] #input Fstat F
	MRoutput[2] <- summary(Fstat_lm)[8] #input Fstat R2

	#OLS
	MRoutput[3:5] <- coef(summary(OLS_lm))[2,c(1,2,4)] #input OLS_B, OLS_se, OLS_P

	#TSLS
	MRoutput[6] <- coef(TSLS_tsls)[2] #input TSLS_B
	MRoutput[7] <- sqrt(diag(TSLS_tsls$V))[2] #input TSLS_se
	MRoutput[8] <- 2 * (1 - pt(abs(TSLS_tsls$coefficients[2]/sqrt(diag(TSLS_tsls$V))[2]), TSLS_tsls$n - TSLS_tsls$p)) #input TSLS_P

	#Daurbin-Hausman
	dh_test <- wald.test(b=coef(Reg2_lm),Sigma=vcov(Reg2_lm),Terms=2,df=Reg2_lm$df.residual)
	MRoutput[9] <- dh_test[[6]][[2]][1] #input DH_F
	MRoutput[10] <- dh_test[[6]][[2]][4] #input DH_P
	return ((MRoutput))
	
	#delete it all so you can start again
	rm(Fstat_lm, OLS_lm, TSLS_tsls, Reg2_lm)
}

print(summary(TSLS_tsls))

combfhsgen3andoffforR2$BMDBMI <-as.character(combfhsgen3andoffforR2$BMDBMI)
combfhsgen3andoffforR2$BMDBMI <-as.numeric(combfhsgen3andoffforR2$BMDBMI)



# All crude total
Reg1 <- lm(URATE~SU.SUM,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
combfhsgen3andoffforR2$UA_Res[!is.na(combfhsgen3andoffforR2$SU.SUM) & !is.na(combfhsgen3andoffforR2$URATE) & combfhsgen3andoffforR2$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[1,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMTOT~URATE,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMTOT~URATE,~SU.SUM,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMTOT~UA_Res+URATE,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
)


# All Adjusted total
Reg1 <- lm(URATE~SU.SUM,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
combfhsgen3andoffforR2$UA_Res[!is.na(combfhsgen3andoffforR2$SU.SUM) & !is.na(combfhsgen3andoffforR2$URATE) & combfhsgen3andoffforR2$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[2,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMTOT~URATE+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMTOT~URATE+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,~SU.SUM+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMTOT~UA_Res+URATE+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
)

# All crude neck
Reg1 <- lm(URATE~SU.SUM,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
combfhsgen3andoffforR2$UA_Res[!is.na(combfhsgen3andoffforR2$SU.SUM) & !is.na(combfhsgen3andoffforR2$URATE) & combfhsgen3andoffforR2$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[3,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMNECK~URATE,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMNECK~URATE,~SU.SUM,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMNECK~UA_Res+URATE,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
)


# All Adjusted neck
Reg1 <- lm(URATE~SU.SUM,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
combfhsgen3andoffforR2$UA_Res[!is.na(combfhsgen3andoffforR2$SU.SUM) & !is.na(combfhsgen3andoffforR2$URATE) & combfhsgen3andoffforR2$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[4,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMNECK~URATE+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMNECK~URATE+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,~SU.SUM+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMNECK~UA_Res+URATE+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
)

# All crude spine
Reg1 <- lm(URATE~SU.SUM,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
combfhsgen3andoffforR2$UA_Res[!is.na(combfhsgen3andoffforR2$SU.SUM) & !is.na(combfhsgen3andoffforR2$URATE) & combfhsgen3andoffforR2$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[5,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(SPINE~URATE,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(SPINE~URATE,~SU.SUM,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	Reg2_lm = lm(SPINE~UA_Res+URATE,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
)


# All Adjusted spine
Reg1 <- lm(URATE~SU.SUM,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
combfhsgen3andoffforR2$UA_Res[!is.na(combfhsgen3andoffforR2$SU.SUM) & !is.na(combfhsgen3andoffforR2$URATE) & combfhsgen3andoffforR2$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[6,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(SPINE~URATE+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(SPINE~URATE+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,~SU.SUM+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3"),
	Reg2_lm = lm(SPINE~UA_Res+URATE+MENOBMD+BMDAGE+BMDBMI+PCA1+PCA2,data=combfhsgen3andoffforR2, subset=COHORT=="GEN3")
)

males <-subset(combfhsgen3andoffforR2, (SEX==1))

# males crude total
Reg1 <- lm(URATE~SU.SUM,data=males, subset=COHORT=="GEN3")
males$UA_Res[!is.na(males$SU.SUM) & !is.na(males$URATE) & males$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[7,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMTOT~URATE,data=males, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMTOT~URATE,~SU.SUM,data=males, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMTOT~UA_Res+URATE,data=males, subset=COHORT=="GEN3")
)


# males Adjusted total
Reg1 <- lm(URATE~SU.SUM,data=males, subset=COHORT=="GEN3")
males$UA_Res[!is.na(males$SU.SUM) & !is.na(males$URATE) & males$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[8,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMTOT~URATE+BMDAGE+BMDBMI+PCA1+PCA2,data=males, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMTOT~URATE+BMDAGE+BMDBMI+PCA1+PCA2,~SU.SUM+BMDAGE+BMDBMI+PCA1+PCA2,data=males, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMTOT~UA_Res+URATE+BMDAGE+BMDBMI+PCA1+PCA2,data=males, subset=COHORT=="GEN3")
)

# males crude neck
Reg1 <- lm(URATE~SU.SUM,data=males, subset=COHORT=="GEN3")
males$UA_Res[!is.na(males$SU.SUM) & !is.na(males$URATE) & males$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[9,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMNECK~URATE,data=males, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMNECK~URATE,~SU.SUM,data=males, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMNECK~UA_Res+URATE,data=males, subset=COHORT=="GEN3")
)


# males Adjusted neck
Reg1 <- lm(URATE~SU.SUM,data=males, subset=COHORT=="GEN3")
males$UA_Res[!is.na(males$SU.SUM) & !is.na(males$URATE) & males$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[10,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMNECK~URATE+BMDAGE+BMDBMI+PCA1+PCA2,data=males, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMNECK~URATE+BMDAGE+BMDBMI+PCA1+PCA2,~SU.SUM+BMDAGE+BMDBMI+PCA1+PCA2,data=males, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMNECK~UA_Res+URATE+BMDAGE+BMDBMI+PCA1+PCA2,data=males, subset=COHORT=="GEN3")
)

# males crude spine
Reg1 <- lm(URATE~SU.SUM,data=males, subset=COHORT=="GEN3")
males$UA_Res[!is.na(males$SU.SUM) & !is.na(males$URATE) & males$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[11,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(SPINE~URATE,data=males, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(SPINE~URATE,~SU.SUM,data=males, subset=COHORT=="GEN3"),
	Reg2_lm = lm(SPINE~UA_Res+URATE,data=males, subset=COHORT=="GEN3")
)


# males Adjusted spine
Reg1 <- lm(URATE~SU.SUM,data=males, subset=COHORT=="GEN3")
males$UA_Res[!is.na(males$SU.SUM) & !is.na(males$URATE) & males$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[12,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(SPINE~URATE+BMDAGE+BMDBMI+PCA1+PCA2,data=males, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(SPINE~URATE+BMDAGE+BMDBMI+PCA1+PCA2,~SU.SUM+BMDAGE+BMDBMI+PCA1+PCA2,data=males, subset=COHORT=="GEN3"),
	Reg2_lm = lm(SPINE~UA_Res+URATE+BMDAGE+BMDBMI+PCA1+PCA2,data=males, subset=COHORT=="GEN3")
)

female <-subset(combfhsgen3andoffforR2, (SEX==2))
female$MENOBMD <-as.character(female$MENOBMD)
female$MENOBMD <-as.factor(female$MENOBMD)

# female crude total
Reg1 <- lm(URATE~SU.SUM,data=female, subset=COHORT=="GEN3")
female$UA_Res[!is.na(female$SU.SUM) & !is.na(female$URATE) & female$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[13,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMTOT~URATE,data=female, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMTOT~URATE,~SU.SUM,data=female, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMTOT~UA_Res+URATE,data=female, subset=COHORT=="GEN3")
)


# female Adjusted total
Reg1 <- lm(URATE~SU.SUM,data=female, subset=COHORT=="GEN3")
female$UA_Res[!is.na(female$SU.SUM) & !is.na(female$URATE) & female$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[14,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMTOT~URATE+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,data=female, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMTOT~URATE+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,~SU.SUM+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,data=female, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMTOT~UA_Res+URATE+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,data=female, subset=COHORT=="GEN3")
)

# female crude neck
Reg1 <- lm(URATE~SU.SUM,data=female, subset=COHORT=="GEN3")
female$UA_Res[!is.na(female$SU.SUM) & !is.na(female$URATE) & female$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[15,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMNECK~URATE,data=female, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMNECK~URATE,~SU.SUM,data=female, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMNECK~UA_Res+URATE,data=female, subset=COHORT=="GEN3")
)


# female Adjusted neck
Reg1 <- lm(URATE~SU.SUM,data=female, subset=COHORT=="GEN3")
female$UA_Res[!is.na(female$SU.SUM) & !is.na(female$URATE) & female$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[16,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(FEMNECK~URATE+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,data=female, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(FEMNECK~URATE+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,~SU.SUM+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,data=female, subset=COHORT=="GEN3"),
	Reg2_lm = lm(FEMNECK~UA_Res+URATE+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,data=female, subset=COHORT=="GEN3")
)

# female crude spine
Reg1 <- lm(URATE~SU.SUM,data=female, subset=COHORT=="GEN3")
female$UA_Res[!is.na(female$SU.SUM) & !is.na(female$URATE) & female$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[17,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(SPINE~URATE,data=female, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(SPINE~URATE,~SU.SUM,data=female, subset=COHORT=="GEN3"),
	Reg2_lm = lm(SPINE~UA_Res+URATE,data=female, subset=COHORT=="GEN3")
)


# female Adjusted spine
Reg1 <- lm(URATE~SU.SUM,data=female, subset=COHORT=="GEN3")
female$UA_Res[!is.na(female$SU.SUM) & !is.na(female$URATE) & female$COHORT=="GEN3"] <- residuals(Reg1)
output_OS[18,2:11] = MendRand(
	Fstat_lm = Reg1,
	OLS_lm = lm(SPINE~URATE+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,data=female, subset=COHORT=="GEN3"),
	TSLS_tsls= tsls(SPINE~URATE+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,~SU.SUM+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,data=female, subset=COHORT=="GEN3"),
	Reg2_lm = lm(SPINE~UA_Res+URATE+BMDAGE+MENOBMD+BMDBMI+PCA1+PCA2,data=female, subset=COHORT=="GEN3")
)







####
## after you've done everything
write.table(output_OS,file="Offspring UNweighted urate GRS.txt",quote=F,row.names=F,sep="\t")
