#############
setwd(work_dir)

################################
# 1. load clock CpG list
################################

hannum <- read.csv("Hannum71cpg", header = T)
levine <- read.csv("Levine513cpg", header = T)
horvath <- read.csv("Horvath353cpg", header = T)
dunedin <- read.csv("Dunedin46cpg", header = T)

################################
# 2. Dunedin pace of aging 
################################

index <- which(dunedin$id %in% names(cpg))
cpg_beta <- dunedin$coef[index]
names(cpg_beta) <- dunedin$id[index]
cpg_score <- cpg[,names(cpg_beta)]
age_score <- apply(cpg_score,1,function(x){ sum(x * cpg_beta)}) + (-0.06)
cor(age_score,age)
summary(lm(age~age_score))$r.squared

################################
# 3.  First generation clock Hannum 
################################

index <- which(hannum$Marker %in% names(cpg))
cpg_beta <- hannum$Coefficient[index]
names(cpg_beta) <- hannum$Marker[index]
cpg_score <- cpg[,names(cpg_beta)]
age_score <- apply(cpg_score,1,function(x){ sum(x * cpg_beta)})
cor(age_score,age)
summary(lm(age~age_score))$r.squared

################################
# 4. First generation clock Horvath 
################################

index <- which(horvath$CpGmarker %in% names(cpg))
cpg_beta <- horvath$CoefficientTraining[index]
names(cpg_beta) <- horvath$CpGmarker[index]
cpg_score <- cpg[,names(cpg_beta)]
age_score <- apply(cpg_score,1,function(x){ sum(x * cpg_beta)}) + 0.695507258
cor(age_score,age)
summary(lm(age~age_score))$r.squared

################################
# 5. Second generation clock Levine 
################################

index <- which(levine$CpG %in% names(cpg))
cpg_beta <- levine$Weight[index]
names(cpg_beta) <- levine$CpG[index]
cpg_score <- cpg[,names(cpg_beta)]
age_score <- apply(cpg_score,1,function(x){ sum(x * cpg_beta)}) + 60.66400
cor(age_score,age)
summary(lm(age~age_score))$r.squared
