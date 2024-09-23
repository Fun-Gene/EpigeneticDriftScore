# EDS_NEG
options(stringsAsFactors= F)

# ref data
ref_stat <- readRDS("ref_neg_cpgbeta_meansd.rds")
ref_drift <- readRDS("ref_neg_drift_meansd.rds")
ref_coef <- readRDS("ref_neg_drift_coef.rds")

EDS_NEG <- function(input, ref_stat, ref_drift, ref_coef){
		# overlap and missing
		input_cpg <- colnames(input)
		score_cpg <- ref_coef$cpg
		opcpg <- as.character(score_cpg[score_cpg %in% input_cpg])
		nacpg <- as.character(score_cpg[!(score_cpg %in% input_cpg)])
		oppub <- as.matrix(input[, opcpg])
		misscpg <- colnames(oppub)[unique(which(is.na(oppub),arr.ind=T)[,2])]

		ref_index <- which(ref_stat[,1] %in% opcpg)
		usedmean <- as.numeric(t(ref_stat[,2]))[ref_index]
		names(usedmean) <- ref_stat[ref_index, 1]
		usedsd <- as.numeric(t(ref_stat[,3]))[ref_index]
		names(usedsd) <- ref_stat[ref_index, 1]

		# cpg_score
		cpg_score <- ((as.matrix(oppub[,names(usedmean)]) - matrix(usedmean,nrow=nrow(oppub),ncol=length(usedmean),byrow=T)))^2

		# missing inputation 
		if (sum(is.na(cpg_score)) > 0){
			# missing sample
			cpg_score2 <- cpg_score
			for(i in colnames(cpg_score2)){
			cpg_score2[,i][which(is.na(cpg_score2[,i]))] <-  ref_drift[,2][which(ref_drift$cpg==i)]
			}
			cpg_score <- cpg_score2
		}

		if (length(nacpg) > 0){
			# missing CpG
			mean_score <- ref_drift$res2mean
			index <- which(ref_drift$cpg %in% nacpg)
			fillscore <- data.frame(t(mean_score))[index]
			names(fillscore) <- nacpg
			combine_score <- data.frame(cbind(cpg_score, fillscore))

			} else{ # combine
			combine_score <- cpg_score
		}

		# weight
		pub_score <- (as.matrix(combine_score[,ref_coef$cpg]) %*% ref_coef$x)[,1]

		# scale 
		range01 <- function(x){(x-(-32.72674))/(170.824-(-32.72674))}

		eds <- range01(pub_score)

return(eds)
}