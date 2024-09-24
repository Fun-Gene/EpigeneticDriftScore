# compute the cumulative score and correlation
get_cumu_rep_cor <- function(cpg, nspt_cpg, age, step, rep ) {
                      # number of CpGs
                      sample_num <- c(seq(0,length(nspt_cpg),step),length(nspt_cpg))
                      # set cor
                      cor <- data.frame(matrix(NA, rep, length(sample_num)))
                      colnames(cor) <- paste0("CpG",sample_num)

                      
                      for(i in 1:length(sample_num)) {
                               # step
                               if(i==1 & sample_num[i]==0) {
                                    cor[,i] <- 0
                                    }
                               else if (i>1 & sample_num[i]>1){ #stepwise
                                    #set correlation
                                    cor_rep <- c(NA)
                                    for(j in 1:rep){
                                                #sample CpG
                                                cpg_index <- sample(1:length(nspt_cpg), sample_num[i])
                                                #calculate score
                                                out <- nnls::nnls(as.matrix(cpg[,cpg_index]),age)
                                                # weight
                                                sum(out$x>0)
                                                # predict
                                                cpg_score_age <- as.matrix(cpg[,cpg_index]) %*% out$x
                                                # correlation
                                                cor_rep[j] <- cor(cpg_score_age, age)
                                                }
                                    cor[,i] <- cor_rep
                                  }
                            }
                       return(cor)
}
