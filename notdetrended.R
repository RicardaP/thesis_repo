
data <- data_daily_imputed

# empty data frame to store path estimates in long format (two rows per participant, 1 for T1, 1 for T2)
networks_notdetrended <- data.frame("ID" = rep(1:N, each=2),
                       "timeperiod" = rep(c("T1", "T2"), N),
                       # 3:17 Contemporaneous network: lower triangle in col by col format
                       "PCC_PA_NA" = NA,
                       "PCC_PA_Stress" = NA,
                       "PCC_PA_Task" = NA,
                       "PCC_PA_Rank1" = NA,
                       "PCC_PA_Rank2" = NA,
                       "PCC_NA_Stress" = NA,
                       "PCC_NA_tasks" = NA,
                       "PCC_NA_Rank1" = NA,
                       "PCC_NA_Rank2" = NA,
                       "PCC_Stress_tasks" = NA,
                       "PCC_Stress_Rank1" = NA,
                       "PCC_Stress_Ranks2" = NA,
                       "PCC_Tasks_Rank1" = NA,
                       "PCC_Tasks_Rank2" = NA,
                       "PCC_Rank1_Rank2" = NA,
                       # 18:53 Temporal networks: full matrix in col by col format
                       "PDC_PA_PA" = NA,
                       "PDC_PA_NA" = NA,
                       "PDC_PA_Stress" = NA,
                       "PDC_PA_Tasks" = NA,
                       "PDC_PA_Rank1" = NA,
                       "PDC_PA_Rank2" = NA,
                       "PDC_NA_PA" = NA,
                       "PDC_NA_NA" = NA,
                       "PDC_NA_Stress" = NA,
                       "PDC_NA_Tasks" = NA,
                       "PDC_NA_Rank1" = NA,
                       "PDC_NA_Rank2" = NA,
                       "PDC_Stress_PA" = NA,
                       "PDC_Stress_NA" = NA,
                       "PDC_Stress_Stress" = NA,
                       "PDC_Stress_Tasks" = NA,
                       "PDC_Stress_Rank1" = NA,
                       "PDC_Stress_Rank2" = NA,
                       "PDC_Tasks_PA" = NA,
                       "PDC_Tasks_NA" = NA,
                       "PDC_Tasks_Stress" = NA,
                       "PDC_Tasks_Tasks" = NA,
                       "PDC_Tasks_Rank1" = NA,
                       "PDC_Tasks_Rank2" = NA,
                       "PDC_Rank1_PA" = NA,
                       "PDC_Rank1_NA" = NA,
                       "PDC_Rank1_Stress" = NA,
                       "PDC_Rank1_Tasks" = NA,
                       "PDC_Rank1_Rank1" = NA,
                       "PDC_Rank1_Rank2"   = NA,
                       "PDC_Rank2_PA" = NA,
                       "PDC_Rank2_NA" = NA,
                       "PDC_Rank2_Stress" = NA,
                       "PDC_Rank2_Tasks" = NA,
                       "PDC_Rank2_Rank1" = NA,
                       "PDC_Rank2_Rank2" = NA,
                       # 54
                       "BIC" = NA
)

# index columns for contemporaneous network (PPC: 3:17) and temporal networks (PDC 18:53) in Networks df
idx_PCC <- 3:17
idx_PDC <- 18:53
# data frame for main outcome statistics
comparisons_notdetrended <- data.frame(participantID = 1:N)

# Specification a la BeckJohnson2020, gamma = 0, Lambda range 0.025 to 1
for (p in 1:N){
  # if( p %in% c(4)) { next }
  # filter data by participant
  d_T1 <- filter(data, participantID==p, day.response %in% 1:50)
  d_T2 <- filter(data, participantID==p, day.response %in% 51:100)
  #  select two highest ranking DPDS/behavior variables, higher rank stat = better
  selectVar  <- filter(EMA_descriptives, 
                       EMA_descriptives$ID == p, 
                       EMA_descriptives$emaLabels %in% 
                         emaLabels[c(11:23, 31:60)]
  )
  selectVar <- selectVar[order(selectVar$rank_stat_detrended, decreasing = TRUE), ]
  Rank1       <- selectVar[1, 3]
  Rank2       <- selectVar[2, 3]
  # log name of DPDS Variables
  comparisons_notdetrended$Rank1[p] <- Rank1
  comparisons_notdetrended$Rank2[p] <- Rank2
  # log total nr imputed data points
  comparisons_notdetrended$nr_imputed[p] <- filter(data_daily_raw, participantID==p, day.response %in% 1:100) %>%
    select(PosAffect, NegAffect, StressSev, tasks, all_of(Rank1), all_of(Rank2)) %>%
    is.na() %>%
    sum()
  # define individualized networks
  d_T1 <- select(d_T1, 
                 day.response, 
                 PosAffect, 
                 NegAffect, 
                 StressSev, 
                 tasks, 
                 all_of(Rank1), 
                 all_of(Rank2)
  )
  d_T2 <- select(d_T2, 
                 day.response, 
                 PosAffect, 
                 NegAffect, 
                 StressSev, 
                 tasks, 
                 all_of(Rank1), 
                 all_of(Rank2)
  ) 
  # BeckJohnson2020 used gamma = 0, Lambda 0.025:1 by 0.025
  netw_T1_notdetrended <- graphicalVAR(d_T1, 
                          beepvar = "day.response",
                          gamma = 0,
                          lambda_beta = seq(.025, 1, .0125)
  )
  cat(paste("T1 participant",p ,"converged"))
  netw_T2_notdetrended <- graphicalVAR(d_T2, 
                          beepvar = "day.response",
                          gamma = 0,
                          lambda_beta = seq(.025, 1, .0125)
  )
  cat(paste("T2 participant",p ,"converged"))
  # store full network objects in global environment for later inspection
  assign(paste0("ID",p,"_T1_notdetrended"), netw_T1_notdetrended)
  assign(paste0("ID",p,"_T2_notdetrended"), netw_T2_notdetrended)
  # save network objects for detailed inspection / plotting
  saveRDS(netw_T1_notdetrended, here("data", "network_objects_notdetrended", paste0("ID",p,"_T1_notdetrended.RDS")))
  saveRDS(netw_T2_notdetrended, here("data", "network_objects_notdetrended", paste0("ID",p,"_T2_notdetrended.RDS")))
  # extract path estimates for easy comparison
  networks_notdetrended[networks_notdetrended$ID == p, idx_PCC] <-
    rbind(netw_T1_notdetrended$PCC[lower.tri(netw_T1_notdetrended$PCC)] %>%
            as.numeric()%>%
            as.vector() %>%
            round(digits=4),
          netw_T2_notdetrended$PCC[lower.tri(netw_T2_notdetrended$PCC)] %>%
            as.numeric()%>%
            as.vector() %>%
            round(digits=4)
    )
  networks_notdetrended[networks_notdetrended$ID == p, idx_PDC] <- 
    rbind(netw_T1_notdetrended$PDC %>%
            as.numeric()%>%
            as.vector()%>%
            round(digits = 4),
          netw_T2_notdetrended$PDC %>%
            as.numeric()%>%
            as.vector()%>%
            round(digits = 4)
    )
  networks_notdetrended$BIC[networks_notdetrended$ID == p] <- rbind( netw_T1_notdetrended$EBIC, netw_T2_notdetrended$EBIC)
}

networks_notdetrended$nr_zero_PCC        <- rowSums(networks_notdetrended[ ,idx_PCC] == 0)
networks_notdetrended$nr_zero_PDC        <- rowSums(networks_notdetrended[ ,idx_PDC] == 0)
networks_notdetrended$nr_estimated_PCC   <- rowSums(networks_notdetrended[ ,idx_PCC] != 0)
networks_notdetrended$nr_estimated_PDC   <- rowSums(networks_notdetrended[ ,idx_PDC] != 0)
networks_notdetrended$pr_zero_PCC        <- rowMeans(networks_notdetrended[ ,idx_PCC] == 0) * 100 %>% as.numeric() %>%round( digits = 2)
networks_notdetrended$pr_zero_PDC        <- rowMeans(networks_notdetrended[ ,idx_PDC] == 0) * 100 %>% as.numeric() %>%round( digits = 2)
networks_notdetrended$nr_pos_PCC         <- rowSums(networks_notdetrended[ ,idx_PCC] > 0)
networks_notdetrended$nr_pos_PDC         <- rowSums(networks_notdetrended[ ,idx_PDC] > 0)
networks_notdetrended$nr_neg_PCC         <- rowSums(networks_notdetrended[ ,idx_PCC] < 0)
networks_notdetrended$nr_neg_PDC         <- rowSums(networks_notdetrended[ ,idx_PDC] < 0)
networks_notdetrended$connectivity_PCC   <- round( networks_notdetrended$nr_estimated_PCC /length(idx_PCC) * 100, 2)
networks_notdetrended$connectivity_PDC   <- round( networks_notdetrended$nr_estimated_PDC /length(idx_PDC) * 100, 2)

sum(networks_notdetrended$pr_zero_PCC == 100) # 1 empty contemp net total
sum(networks_notdetrended$pr_zero_PDC == 100) # 46 empty temp net total

saveRDS(networks_notdetrended, file = here("data", "networks_notdetrended.RDS"))

for (p in 1:N) {
  # dynamic index for T1 and T2 rows per person in 'networks_notdetrended' data frame
  idx_p_T1  <- (p-1)*2+1
  idx_p_T2  <- (p-1)*2+2
  PCC_T1    <- networks_notdetrended[idx_p_T1, idx_PCC] %>% as.numeric()
  PCC_T2    <- networks_notdetrended[idx_p_T2, idx_PCC] %>% as.numeric()
  PDC_T1    <- networks_notdetrended[idx_p_T1, idx_PDC] %>% as.numeric()
  PDC_T2    <- networks_notdetrended[idx_p_T2, idx_PDC] %>% as.numeric()
  # correlation of edge weights
  comparisons_notdetrended$PCCs_cor[p] <- cor(PCC_T1, PCC_T2, method = "spearman")  
  comparisons_notdetrended$PDCs_cor[p] <- cor(PDC_T1, PDC_T2, method = "spearman") 
  # Jaccard similarity, https://en.wikipedia.org/wiki/Jaccard_index#Similarity_of_asymmetric_binary_attribute
  # proportion of shared non-zero edges relative to number of edges which are non-zero in either, but not both, networks 
  comparisons_notdetrended$PCCs_Jaccard[p] <-   sum(PCC_T1!=0 & PCC_T2!=0) / sum(PCC_T1!=0 | PCC_T2!=0)
  comparisons_notdetrended$PDCs_Jaccard[p] <-   sum(PDC_T1!=0 & PDC_T2!=0) / sum(PDC_T1!=0 | PDC_T2!=0)
  # proportion empty edges both netw
  comparisons_notdetrended$PCCs_prop_empty[p]     <- mean( PCC_T1 == 0 & PCC_T2 == 0)
  comparisons_notdetrended$PDCs_prop_empty[p]     <- mean( PDC_T1 == 0 & PDC_T2 == 0)
  # proportion of edges with equal sign or both zeroin network
  comparisons_notdetrended$PCCs_prop_equ_sign[p]  <- mean( (PCC_T1 == 0 & PCC_T2 == 0) |
                                                (PCC_T1 > 0 & PCC_T2 > 0) |
                                                (PCC_T1 < 0 & PCC_T2 < 0))
  comparisons_notdetrended$PDCs_prop_equ_sign[p]  <- mean( (PDC_T1 == 0 & PDC_T2 == 0) |
                                                (PDC_T1 > 0 & PDC_T2 > 0) |
                                                (PDC_T1 < 0 & PDC_T2 < 0))
  comparisons_notdetrended$avg_BIC[p] <- networks_notdetrended[networks_notdetrended$ID==p, 54] %>% mean()
}

saveRDS(comparisons_notdetrended, here("data", "comparisons_notdetrended.RDS"))


######################################################################
#### rename participant ID with same pipeline used for daily data ####
data <- read.csv(here("data", "dd100_baseline.csv")) %>% select(-X)
# included IDs
idx_included <- idx_missing %>% as.data.frame() %>% filter(Include == TRUE)
# exclude participants
data <- filter(data, participantID %in% idx_included$ID)
# rename IDs for easy looping
data$participantID <- factor(data$participantID, 
                             levels = unique(data$participantID), 
                             labels = 1:length(unique(data$participantID))) %>% 
  as.numeric()
# N: nr participants retained
N1 <- length(unique(data$participantID)) #79
data_baseline <- data
######################################################################

data_merged_notdetrended <- inner_join(data_baseline, comparisons_notdetrended, by = "participantID")  

fit_notdetrended <- lm(formula = PCCs_cor ~ 1 + gender + age + happy + swlsMean + neoN + neoE + neoO + neoA + neoC + nr_imputed + avg_BIC, data = data_merged
)
summary(fit_notdetrended)
