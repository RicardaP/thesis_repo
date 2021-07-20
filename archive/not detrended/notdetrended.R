
data <- data_daily_imputed
N <- length(unique(data$participantID))

# empty data frame to store path estimates in long format (two rows per participant, 1 for T1, 1 for T2)
networks_notdetrended <- data.frame("ID" = rep(1:N, each=2),
                       "timeperiod" = rep(c("T1", "T2"), N),
                       # 3:17 Contemporaneous network: lower triangle in col by col format
                       "PCC_RankA_RankB" = NA,
                       "PCC_RankA_RankC" = NA,
                       "PCC_RankA_Task" = NA,
                       "PCC_RankA_RankE" = NA,
                       "PCC_RankA_RankF" = NA,
                       "PCC_RankB_RankC" = NA,
                       "PCC_RankB_RankD" = NA,
                       "PCC_RankB_RankE" = NA,
                       "PCC_RankB_RankF" = NA,
                       "PCC_RankC_RankD" = NA,
                       "PCC_RankC_RankE" = NA,
                       "PCC_RankC_Ranks2" = NA,
                       "PCC_RankD_RankE" = NA,
                       "PCC_RankD_RankF" = NA,
                       "PCC_RankE_RankF" = NA,
                       # 18:53 Temporal networks: full matrix in col by col format
                       "PDC_RankA_RankA" = NA,
                       "PDC_RankA_RankB" = NA,
                       "PDC_RankA_RankC" = NA,
                       "PDC_RankA_RankD" = NA,
                       "PDC_RankA_RankE" = NA,
                       "PDC_RankA_RankF" = NA,
                       "PDC_RankB_RankA" = NA,
                       "PDC_RankB_RankB" = NA,
                       "PDC_RankB_RankC" = NA,
                       "PDC_RankB_RankD" = NA,
                       "PDC_RankB_RankE" = NA,
                       "PDC_RankB_RankF" = NA,
                       "PDC_RankC_RankA" = NA,
                       "PDC_RankC_RankB" = NA,
                       "PDC_RankC_RankC" = NA,
                       "PDC_RankC_RankD" = NA,
                       "PDC_RankC_RankE" = NA,
                       "PDC_RankC_RankF" = NA,
                       "PDC_RankD_RankA" = NA,
                       "PDC_RankD_RankB" = NA,
                       "PDC_RankD_RankC" = NA,
                       "PDC_RankD_RankD" = NA,
                       "PDC_RankD_RankE" = NA,
                       "PDC_RankD_RankF" = NA,
                       "PDC_RankE_RankA" = NA,
                       "PDC_RankE_RankB" = NA,
                       "PDC_RankE_RankC" = NA,
                       "PDC_RankE_RankD" = NA,
                       "PDC_RankE_RankE" = NA,
                       "PDC_RankE_RankF"   = NA,
                       "PDC_RankF_RankA" = NA,
                       "PDC_RankF_RankB" = NA,
                       "PDC_RankF_RankC" = NA,
                       "PDC_RankF_RankD" = NA,
                       "PDC_RankF_RankE" = NA,
                       "PDC_RankF_RankF" = NA,
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
  if( p %in% c(10, 40, 43)) { next } # partx did not converge
  # filter data by participant
  d_T1 <- filter(data, participantID==p, day.response %in% 1:50)
  d_T2 <- filter(data, participantID==p, day.response %in% 51:100)
  #  select two highest ranking DPDS/behavior variables, higher rank stat = better
  selectVar  <- filter(EMA_descriptives, 
                       EMA_descriptives$ID == p, 
                       EMA_descriptives$emaLabels %in% 
                         emaLabels[-c(61:63)]
  )
  selectVar <- selectVar[order(selectVar$rank_stat_detrended, decreasing = TRUE), ]
  RankA       <- selectVar[1, 3]
  RankB       <- selectVar[2, 3]
  RankC       <- selectVar[3, 3]
  RankD       <- selectVar[4, 3]
  RankE       <- selectVar[5, 3]
  RankF       <- selectVar[6, 3]
  # log name of DPDS Variables
  comparisons_notdetrended$RankA_notdetrended[p] <- RankA
  comparisons_notdetrended$RankB_notdetrended[p] <- RankB
  comparisons_notdetrended$RankC_notdetrended[p] <- RankC
  comparisons_notdetrended$RankD_notdetrended[p] <- RankD
  comparisons_notdetrended$RankE_notdetrended[p] <- RankE
  comparisons_notdetrended$RankF_notdetrended[p] <- RankF
  # log total nr imputed data points
  comparisons_notdetrended$nr_imputed[p] <- filter(data_daily_raw, participantID==p, day.response %in% 1:100) %>%
    select( all_of(RankA), 
            all_of(RankB),
            all_of(RankC), 
            all_of(RankD),
            all_of(RankE), 
            all_of(RankF)) %>%
    is.na() %>%
    sum()
  # define individualized networks
  d_T1 <- select(d_T1, 
                 day.response, 
                 all_of(RankA), 
                 all_of(RankB),
                 all_of(RankC), 
                 all_of(RankD),
                 all_of(RankE), 
                 all_of(RankF)
  )
  d_T2 <- select(d_T2, 
                 day.response, 
                 all_of(RankA), 
                 all_of(RankB),
                 all_of(RankC), 
                 all_of(RankD),
                 all_of(RankE), 
                 all_of(RankF)
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
  networks_notdetrended$BIC_notdetrended[networks_notdetrended$ID == p] <- rbind( netw_T1_notdetrended$EBIC, netw_T2_notdetrended$EBIC)
}

networks_notdetrended$nr_zero_PCC_notdetrended        <- rowSums(networks_notdetrended[ ,idx_PCC] == 0)
networks_notdetrended$nr_zero_PDC_notdetrended        <- rowSums(networks_notdetrended[ ,idx_PDC] == 0)
networks_notdetrended$nr_estimated_PCC_notdetrended   <- rowSums(networks_notdetrended[ ,idx_PCC] != 0)
networks_notdetrended$nr_estimated_PDC_notdetrended   <- rowSums(networks_notdetrended[ ,idx_PDC] != 0)
networks_notdetrended$pr_zero_PCC_notdetrended        <- rowMeans(networks_notdetrended[ ,idx_PCC] == 0) * 100 %>% as.numeric() %>%round( digits = 2)
networks_notdetrended$pr_zero_PDC_notdetrended        <- rowMeans(networks_notdetrended[ ,idx_PDC] == 0) * 100 %>% as.numeric() %>%round( digits = 2)
networks_notdetrended$nr_pos_PCC_notdetrended         <- rowSums(networks_notdetrended[ ,idx_PCC] > 0)
networks_notdetrended$nr_pos_PDC_notdetrended         <- rowSums(networks_notdetrended[ ,idx_PDC] > 0)
networks_notdetrended$nr_neg_PCC_notdetrended         <- rowSums(networks_notdetrended[ ,idx_PCC] < 0)
networks_notdetrended$nr_neg_PDC_notdetrended         <- rowSums(networks_notdetrended[ ,idx_PDC] < 0)
networks_notdetrended$connectivity_PCC_notdetrended   <- round( networks_notdetrended$nr_estimated_PCC /length(idx_PCC) * 100, 2)
networks_notdetrended$connectivity_PDC_notdetrended   <- round( networks_notdetrended$nr_estimated_PDC /length(idx_PDC) * 100, 2)

sum(networks_notdetrended$pr_zero_PCC_notdetrended == 100) # 1 empty contemp net total
sum(networks_notdetrended$pr_zero_PDC_notdetrended == 100) # 46 empty temp net total

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
  comparisons_notdetrended$PCCs_cor_notdetrended[p] <- cor(PCC_T1, PCC_T2, method = "spearman")  
  comparisons_notdetrended$PDCs_cor_notdetrended[p] <- cor(PDC_T1, PDC_T2, method = "spearman") 
  # Jaccard similarity, https://en.wikipedia.org/wiki/Jaccard_index#Similarity_of_asymmetric_binary_attribute
  # proportion of shared non-zero edges relative to number of edges which are non-zero in either, but not both, networks 
  comparisons_notdetrended$PCCs_Jaccard_notdetrended[p] <-   sum(PCC_T1!=0 & PCC_T2!=0) / sum(PCC_T1!=0 | PCC_T2!=0)
  comparisons_notdetrended$PDCs_Jaccard_notdetrended[p] <-   sum(PDC_T1!=0 & PDC_T2!=0) / sum(PDC_T1!=0 | PDC_T2!=0)
  # proportion empty edges both netw
  comparisons_notdetrended$PCCs_prop_empty_notdetrended[p]     <- mean( PCC_T1 == 0 & PCC_T2 == 0)
  comparisons_notdetrended$PDCs_prop_empty_notdetrended[p]     <- mean( PDC_T1 == 0 & PDC_T2 == 0)
  # proportion of edges with equal sign or both zeroin network
  comparisons_notdetrended$PCCs_prop_equ_sign_notdetrended[p]  <- mean( (PCC_T1 == 0 & PCC_T2 == 0) |
                                                (PCC_T1 > 0 & PCC_T2 > 0) |
                                                (PCC_T1 < 0 & PCC_T2 < 0))
  comparisons_notdetrended$PDCs_prop_equ_sign_notdetrended[p]  <- mean( (PDC_T1 == 0 & PDC_T2 == 0) |
                                                (PDC_T1 > 0 & PDC_T2 > 0) |
                                                (PDC_T1 < 0 & PDC_T2 < 0))
  comparisons_notdetrended$avg_BIC_notdetrended[p] <- networks_notdetrended[networks_notdetrended$ID==p, 54] %>% mean()
}

saveRDS(comparisons_notdetrended, here("data", "comparisons_notdetrended.RDS"))


######################################################################
#### rename participant ID with same pipeline used for daily data ####
data2 <- read.csv(here("data", "dd100_baseline.csv")) %>% select(-X)
# included IDs
idx_included <- idx_missing %>% as.data.frame() %>% filter(Include == TRUE)
# exclude participants
data2 <- filter(data2, participantID %in% idx_included$ID)
# rename IDs for easy looping
data2$participantID <- factor(data2$participantID, 
                             levels = unique(data2$participantID), 
                             labels = 1:length(unique(data2$participantID))) %>% 
  as.numeric()
# N: nr participants retained
N1 <- length(unique(data2$participantID)) #79
data_baseline <- data2
######################################################################

data_merged_notdetrended <- inner_join(data_baseline, comparisons_notdetrended, by = "participantID")  

fit_notdetrended <- lm(formula = PCCs_cor_notdetrended ~ 1 + gender + age + happy + swlsMean + neoN + neoE + neoO + neoA + neoC + nr_imputed, data = data_merged_notdetrended
)
summary(fit_notdetrended)
plot(data_merged_notdetrended$PCCs_cor_notdetrended, data_merged_notdetrended$neoN)
hist(data_merged_notdetrended$PCCs_cor_notdetrended)


cor(data_merged_notdetrended$PCCs_cor_notdetrended, data_merged$PCCs_cor, use = "pairwise")
