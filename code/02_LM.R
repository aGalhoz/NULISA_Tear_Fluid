source("00_data.R")
source("ML_functions.R")

# data basal & reflex
protein_data_basalvsreflex =  protein_data %>%
  filter(Tears %in% c("basal","reflex") & SampleMatrixType == "OTHER" & !is.na(Diagnosis)) %>%
  dplyr::select(Target,UniProtID,Patient,Tears,Diagnosis,normNPQ_median)
protein_data_basalvsreflex_old <- protein_data_basalvsreflex

protein_data_basalvsreflex_new <- protein_data_basalvsreflex %>% 
  dplyr::select(Target,Patient,Tears,normNPQ_median) %>% 
  distinct() %>%
  filter(!is.na(normNPQ_median)) 

## alter the data (basal vs rflx is called status)

protein_data_basalvsreflex_new <- protein_data_basalvsreflex_new %>%
  pivot_wider(names_from = Target,values_from = normNPQ_median)
names_patients <- protein_data_basalvsreflex_new$Patient 
protein_data_basalvsreflex_new <- protein_data_basalvsreflex_new %>%
  dplyr::select(-Patient) %>% 
  dplyr::rename(status = Tears)
protein_data_basalvsreflex_new <- protein_data_basalvsreflex_new[,c(2:ncol(protein_data_basalvsreflex_new),1)]
df_ml = protein_data_basalvsreflex_new %>%
  mutate(status = ifelse(status == "basal",1,0))

# Run Lasso
lm = runML(df_ml,'lm',BS_number = 500)
saveRDS(lm, file = 'data_output/linearModel.rds')

# Run Lasso with random selection of balanced groups
lm_balanced = runML_balanced(df_ml,'lm',BS_number = 500,
                             condition_1 = 1,
                             condition_2 = 0)
saveRDS(lm_balanced, file = 'data_output/linearModel_balanced.rds')

# plot averaged ROC curve
ROC_curve = calculateROC(lm,df_ml,'plots/rocc_lm.pdf')
write.csv(ROC_curve$roc_data, file = 'data_output/rocc_lm.csv')

ROC_curve_balanced = calculateROC_balanced(lm_balanced,df_ml,'plots/rocc_lm_balanced.pdf')
write.csv(ROC_curve_balanced$roc_data, file = 'data_output/rocc_lm_balanced.csv')

# extract weights + plot averaged
lm_weights = plotWeights(lm,"plots/weights_lm.pdf",number = 37)
write.csv(lm_weights, file = 'data_output/weights_lm.csv')

lm_weights_balanced = plotWeights(lm_balanced,"plots/weights_lm_balanced.pdf",number = 25)
write.csv(lm_weights_balanced, file = 'data_output/weights_lm_balanced.csv')

# Lasso output and DE from Basal vs reflex
DE_basal_reflex_LM <- DE_basal_reflex %>%
  left_join(lm_weights %>% mutate(protein = row.names(.)))
writexl::write_xlsx(DE_basal_reflex_LM, 'data_output/DE_basal_reflex_LM.xlsx')

DE_basal_reflex_LM_balanced <- DE_basal_reflex %>%
  left_join(lm_weights_balanced %>% mutate(protein = row.names(.)))
writexl::write_xlsx(DE_basal_reflex_LM_balanced, 
                    'data_output/DE_basal_reflex_LM_balanced.xlsx')

# plot of weights average vs frequency
lm_weights <- lm_weights %>% mutate(label = row.names(.))
pdf("plots/LM_meanweight_frequency.pdf", width = 8, height = 6)
ggplot(lm_weights, aes(x = avg, y = freq, label = label)) +
  geom_point(aes(size = abs(avg), color = avg > 0), alpha = 0.7) +
  scale_color_manual(values = c("red", "blue"), labels = c("Negative", "Positive")) +
  geom_text(data = subset(lm_weights, freq > 0.1 | abs(avg) > 7),
            vjust = -0.5, hjust = 0.5, size = 3) +
  labs(title = "Coefficients' average weight vs frequency",
       y = "Selection Frequency",
       x = "Average weight",
       color = "Direction",
       size = "|Average weight|") +
  theme_minimal()
dev.off()

lm_weights_balanced <- lm_weights_balanced %>% mutate(label = row.names(.))
pdf("plots/LM_meanweight_frequency_balanced.pdf", width = 8, height = 6)
ggplot(lm_weights_balanced, aes(x = avg, y = freq, label = label)) +
  geom_point(aes(size = abs(avg), color = avg > 0), alpha = 0.7) +
  scale_color_manual(values = c("red", "blue"), labels = c("Negative", "Positive")) +
  geom_text(data = subset(lm_weights_balanced, freq > 0.06 | abs(avg) > 7),
            vjust = -0.5, hjust = 0.5, size = 3) +
  labs(title = "Coefficients' average weight vs frequency",
       y = "Selection Frequency",
       x = "Average weight",
       color = "Direction",
       size = "|Average weight|") +
  theme_minimal()
dev.off()



  