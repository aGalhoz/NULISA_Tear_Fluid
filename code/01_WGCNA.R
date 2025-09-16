source("00_data.R")

# data ALS & CTR
protein_data_ALSvsCTR =  protein_data %>%
  filter(Diagnosis %in% c("ALS","CTRL") & Tears == "basal" & SampleMatrixType == "OTHER") %>%
  select(Target,UniProtID,Patient,Diagnosis,normNPQ_median)
protein_data_ALSvsCTR_old <- protein_data_ALSvsCTR

protein_data_ALSvsCTR_new <- protein_data_ALSvsCTR %>% 
  select(Target,Patient,normNPQ_median) %>% 
  distinct() %>%
  filter(!is.na(normNPQ_median)) 

protein_data_ALSvsCTR_new <- protein_data_ALSvsCTR_new %>%
  pivot_wider(names_from = Patient,values_from = normNPQ_median)

# expression data for WGCNA
expression_data = protein_data_ALSvsCTR_new[,2:ncol(protein_data_ALSvsCTR_new)]
rownames(expression_data) = protein_data_ALSvsCTR_new[,1] %>% pull()
expression_data = t(expression_data)

## WGCNA
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Network topology analysis function
sft = pickSoftThreshold(
  expression_data,  
  powerVector = powers,
  verbose = 5,
  dataIsExpr = T
)

# Plot scale independence and mean connectivity together
pdf("plots/WGCNA_thresholds.pdf",width = 10, height = 8)
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(1, 1.5))  # Make second row taller
cex1 = 0.9;

# 1st plot: scale independence
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = expression("Scale-free Topology Fit Index (R"^2*")"),
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red")

#2nd plot: mean connectivity
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# 3rd plot: all together
plot(sft$fitIndices[, 1], sft$fitIndices[, 2], type = "o", col = "blue", pch = 16,
     xlab = "Soft-threshold (power)",
     ylab = expression("Scale-free Topology Fit Index (R"^2*")"),
     ylim = c(0, 1))
abline(h = 0.8, col = "gray", lty = 2)
par(new = TRUE)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "o", col = "red", pch = 17,
     axes = FALSE, xlab = "", ylab = "", log = "")
axis(side = 4, col = "red", col.axis = "red")
mtext("Mean Connectivity", side = 4, line = 3, col = "red")
title("WGCNA Soft-threshold Selection")
grid()
legend("topright", legend = c("Scale of independence", "Mean Connectivity"),
       col = c("blue", "red"), pch = c(16, 17), lty = 1)
dev.off()

## comment: ideally pick 8-10 powers, tested all the powers, always returns 3 modules

# Construct WGCNA network
picked_power = 8
temp_cor <- WGCNA::cor        
network <- blockwiseModules(expression_data,   
                          # == Adjacency Function ==
                          power = picked_power,               
                          networkType = "signed",
                          TOMType = "signed",
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 20,
                          maxBlockSize = 4000,
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor

# Dendrogram with module colors
color_palette <- c("#1B9E77","#D95F02","#7570B3")
mergedColors <- color_palette[network$colors + 1]
pdf(file = "plots/dendrogram_WGCNA.pdf",width = 8, height = 6)
plotDendroAndColors(
  network$dendrograms[[1]],
  mergedColors[network$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
dev.off()

# Relate modules with ALS/CTR groups
module_df <- data.frame(
  protein_id = names(network$colors),
  colors = c("green","orange","purple")[network$colors+1]
)

writexl::write_xlsx(module_df,"data_output/protein_modules.xlsx")

# Get module eingenproteins per cluster
MEs0 <- moduleEigengenes(expression_data,mergedColors)$eigengenes
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
MEs0$patient = row.names(MEs0)

# correlation of patients with a specific module
mME = MEs0 %>%
  pivot_longer(-patient) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)) %>%
  left_join(protein_data_ALSvsCTR %>% 
              select(Patient,Diagnosis) %>%
              rename(patient = Patient) %>%
              distinct()) %>%
  left_join(data.frame(name = c("#1B9E77","#D95F02","#7570B3"),
                       color = c("green","orange","purple")))

mME_summary <- mME %>%
  group_by(color, Diagnosis) %>%
  summarise(
  mean_value = mean(value),
   sd_value = sd(value),
  count = n(),
  .groups = 'drop'
   )
writexl::write_xlsx(mME_summary,"data_output/correlation_summary_module.xlsx")

pdf("plots/Module-patient relationship_2.pdf")
mME %>% arrange(Diagnosis) %>% ggplot(., aes(x=patient, y=color, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis_c(limit = c(-0.45,0.45)) +
  # scale_fill_gradient2(
  #   low = "blue",
  #   high = "red",
  #   mid = "white",
  #   midpoint = 0,
  #   limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")
dev.off()

pdf("plots/Module-group relationship_2.pdf")
mME_summary %>% ggplot(., aes(x=Diagnosis, y=color, fill=mean_value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_viridis_c(limit = c(-0.45,0.45)) +
  # scale_fill_gradient2(
  #   low = "blue",
  #   high = "red",
  #   mid = "white",
  #   midpoint = 0,
  #   limit = c(-0.1,0.1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="Mean corr")
dev.off()

submod_df = data.frame(t(expression_data)) %>%
  mutate(
    protein_id = row.names(.)
  ) %>%
  pivot_longer(-protein_id) %>%
  left_join(module_df %>% rename(module = colors)) %>%
  left_join(protein_data_ALSvsCTR %>% 
              select(Patient,Diagnosis) %>%
              rename(name = Patient) %>%
              distinct())

pdf("plots/expression by treatment groups.pdf")
submod_df %>% ggplot(., aes(x=Diagnosis, y=value, group=protein_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  scale_color_manual(values = c("purple" = "#7570B3",  # blue
    "green" = "#1B9E77",  # orange
    "orange" = "#D95F02"))  + # green 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment group",
       y = "normalized expression")
dev.off()

pdf("plots/expression by patients.pdf")
submod_df %>% ggplot(., aes(x=name, y=value, group=protein_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  scale_color_manual(values = c("purple" = "#7570B3",  # blue
                                "green" = "#1B9E77",  # orange
                                "orange" = "#D95F02"))  + # green 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment group",
       y = "normalized expression")
dev.off()

# Get network to plot in Cytoscape
TOM = TOMsimilarityFromExpr(expression_data,power = picked_power)
row.names(TOM) = colnames(TOM) = colnames(expression_data)

# edge list
edge_list = data.frame(TOM) %>%
  mutate(
    protein1 = row.names(.)
  ) %>%
  pivot_longer(-protein1) %>%
  dplyr::rename(protein2 = name, correlation = value) %>%
  unique() %>%
  subset(!(protein1==protein2)) %>%
  left_join(module_df %>% rename(protein1 = protein_id,
                                 module1 = colors)) %>%
  left_join(module_df %>% rename(protein2 = protein_id,
                                 module2 = colors))
 
writexl::write_xlsx(edge_list,"data_output/edge_list.xlsx")

# check module structure and expression level based on DE
module_DE_df <- module_df %>%
  rename(protein = protein_id) %>%
  left_join(DE_ALS_CTR) %>%
  arrange(`p-value`)

writexl::write_xlsx(module_DE_df,"data_output/DE_modules.xlsx")

# enrichment of genes per module with clusterProfiler
gse_module_green <- enrichGO(gene= module_df %>% filter(colors == "green") %>% pull(protein_id), #51 proteins
                    ont ="ALL", universe = module_df %>% pull(protein_id),
                    OrgDb = organism,
                    keyType = "SYMBOL",
                    # keyType = 'UNIPROT',
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.1)
# significance in BP
writexl::write_xlsx(gse_module_green@result,"data_output/GO_module_green.xlsx")
gse_module_purple <- enrichGO(gene= module_df %>% filter(colors == "purple") %>% pull(protein_id), #25 proteins
                              ont ="ALL", universe = module_df %>% pull(protein_id),
                              OrgDb = organism,
                              keyType = "SYMBOL",
                              # keyType = 'UNIPROT',
                              minGSSize = 3, 
                              maxGSSize = 800, 
                              pvalueCutoff = 0.1)
# result: no significance
writexl::write_xlsx(gse_module_purple@result,"data_output/GO_module_purple.xlsx")
gse_module_orange <- enrichGO(gene= module_df %>% filter(colors == "orange") %>% pull(protein_id), #49 proteins
                              ont ="ALL", universe = module_df %>% pull(protein_id),
                              OrgDb = organism,
                              keyType = "SYMBOL",
                              # keyType = 'UNIPROT',
                              minGSSize = 3, 
                              maxGSSize = 800, 
                              pvalueCutoff = 0.1)
# significance in CC and MF
writexl::write_xlsx(gse_module_orange@result,"data_output/GO_module_orange.xlsx")

# GO: Biological Processes
analysis <- gse_module_green@result
GP_BP_green = analysis %>% filter(ONTOLOGY == "BP")
GP_BP_green = GP_BP_green %>% mutate(log10_FDR = -log10(p.adjust)) %>%
  arrange(p.adjust)
GP_BP_green$Description <- factor(GP_BP_green$Description,levels = GP_BP_green$Description)
pdf("plots/GO_BP_green_module.pdf") 
ggplot(GP_BP_green,aes(Description,log10_FDR,fill = log10_FDR)) + 
  geom_bar(stat = "identity")  + 
  scale_fill_gradient(low = "blue", high = "red") + coord_flip() + theme_minimal() +
  ylab("-log10(FDR)") + xlab("Biological Process Term")
dev.off()

# GO: Cellular Component
analysis <- gse_module_orange@result
GP_CC_orange = analysis %>% filter(ONTOLOGY == "CC")
GP_CC_orange = GP_CC_orange %>% mutate(log10_FDR = -log10(p.adjust)) %>%
  arrange(p.adjust)
GP_CC_orange$Description <- factor(GP_CC_orange$Description,levels = GP_CC_orange$Description)
pdf("plots/GP_CC_orange_module.pdf") 
ggplot(GP_CC_orange,aes(Description,log10_FDR,fill = log10_FDR)) + 
  geom_bar(stat = "identity")  + 
  scale_fill_gradient(low="blue", high="red") + coord_flip() + theme_minimal() +
  ylab("-log10(FDR)") + xlab("Cellular Component Term")
dev.off()

# GO: Molecular Function
GP_MF_orange = analysis %>% filter(ONTOLOGY == "MF")
GP_MF_orange = GP_MF_orange %>% mutate(log10_FDR = -log10(p.adjust)) %>%
  arrange(p.adjust)
GP_MF_orange$Description <- factor(GP_MF_orange$Description,levels = GP_MF_orange$Description)
pdf("plots/GP_MF_orange_module.pdf") 
ggplot(GP_MF_orange,aes(Description,log10_FDR,fill = log10_FDR)) + 
  geom_bar(stat = "identity")  + 
  scale_fill_gradient(low="blue", high="red") + coord_flip() + theme_minimal() +
  ylab("-log10(FDR)") + xlab("Molecular Function Term")
dev.off()

# ------------------------------------------------
## Correlation of modules with clinical variables 
# -> NfL or pNfH
clinical_data$NfL # NfL with less missing but there is no intersection between them
clinical_data$pNfH
# -> Progression rate group (high or low)
median_PR <- median(clinical_data$PR,na.rm = T)
clinical_data$PR_group = ifelse(clinical_data$PR <= median_PR,"low","high")
# Age at onset group (high or low)
median_AO = median(clinical_data$AAO,na.rm = T)
clinical_data$AAO_group = ifelse(clinical_data$AAO <= median_AO,"low","high")

mME_clinical = MEs0 %>%
  pivot_longer(-patient) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)) %>%
  left_join(clinical_data %>% 
              dplyr::select(Patient,Diagnosis,Site,ALSFRSR_1,
                     VC_percent,NfL,PR_group,AAO_group,RL_total,
                     Eye_disease_cat,Eye_medication,CL) %>%
              dplyr::rename(patient = Patient) %>%
              distinct()) %>%
  left_join(data.frame(name = c("#1B9E77","#D95F02","#7570B3"),
                       color = c("green","orange","purple")))

module_matrix <- mME_clinical %>%
  dplyr::select(patient, color, value) %>%
  pivot_wider(names_from = color, values_from = value) %>%
  column_to_rownames(var = "patient")

clinical <- mME_clinical %>%
  dplyr::select(patient, Diagnosis, Site, ALSFRSR_1, VC_percent, NfL,
                PR_group,AAO_group,RL_total,
                Eye_disease_cat, Eye_medication, CL) %>%
  distinct(patient, .keep_all = TRUE) %>%
  column_to_rownames(var = "patient")

# Ensure the rows match
clinical <- clinical[rownames(module_matrix), ]

# -----------------------------
# correlation and significance matrices
cor_matrix <- matrix(NA, nrow=ncol(module_matrix), ncol=ncol(clinical))
pval_matrix <- matrix(NA, nrow=ncol(module_matrix), ncol=ncol(clinical))
rownames(cor_matrix) <- colnames(module_matrix)
rownames(pval_matrix) <- colnames(module_matrix)
colnames(cor_matrix) <- colnames(clinical)
colnames(pval_matrix) <- colnames(clinical)

for(mod in colnames(module_matrix)){
  ME <- module_matrix[[mod]]
  for(clin_var in colnames(clinical)){
    cv <- clinical[[clin_var]]
    
    if(is.numeric(cv)){
      # Spearman correlation
      tmp_cor <- cor(ME, cv, use="pairwise.complete.obs", method="spearman")
      test <- cor.test(ME, cv, method="spearman", exact=FALSE)
      cor_matrix[mod, clin_var] <- tmp_cor
      pval_matrix[mod, clin_var] <- test$p.value
    } else {
      # Multi-level categorical → Kruskal-Wallis test
      test <- kruskal.test(ME ~ as.factor(cv))
      pval_matrix[mod, clin_var] <- test$p.value
      cor_matrix[mod, clin_var] <- NA
    }
  }
}

cor_matrix <- as.data.frame(cor_matrix)
pval_matrix <- as.data.frame(pval_matrix)

# heatmaps of correlation and significance
cor_matrix_numeric <- cor_matrix[, colSums(is.na(cor_matrix)) == 0]
pdf("plots/heatmap_WGCNA_correlation_clinical.pdf")
pheatmap(cor_matrix_numeric, cluster_rows=TRUE, cluster_cols=TRUE, display_numbers=TRUE,
         main="Module Eigengene Correlations",
         color = colorRampPalette(c("blue", "white", "darkred"))(100),
         breaks = seq(-1, 1, length.out = 101))
dev.off()

pdf("plots/heatmap_WCGNA_pvalue_clinical.pdf")
pheatmap(-log10(pval_matrix), cluster_rows=TRUE, cluster_cols=TRUE, display_numbers=TRUE,
         main="-log10(p-values) for Module vs Clinical",
         color = colorRampPalette(c("white", "red"))(100))
dev.off()

# -----------------------------
# boxplots/linear slope of numerical and categorical variables per module
modules <- rownames(cor_matrix_numeric)
clin_vars <- colnames(clinical)  
MEs = MEs0 
colnames(MEs)= c("orange","green","purple","patient")

# Loop over all module–clinical variable combinations
for(module in modules){
  
  # Identify numeric vs categorical variables
  numeric_vars <- clin_vars[sapply(clinical[, clin_vars], is.numeric)]
  categorical_vars <- setdiff(colnames(clinical), numeric_vars)
  
  if(length(numeric_vars) > 0){
  # Data and plot for numerical variables
  df_num = clinical %>%
    rownames_to_column("patient") %>%
    dplyr::select(patient, all_of(numeric_vars)) %>%
    left_join(MEs %>% dplyr::select(patient, all_of(module)) %>%
                rename_with(~"ModuleEigengene", all_of(module)),
              by = "patient") %>%
    pivot_longer(cols = -c(patient, ModuleEigengene),
                 names_to = "ClinicalVar", values_to = "Value_num") %>%
    filter(!is.na(Value_num)) 
  
  annotations_num <- df_num %>%
    group_by(ClinicalVar) %>%
    summarise(
      rho = round(cor(ModuleEigengene, Value_num, method="spearman"), 2),
      pval = signif(cor.test(ModuleEigengene, Value_num, method="spearman")$p.value, 4),
      slope = format(coef(lm(ModuleEigengene ~ Value_num))[2], scientific = TRUE, digits = 2),
      label = paste0("Slope = ", slope, "\nSpearman r = ", rho, "\nP = ", pval),
      x_pos = median(Value_num, na.rm=TRUE),
      y_pos = max(ModuleEigengene, na.rm=TRUE) * 1.2,
      .groups = "drop"
    )
  
  p_num <- ggplot(df_num, aes(x=Value_num, y=ModuleEigengene)) +
    geom_point(alpha=0.7, color="black") +
    geom_smooth(method="lm", se=TRUE, color="blue") +
    geom_text(data=annotations_num, aes(x=x_pos, y=y_pos, label=label),
              inherit.aes=FALSE, hjust=0.5, vjust=0, size=3) +
    facet_wrap(~ ClinicalVar, scales="free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          plot.margin = margin(1,1,2,1,"cm")) +
    labs(title=paste(module, "module across numeric clinical variables"),
         y="Module Eigengene", x="")
  
  ggsave(filename=paste0("plots/", module, "_numeric_variables.pdf"),
         plot=p_num, width=12, height=8)
  }
  
  if(length(categorical_vars) > 0){
  # data for categorical variables 
    df_cat <- clinical %>%
      rownames_to_column("patient") %>%
      dplyr::select(patient, all_of(categorical_vars)) %>%
      left_join(MEs %>% dplyr::select(patient, all_of(module)) %>% 
                  rename_with(~"ModuleEigengene", all_of(module)),
                by="patient") %>%
      pivot_longer(cols = all_of(categorical_vars), 
                   names_to = "ClinicalVar", values_to = "Value_fac") %>%
      filter(!is.na(Value_fac) & Value_fac != "0") %>%
      mutate(Value_fac = as.factor(Value_fac))
  
    annotations_cat <- df_cat %>%
      group_by(ClinicalVar) %>%
      summarise(
        pval = signif(kruskal.test(ModuleEigengene ~ Value_fac)$p.value, 4),
        label = paste0("P = ", pval),
        x_pos = 1,
        y_pos = max(ModuleEigengene, na.rm=TRUE) * 1.05,
        .groups = "drop"
      )
    
    # Categorical plot
    p_cat <- ggplot(df_cat, aes(x=Value_fac, y=ModuleEigengene, fill=Value_fac)) +
      geom_boxplot(alpha=0.7) +
      geom_jitter(aes(color=Value_fac), width=0.2, alpha=0.7) +
      scale_fill_viridis_d(option = "turbo")  +
      scale_color_viridis_d(option = "turbo") +
      geom_text(data=annotations_cat, aes(x=x_pos, y=y_pos, label=label),
                inherit.aes=FALSE, hjust=0.5, vjust=0, size=3) +
      facet_wrap(~ ClinicalVar, scales="free") +
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(title=paste(module, "module vs categorical clinical variables"),
           y="Module Eigengene", x="")
    
    ggsave(filename=paste0("plots/", module, "_categorical_variables.pdf"),
           plot=p_cat, width=12, height=8)
  }
  
}

# -----------------------------
# clinical variables and modules of interest: RL_total, eye disease (orange),
# RL_total, NfL (purple) 

# orange
proteins_orange <- module_df %>%
  filter(colors == "orange") %>%
  pull(protein_id)

# purple
proteins_purple <- module_df %>%
  filter(colors == "purple") %>%
  pull(protein_id)

# green
proteins_green <- module_df %>%
  filter(colors == "green") %>%
  pull(protein_id)

# module clinical pairs of interest
module_clinical_pairs = list(
  list(module = "purple", proteins = proteins_purple, clinical = "NfL"),
  list(module = "purple", proteins =  proteins_orange, clinical = "RL_total"),
  list(module = "orange", proteins =  proteins_orange, clinical = "RL_total"),
  list(module = "orange", proteins =  proteins_green, clinical = "Eye_disease_cat"))

plot_module_vs_clinical <- function(module_name, clinical_var_name,
                                    expression_data, module_df, clinical_data,
                                    outfolder = "plots") {
  
  if(!dir.exists(outfolder)) dir.create(outfolder)
  
  # Subset expression data for proteins in this module
  module_proteins <- colnames(expression_data)[module_df$colors == module_name]
  if(length(module_proteins) == 0){
    message("No proteins in module ", module_name, ". Skipping.")
    return(NULL)
  }
  
  df <- as.data.frame(expression_data[, module_proteins, drop=FALSE])
  df$Patient <- rownames(expression_data)
  
  # Pivot to long format
  df_long <- df %>%
    tidyr::pivot_longer(cols = -Patient, names_to = "Protein", values_to = "Expression") %>%
    dplyr::left_join(clinical_data %>% dplyr::select(Patient, all_of(clinical_var_name)),
                     by = "Patient") %>%
    dplyr::filter(!is.na(Expression) & !is.na(.data[[clinical_var_name]]) &
                    .data[[clinical_var_name]]!="0")
  
  if(nrow(df_long) == 0){
    message("No valid data for module ", module_name, " vs ", clinical_var_name)
    return(NULL)
  }
  
  is_numeric <- is.numeric(df_long[[clinical_var_name]])
  
  # Compute per-protein stats
  annotations <- df_long %>%
    dplyr::group_by(Protein) %>%
    dplyr::group_modify(~ {
      dat <- .x  
      n <- nrow(dat)
      if(n < 2) return(tibble(label = NA_character_, y_pos = NA_real_))
      
      y_var <- dat[[clinical_var_name]]
      expr_vals <- dat[!is.na(y_var) & y_var != "0",]$Expression
      y_var <- y_var[!is.na(y_var) & y_var != "0"]
      print(y_var)
      
      if(is.numeric(y_var)){
        rho <- round(cor(expr_vals, y_var, method = "spearman", use="pairwise.complete.obs"), 2)
        pval <- signif(cor.test(expr_vals, y_var, method = "spearman")$p.value, 3)
        slope <- round(coef(lm(expr_vals ~ y_var))[2], 4)
        label <- paste0("Slope = ", slope, "\nSpearman r = ", rho, "\nP = ", pval)
      } else {
        # Categorical variable: Kruskal-Wallis
        pval <- signif(kruskal.test(expr_vals ~ as.factor(y_var))$p.value, 3)
        label <- paste0("Kruskal-Wallis P = ", pval)
      }
      
      y_pos <- max(expr_vals, na.rm = TRUE) * 1.05
      tibble(label = label, y_pos = y_pos)
    }) %>%
    dplyr::ungroup()
  
  # Start ggplot
  p <- ggplot2::ggplot(df_long, aes(y = Expression))
  
  if(is_numeric){
    p <- p +
      geom_point(aes(x = .data[[clinical_var_name]]), alpha=0.7) +
      geom_smooth(aes(x = .data[[clinical_var_name]]), method="lm", color="blue", se=TRUE)
    annotations$x_pos <- median(df_long[[clinical_var_name]], na.rm=TRUE)
  } else {
    p <- p +
      geom_boxplot(aes(x = as.factor(.data[[clinical_var_name]]), fill = as.factor(.data[[clinical_var_name]])), alpha=0.7) +
      geom_jitter(aes(x = as.factor(.data[[clinical_var_name]]), color = as.factor(.data[[clinical_var_name]])), width=0.2, alpha=0.7) +
      scale_fill_brewer(palette="Set2") +
      scale_color_brewer(palette="Dark2")
    # x-position for annotations = first level of factor
    annotations$x_pos <- levels(as.factor(df_long[[clinical_var_name]]))[1]
  }
  
  p <- p +
    facet_wrap(~ Protein, scales = "free_y") +
    geom_text(data = annotations,
              aes(x = x_pos, y = y_pos, label = label),
              inherit.aes = FALSE, hjust=0.5, vjust=0, size=3) +
    labs(title = paste(module_name, "module vs", clinical_var_name),
         x = clinical_var_name, y = "Expression") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  
  outname <- paste0(outfolder, "/", module_name, "_", clinical_var_name, "_all_proteins.pdf")
  message("Saved PDF: ", outname)
  
  return(p)
}

# Generate and save PDFs
for(pair in module_clinical_pairs){
  p <- plot_module_vs_clinical(pair$module, 
                               pair$clinical,
                               expression_data,
                               module_df,
                               clinical_data)
  ggsave(filename = paste0("plots/", pair$module, "_", pair$clinical, ".pdf"),
         plot = p, width = 16, height = 14)
}


