library(readr)
library(Maaslin2)

Alldiet1 <- read_csv("DietData_Imputed1.csv")
finalMicrobiome_1_ <- read_csv("finalMicrobiome (1).csv")
input_data <- read.table("data_genus.tsv", header = TRUE, sep = "\t",row.names = 1)
input_metadata <- read.table("metadata_1.tsv", header = TRUE, sep = "\t",row.names = 1)
dim(Alldiet)
dim(finalMicrobiome_1_)
dim(input_metadata)
finalMicrobiome_1_<-finalMicrobiome_1_[,6:ncol(finalMicrobiome_1_)]
finalMicrobiome_df <- as.data.frame(finalMicrobiome_1_)
rownames(finalMicrobiome_df) <- rownames(input_metadata)

input_metadata$Body.Mass.Index <- as.numeric(input_metadata$Body.Mass.Index)
input_metadata$Age..in.years. <- as.numeric(input_metadata$Age..in.years.)
input_metadata$AlcoholServings<- as.numeric(input_metadata$AlcoholServings)
input_metadata$fiber<- as.numeric(input_metadata$fiber)

#include diabet_condition
fit_data_Microbiome <- Maaslin2(input_data = finalMicrobiome_df, 
                                input_metadata = input_metadata, 
                                output = "maaslin2 results_finalllllll", 
                                fixed_effects = c("CaseControl","Age..in.years.","Body.Mass.Index"),
                                random_effects = c("BatchEffect","Gender","Race"),
                                reference = c("CaseControl,Control"),
                                max_significance=0.5,
                                min_prevalence = 0.3,
                                cores = 10
)
taxa_maslin_results11111111<-fit_data_Microbiome$results
taxa_maslin_results<-fit_data_Microbiome$results
sorted_taxa_maslin <- taxa_maslin_results[order(taxa_maslin_results$pval), ]
sorted_taxa_maslin <- sorted_taxa_maslin[sorted_taxa_maslin$pval < 0.05, ]



boxplot <- ggplot(sorted_taxa_maslin, aes(x = sorted_taxa_maslin$feature, y = sorted_taxa_maslin$coef, fill = sorted_taxa_maslin$metadata)) +
  geom_boxplot() +
  labs(x = "Group", y = "Abundance", title = "Abundance of Microbial Taxa Across Groups (Boxplot)") +
  theme_minimal()

taxa_maaslin<-unique(sorted_taxa_maslin$feature)
taxa_maaslin<-finalMicrobiome_df[c("Haemophilus", "Gordonibacter","Veillonella", "Allisonella", "Dielma", 
                                  "Rhodococcus", "Bacteroides", "Marvinbryantia", 
                                  "Dialister", "Negativibacillus", "Lachnospiraceae_FCS020_group", 
                                  "Turicibacter", "Parabacteroides", "Butyricimonas", "Romboutsia", "Stenotrophomonas", "Coprococcus_1", "Granulicatella", 
                                  "Prevotella_9", "Streptococcus", "Ruminiclostridium_6", 
                                  "Lachnospiraceae_NK4A136_group","Ruminococcaceae_UCG-014","Lachnospiraceae_UCG-010","Ruminococcaceae_UCG-009")]

summary(taxa_maaslin)
Alldiet<-Alldiet1[,9:ncol(Alldiet1)]
Alldiet <- as.data.frame(Alldiet)
rownames(Alldiet) <- rownames(input_metadata)
fit_data_Alldiet <- Maaslin2(input_data = Alldiet, 
                                input_metadata = input_metadata, 
                                output = "diet_maaslin2", 
                                fixed_effects = c("CaseControl","Age..in.years.","Body.Mass.Index"),
                                random_effects = c("BatchEffect","Gender","Race"),
                                reference = c("CaseControl,Control"),
                                max_significance=0.9,
                                min_prevalence = 0.1,
                                cores = 10
)

#for output
fit_data_Alldiet <- Maaslin2(input_data = sorted_diet, 
                             input_metadata = input_metadata, 
                             output = "missingtotalfat_maaslin2", 
                             fixed_effects = c("CaseControl","Age..in.years.","Body.Mass.Index"),
                             random_effects = c("BatchEffect","Gender","Race"),
                             reference = c("CaseControl,Control"),
                             max_significance=0.9,
                             min_prevalence = 0.1,
                             cores = 10
)
#for output
fit_data_Alldiet <- Maaslin2(input_data = top_20_taxa, 
                             input_metadata = input_metadata, 
                             output = "top20taxa_plots", 
                             fixed_effects = c("CaseControl","Age..in.years.","Body.Mass.Index"),
                             random_effects = c("BatchEffect","Gender","Race"),
                             reference = c("CaseControl,Control"),
                             max_significance=0.9,
                             min_prevalence = 0.1,
                             cores = 10
)
fit_data_Alldiet <- Maaslin2(input_data = taxa_maaslin, 
                             input_metadata = input_metadata, 
                             output = "taxa_maaslin2plots", 
                             fixed_effects = c("CaseControl","Age..in.years.","Body.Mass.Index"),
                             random_effects = c("BatchEffect","Gender","Race"),
                             reference = c("CaseControl,Control"),
                             max_significance=0.9,
                             min_prevalence = 0.1,
                             cores = 10
)
diet_maslin<-fit_data_Alldiet$results
sorted_diet_maslin <- diet_maslin[order(diet_maslin$pval), ]
sorted_diet_maslin <- sorted_diet_maslin[sorted_diet_maslin$pval < 0.05, ]

sorted_diet_maslin<-unique(sorted_diet_maslin$feature)


sorted_diet<-Alldiet[c("NonFriedFishServings", "M_FISH_HI", "pfa184", 
  "FishServings", "pfa205", "pfa226", "vitd", "M_FRANK", "sfa100","sfa40", "M_FISH_LO", "sfa80", 
  "protanim", "vitd3", "sfa120","totaltfa")]
maaslin_diet<-Alldiet[,sorted_diet_maslin]







maaslin_diet_names<-colnames(sorted_diet)
maaslin_diet_names<-as.data.frame(maaslin_diet_names)
maaslin_taxa_names<-colnames(sorted_taxa)
maaslin_taxa_names<-as.data.frame(maaslin_taxa_names)
maaslin_top_20taxa_names<-colnames(top_20_taxa)
maaslin_top_20taxa_names<-as.data.frame(maaslin_top_20taxa_names)


top20_diet1<-Alldiet1[c("sfa40", "M_FISH_LO", "FishServings", "M_FISH_HI", "pfa205", "tfa181t", "sfa170", "sfa100", "NonFriedFishServings", "pfa225", "totaltfa", "sfa220", "alcohol", "AlcoholServings", "omega3", "calories", "pfa204")]
top20_diet<-Alldiet1[c("M_FISH_LO", "sfa40", "sfa100", "M_FISH_HI", "FishServings", "pfa225", "pfa205", "NonFriedFishServings", "vitd3", "pfa204", "omega3", "sfa80", "sfa120", "mfa161", "calories", "tfa181t", "sfa170", "pfa184")]
rf_diet<-cbind(top20_diet,top20_diet1)
names<-names(rf_diet)
rf_diet<-unique(names)
rf_diet<-Alldiet1[,rf_diet]


union_diet<-cbind(sorted_diet,top20_diet,top20_diet1)
names<-names(union_diet)
union_diet<-unique(names)
union_diet<-Alldiet1[,union_diet]

mediators<-cbind(sorted_taxa,top_20_taxa)
names<-colnames(mediators)
names<-as.data.frame(names)
mediators<-unique(names)
mediators<-finalMicrobiome_df[,mediators]


############################ single diet- single microbe
library(openxlsx)
library(mediation)
c<-Alldiet1$CaseControl
covariates <- Alldiet1[, 3:6]
top_20_taxa<-finalMicrobiome_df[c("Veillonella", "Gordonibacter", "Haemophilus", "Parabacteroides", "Agathobacter", "Bacteroides", "Ruminococcaceae_UCG-014", "Anaerotruncus", "Faecalibacterium", "Streptococcus", "Blautia", "Megasphaera", "Allisonella", "Other.1", "Dielma", "Marvinbryantia", "Alistipes", "Lachnoclostridium", "Gemella")]

cdiet <- cbind(c,covariates, exposures, mediators)
cdiet<-as.data.frame(cdiet)


results_list <- list()

for (j in names(cdiet[, 6:31])) {
  evar <- cdiet[[j]]
  expvar <- list()
  
  for (i in names(cdiet[, 32:65])) {
    dvar <- cdiet[[i]]
    
    # Perform mediation analysis
    model.M <- lm(dvar ~ evar, data = cdiet)
    model.Y <- lm(c ~ dvar + evar + Age + Race + Gender + BMI, data = cdiet)
    result <- mediate(model.M, model.Y, treat = "evar", mediator = "dvar", boot = TRUE, sims = 100)
    result_summary <- summary(result)
    
    # Extract relevant information
    acme_coef <- result_summary$d0
    acme_ci <- result_summary$d0.ci
    acme_p_value <- result_summary$d0.p
    ade_coef <- result_summary$z0
    ade_ci <- result_summary$z0.ci
    ade_p_value <- result_summary$z0.p
    te_coef <- result_summary$n0
    te_ci <- result_summary$n0.ci
    te_p_value <- result_summary$n0.p
    prop_coef <- result_summary$tau.coef
    prop_ci <- result_summary$tau.ci
    prop_p_value <- result_summary$tau.p
    
    # Create a dataframe with the extracted information
    expvar[[i]] <- data.frame(
      Taxa = i,
      diet= j
      ACME_Coefficient = acme_coef,
      ACME_Confidence_Interval_Lower = acme_ci[1],
      ACME_Confidence_Interval_Upper = acme_ci[2],
      ACME_P_Value = acme_p_value,
      ADE_Coefficient = ade_coef,
      ADE_Confidence_Interval_Lower = ade_ci[1],
      ADE_Confidence_Interval_Upper = ade_ci[2],
      ADE_P_Value = ade_p_value,
      TE_Coefficient = te_coef,
      TE_Confidence_Interval_Lower = te_ci[1],
      TE_Confidence_Interval_Upper = te_ci[2],
      TE_P_Value = te_p_value,
      Prop_Mediate_Coefficient = prop_coef,
      Prop_Mediate_Confidence_Interval_Lower = prop_ci[1],
      Prop_Mediate_Confidence_Interval_Upper = prop_ci[2],
      Prop_Mediate_P_Value = prop_p_value
      )
    
    combined_results <- do.call(rbind, expvar)
    results_list[[j]] <- combined_results
    combined_results <- do.call(rbind,results_list)
    }
  }

    
}
combined_results$names<-rownames(combined_results)
sorted_combined_results <- combined_results[combined_results$ACME_P_Value < 0.1, ]

#write.xlsx(combined_results, "union_singlediet-singlemed.xlsx")


########################### single diet with alpha diversity
library(readr)
alphaDiversity <- read_csv("alphaDiversity.csv")
shanon<-alphaDiversity$shanon
cdiet <- cbind(cdiet,shanon)
cdiet<-as.data.frame(cdiet)

results_list <- list()

# Loop through each exposure variable
for (exposure_variable in names(cdiet[,6:31])) {
  
  # Extract exposure variable
  exposure <- cdiet[[exposure_variable]]
  
  model.M <- lm(shanon ~ exposure, data = cdiet)
  model.Y <- lm(c ~ shanon + exposure + Age + Race + Gender + BMI, data = cdiet)
  result <- mediate(model.M, model.Y, treat = "exposure", mediator = "shanon", boot = TRUE, sims = 100)
  result_summary <- summary(result)
  
  # Extract relevant information
  acme_coef <- result_summary$d0
  acme_ci <- result_summary$d0.ci
  acme_p_value <- result_summary$d0.p
  ade_coef <- result_summary$z0
  ade_ci <- result_summary$z0.ci
  ade_p_value <- result_summary$z0.p
  te_coef <- result_summary$n0
  te_ci <- result_summary$n0.ci
  te_p_value <- result_summary$n0.p
  prop_coef <- result_summary$tau.coef
  prop_ci <- result_summary$tau.ci
  prop_p_value <- result_summary$tau.p
  
  results_list[[exposure_variable]] <- data.frame(
    Exposure_Variable = exposure_variable,
    ACME_Coefficient = acme_coef,
    ACME_Confidence_Interval_Lower = acme_ci[1],
    ACME_Confidence_Interval_Upper = acme_ci[2],
    ACME_P_Value = acme_p_value,
    ADE_Coefficient = ade_coef,
    ADE_Confidence_Interval_Lower = ade_ci[1],
    ADE_Confidence_Interval_Upper = ade_ci[2],
    ADE_P_Value = ade_p_value,
    TE_Coefficient = te_coef,
    TE_Confidence_Interval_Lower = te_ci[1],
    TE_Confidence_Interval_Upper = te_ci[2],
    TE_P_Value = te_p_value,
    Prop_Mediate_Coefficient = prop_coef,
    Prop_Mediate_Confidence_Interval_Lower = prop_ci[1],
    Prop_Mediate_Confidence_Interval_Upper = prop_ci[2],
    Prop_Mediate_P_Value = prop_p_value
  )
}

# Combine all results
combined_results <- do.call(rbind, results_list)

# Write results to Excel
write.xlsx(combined_results, "alpha_diversity_singlediet_singlemed.xlsx")

######################## HILAMA: SINGLE DIET_MULTI TAXA
library(HILAMA)
library(openxlsx)

Y1 <- cdiet[1]
X <- union_diet
M1 <- 

# Create an empty list to store all results
results_list <- list()

# Iterate over each diet variable
for (variable in names(X)) {
  expo<-X[[variable]]
  current_data <- cbind(covariates,expo)
  X1 <- as.matrix(current_data)
  Y1 <- as.matrix(Y1)
  M1 <- as.matrix(M1)
  
  result_ls <- hilama(X = current_data, M_multi = M1, Y = Y1, is.parallel = TRUE, core_num = 5, verbose = TRUE, is.adap = FALSE)
  sorted<-result_ls$pvalue_nde
  nie_mat_hat <- as.data.frame(result_ls$nie_mat_hat)
  pvalue_Theta_mat <- result_ls$pvalue_Theta
  pvalue_beta <- result_ls$pvalue_beta
  
  # Screening with JST test
  screen_N <- ceiling(0.1 * ncol(X1) * ncol(M1))
  result_screen_jst_ls <- get_nie_pvalue_screen_JST(pvalue_Theta_mat, pvalue_beta, screen_N = screen_N, fdr_level = 0.1)
  pvalue_nie_screen_JST_mat <- as.data.frame(result_screen_jst_ls$pvalue_screen_JST_mat)
  
  # Extract ddcoef values
  ddcoef_values <- matrix(nrow = nrow(result_ls$coef_ls$X2M_ls[[1]]), ncol = length(result_ls$coef_ls$X2M_ls))
  colnames(ddcoef_values) <- colnames(M1)
  for (i in seq_along(result_ls$coef_ls$X2M_ls)) {
    matrix_for_microbe <- result_ls$coef_ls$X2M_ls[[i]]
    ddcoef_values[, i] <- matrix_for_microbe$dd_coef
  }
  
  
  indices <- which(pvalue_nie_screen_JST_mat < 0.1, arr.ind = TRUE)
  #result <- ifelse(indices > 3, "Greater than 3", "Less than or equal to 3")
  if (nrow(indices) > 0) {
    Cofounder<-ifelse(indices[, 1] < 5, rownames(nie_mat_hat)[indices[, 1]], "-")
    #Cofounder<-rownames(nie_mat_hat)[indices[, 1]]
    Diet<-variable
    Microbe <- colnames(nie_mat_hat)[indices[, 2]]
    pvalue <- pvalue_nie_screen_JST_mat[indices[, 1], indices[, 2]]
    pvalue_matrix <- as.matrix(pvalue)
    p_values <- diag(pvalue_matrix)
    effect <- nie_mat_hat[indices[, 1], indices[, 2]]
    effect <- as.matrix(effect)
    NIE <- diag(effect)
    
    theta <- ddcoef_values[indices[, 1], indices[, 2]]
    theta <- as.matrix(theta)
    Theta <- diag(theta)
    
    hilama_result_with_trial_taxa <- data.frame(Diet,Cofounder, Microbe, Theta, NIE, p_values)
    results_list[[variable]] <- hilama_result_with_trial_taxa
  }
  else {
    cat("No significant results found for", variable, "\n")
    next  # Move to next iteration
  }
}
plot(result_ls)

# Mediation analysis plot
mediation_plot(result_ls)
# Combine all results into a single dataframe
combined_results <- do.call(rbind, results_list)

################################### HILAMA MULTI DIET-MULTI TAXA
library(HILAMA)
library(openxlsx)

Y1 <- cdiet[1]
X1 <- maaslin_diet
M1 <- taxa_maaslin
X1 <- as.matrix(X1)
Y1 <- as.matrix(Y1)
M1 <- as.matrix(M1)

result_ls <- hilama(X = X1, M_multi = M1, Y = Y1, is.parallel = TRUE, core_num = 5, verbose = TRUE, is.adap = FALSE)

nie_mat_hat <- as.data.frame(result_ls$nie_mat_hat)
pvalue_Theta_mat <- result_ls$pvalue_Theta
pvalue_beta <- result_ls$pvalue_beta

# Screening with JST test
screen_N <- ceiling(0.1 * ncol(X1) * ncol(M1))
result_screen_jst_ls <- get_nie_pvalue_screen_JST(pvalue_Theta_mat, pvalue_beta, screen_N = screen_N, fdr_level = 0.1)
pvalue_nie_screen_JST_mat <- as.data.frame(result_screen_jst_ls$pvalue_screen_JST_mat)

# Extract ddcoef values
ddcoef_values <- matrix(nrow = nrow(result_ls$coef_ls$X2M_ls[[1]]), ncol = length(result_ls$coef_ls$X2M_ls))
colnames(ddcoef_values) <- colnames(M1)
for (i in seq_along(result_ls$coef_ls$X2M_ls)) {
  matrix_for_microbe <- result_ls$coef_ls$X2M_ls[[i]]
  ddcoef_values[, i] <- matrix_for_microbe$dd_coef
}

indices <- which(pvalue_nie_screen_JST_mat < 0.1, arr.ind = TRUE)

if (nrow(indices) > 0) {
  Cofounder<-ifelse(indices[, 1] < 5, rownames(nie_mat_hat)[indices[, 1]], "-")
  Diet<-variable
  Microbe <- colnames(nie_mat_hat)[indices[, 2]]
  pvalue <- pvalue_nie_screen_JST_mat[indices[, 1], indices[, 2]]
  pvalue_matrix <- as.matrix(pvalue)
  p_values <- diag(pvalue_matrix)
  effect <- nie_mat_hat[indices[, 1], indices[, 2]]
  effect <- as.matrix(effect)
  NIE <- diag(effect)
  theta <- ddcoef_values[indices[, 1], indices[, 2]]
  theta <- as.matrix(theta)
  Theta <- diag(theta)
  
  hilama_result_with_trial_taxa <- data.frame(Diet,Cofounder, Microbe, Theta, NIE, p_values)
  results_list[[variable]] <- hilama_result_with_trial_taxa
}
else {
  cat("No significant results found for", variable, "\n")
  next  # Move to next iteration
}


library(HILAMA)
library(openxlsx)

Y1 <- cdiet[1]
X1 <- maaslin_diet
M1 <- taxa_maaslin
X1 <- as.matrix(X1)
Y1 <- as.matrix(Y1)
M1 <- as.matrix(M1)

result_ls <- hilama(X = X1, M_multi = M1, Y = Y1, is.parallel = TRUE, core_num = 5, verbose = TRUE, is.adap = FALSE)

nie_mat_hat <- as.data.frame(result_ls$nie_mat_hat)
pvalue_Theta_mat <- result_ls$pvalue_Theta
pvalue_beta <- result_ls$pvalue_beta

# Screening with JST test
screen_N <- ceiling(0.1 * ncol(X1) * ncol(M1))
result_screen_jst_ls <- get_nie_pvalue_screen_JST(pvalue_Theta_mat, pvalue_beta, screen_N = screen_N, fdr_level = 0.1)
pvalue_nie_screen_JST_mat <- as.data.frame(result_screen_jst_ls$pvalue_screen_JST_mat)

# Extract ddcoef values
ddcoef_values <- matrix(nrow = nrow(result_ls$coef_ls$X2M_ls[[1]]), ncol = length(result_ls$coef_ls$X2M_ls))
colnames(ddcoef_values) <- colnames(M1)
for (i in seq_along(result_ls$coef_ls$X2M_ls)) {
  matrix_for_microbe <- result_ls$coef_ls$X2M_ls[[i]]
  ddcoef_values[, i] <- matrix_for_microbe$dd_coef
}

# Iterate over each diet variable
for (variable in names(X1)) {
  indices <- which(pvalue_nie_screen_JST_mat < 0.1, arr.ind = TRUE)
  
  if (nrow(indices) > 0) {
    Cofounder <- ifelse(indices[, 1] < 5, rownames(nie_mat_hat)[indices[, 1]], "-")
    Diet <- variable
    Microbe <- colnames(nie_mat_hat)[indices[, 2]]
    pvalue <- pvalue_nie_screen_JST_mat[indices[, 1], indices[, 2]]
    pvalue_matrix <- as.matrix(pvalue)
    p_values <- diag(pvalue_matrix)
    effect <- nie_mat_hat[indices[, 1], indices[, 2]]
    effect <- as.matrix(effect)
    NIE <- diag(effect)
    theta <- ddcoef_values[indices[, 1], indices[, 2]]
    theta <- as.matrix(theta)
    Theta <- diag(theta)
    
    hilama_result_with_trial_taxa <- data.frame(Diet, Cofounder, Microbe, Theta, NIE, p_values)
    print(data.frame(diet_v, Microbe,effect,pvalue))
  } else {
    cat("No significant results found for", variable, "\n")
  }
}

###################################################################################
# Load necessary libraries
library(ggplot2)

# Example data
mediators <- c("Mediator 1", "Mediator 2", "Mediator 3")
exposure <- rep(c("Exposure A", "Exposure B", "Exposure C"), each = 3)
effect <- rep(c("ADE", "ACME", "TE"), each = 3)
coefficients <- c(0.3, 0.5, 0.8, 0.4, 0.8, 1.2, 0.8, 1.2, 1.7)  # Coefficients for ADE, ACME, and TE
se <- c(0.1, 0.2, 0.15)  # Standard errors
p_values_ade <- c(0.03, 0.05, 0.01)  # p-values for ADE
p_values_acme <- c(0.02, 0.08, 0.01)  # p-values for ACME
p_values_te <- c(0.05, 0.01, 0.03)  # p-values for TE

# Create data frame
data <- data.frame(mediators, exposure, effect, coefficients, se, p_values_ade, p_values_acme, p_values_te)

# Create forest plot
ggplot(data, aes(y = mediators, x = coefficients, color = exposure, shape = effect)) +
  geom_point() +
  geom_errorbarh(aes(xmin = coefficients - 1.96 * se, xmax = coefficients + 1.96 * se), height = 0) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Coefficient", y = NULL, title = "Forest Plot of Mediation Analysis") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(shape = guide_legend(title = "Effect")) +
  geom_point(data = subset(data, p_values_ade < 0.05), aes(color = NULL, shape = NULL), size = 3) +
  geom_text(data = subset(data, p_values_ade < 0.05), aes(label = paste("p =", round(p_values_ade, 3))), hjust = -0.1)
