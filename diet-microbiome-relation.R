
```{r}
diet_variables<-df[, 4:ncol(df)]
correlation_results <- list()

for (diet_variable in names(diet_variables)) {
  # Extract the diet variable
  current_diet_variable <- diet_variables[[diet_variable]]
  
  # Check if the diet variable is a vector
  if (is.vector(current_diet_variable)) {
    # Perform Spearman's rank correlation with adjustments
    cor_result <- cor.test(df$shanon, current_diet_variable, method = "spearman")
    
    # Store the correlation result
    correlation_results[[diet_variable]] <- cor_result
  } else {
    warning(paste("Skipping", diet_variable, "as it is not a vector."))
  }
}

# Print the correlation results
for (diet_variable in names(correlation_results)) {
  print(paste("Correlation for", diet_variable, ":"))
  print(correlation_results[[diet_variable]])
}
# Assuming correlation_results is a list containing correlation results
# Replace "correlation_results" with the actual list variable name

# Create a data frame from the correlation results
correlation_df <- data.frame(
  Diet_Variable = names(correlation_results),
  Correlation_Coefficient = sapply(correlation_results, function(x) x$estimate),
  P_Value = sapply(correlation_results, function(x) x$p.value)
)



# Filter only significant variables
correlation_df <- correlation_df[order(correlation_df$P_Value), ]
significant_df <- subset(correlation_df, P_Value < 0.05)
significant_df 
# Plot bar chart for significant correlation coefficients
library(ggplot2)
ggplot(significant_df, aes(x = Diet_Variable, y = Correlation_Coefficient, fill = P_Value < 0.05)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.2f", Correlation_Coefficient)),
            vjust = -0.5, position = position_dodge(width = 0.9), size = 3) +
  labs(title = "Correlation between Diet Variables and Disease Status",
       x = "Diet Variables",
       y = "Correlation Coefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate axis labels
  scale_x_discrete(labels = function(x) abbreviate(x, minlength = 5))

```
```{r}
#install.packages("ggpubr")
library(ggpubr)
library(reshape2)
diet_var <- df[1:ncol(df)]

diet_var$CaseControl[diet_var$CaseControl==1] <- 'AP'
diet_var$CaseControl[diet_var$CaseControl==0] <- 'Control'
fit1 <- lm(df$Verrucomicrobia~diet_var$CaseControl,data = diet_var)

Planctomycetes
Synergistetes


ggplot(data = diet_var, aes(x = CaseControl,y = Synergistetes,col = CaseControl))+
  geom_boxplot()+ stat_compare_means(method = 'wilcox.test',label =  "p.signif", label.x = 1.5)+
  stat_compare_means(method = 'wilcox.test', label.x = 1.5)+
  ylab("Synergistetes")+
  xlab("")+
  #scale_x_discrete(labels = c("Control", "AP"))+
  scale_color_manual(values=c("blue","pink"))+
  geom_jitter(width=0.15)+
  # scale_fill_manual(values=c("#66E91A", "#F29279", "#F50B0B"))+
  # theme_classic()+
  theme(axis.text.y = element_text(size = 10,face = "bold",color = "Black"),
        axis.title.y = element_text(size = 15,face = "bold"),
        axis.text.x = element_text(size = 15,face = "bold",color = "Black"),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  theme(legend.position = "none")

```






```{r}
p_values <- sapply(df2[, 2:ncol(df2)], function(x) {
  lm_model <- lm(x~ CaseControl, data = df2)
  summary(lm_model)$coefficients[2, "Pr(>|t|)"]
})

# Create a dataframe with variable names and their corresponding p-values
p_values_df <- data.frame(Variable = names(p_values), P_Value = p_values)

# Sort variables based on their p-values in ascending order
p_values_df <- p_values_df[order(p_values_df$P_Value), ]

# Select the topmost variables (e.g., top 10)
top_vars <- p_values_df[1:31, ]

# Plot the topmost variables
ggplot(top_vars, aes(x = reorder(Variable, P_Value), y = P_Value)) +
  geom_bar(stat = "identity") +
  labs(title = "Top Variables by P-Value", x = "Variable", y = "P-Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

```
```{r}
install.packages("bama")
library(readr)
#df <- read_csv("C:/Users/Palak/OneDrive/Desktop/Palak Thesis/Diet-microbiome relation/finalDiet1.csv")
df1 <- read_csv("C:/Users/Palak/OneDrive/Desktop/Palak Thesis/Diet-microbiome relation/finalMicrobiome1.csv")
df1 <- read_csv("C:/Users/Palak/OneDrive/Desktop/Palak Thesis/Diet-microbiome relation/microbiome_Phylum1.csv")
dim(abundance_data)
abundance_data<-df1[, 8:23]

abundance_data<- abundance_data[, colMeans(abundance_data) > 0.01]
library(bama)
for (i in vector) {
  
}
Y <- df1$CaseControl
A <- df1$tgrain
# grab the mediators from the example data.frame
M <- as.matrix(abundance_data)
# We just# We just# We just include the intercept term in this example as we have no covariates
C1 <- as.matrix(df1[2:5])
C2 <- as.matrix(df1[2:5])
beta.m <- rep(0, 1000)
alpha.a <- rep(0, 1000)
out <- bama(Y = Y, A = A, M = M, C1 = C1,C2 = C2,  method = "BSLMM", seed = 1234,
            burnin = 1000, ndraws = 1100, weights = NULL, inits = NULL,
            control = list(k = 2, lm0 = 1e-04, lm1 = 1, lma1 = 1, l = 1))
# The package includes a function to summarise output from 'bama'

summary <- summary(out)
write.xlsx(summary, "names1.xlsx")

```

```{r}
library(readr)

library(mediation)
finalPhylum <- read_csv("C:/Users/Palak/OneDrive/Desktop/Palak Thesis/Diet-microbiome relation/microbiome_Phylum.csv")
library(readr)
DietData_Imputed1 <- read_csv("C:/Users/Palak/OneDrive/Desktop/Palak Thesis/Diet-microbiome relation/DietData_Imputed1.csv")
df<-data.frame(DietData_Imputed1)
c<-df$CaseControl
cdf<-cbind(finalPhylum,exp_var)
cdf<-cbind(c,cdf)
cdf <- cbind(covariates, cdf)



library(readr)
df <- read_csv("C:/Users/Palak/OneDrive/Desktop/Palak Thesis/Diet-microbiome relation/finalMicrobiome1.csv")
abundance_data<-df[, 14:324]

df<- abundance_data[, colMeans(abundance_data) > 0.03]
combined_df <- cbind(df, exp_var)
c<-df$CaseControl
combined_df <- cbind(c, combined_df)
covariates<-df[, 2:5]
combined_df <- cbind(covariates, combined_df)




for (j in names(combined_df[,268:279])) {
  evar<-df[[j]]
  wb <- createWorkbook()
  expvar <- list()
  for (i in names(combined_df[, 2:267])) {
    dvar<-df[[i]]
    model.M <- lm(dvar ~ evar, data = combined_df)
    model.Y <- lm(c ~ dvar + evar, data = combined_df)
    
    result <- mediate(model.M, model.Y, treat = "evar", mediator = "dvar", boot = TRUE, sims = 100)
    #expvar[[i]] <-  summary(result)
    result_summary <- summary(result)
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
    
    # Create a dataframe with the extracted information and the variable name
    expvar[[i]] <- data.frame(
      Variable_Name = i,
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
    
    addWorksheet(wb, sheetName = )
    writeData(wb, sheet = j, x = expvar)
  }
}


library(openxlsx)
file_name <- "Exposure_variables.xlsx"
write.xlsx(exp_var, file_name)
```
```{r}
library(openxlsx)

# Create a workbook
wb <- createWorkbook()

# Loop through each column in combined_df[, 268:279]
for (j in names(combined_df[, 272:283])) {
  evar <- combined_df[[j]]
  
  # Create a list to store results
  expvar <- list()
  
  # Loop through each column in combined_df[, 2:267]
  for (i in names(combined_df[, 6:271])) {
    dvar <- combined_df[[i]]
    
    # Perform mediation analysis
    model.M <- lm(dvar ~ evar, data = combined_df)
    model.Y <- lm(c ~ dvar + evar+Age+Race+Gender+BMI, data = combined_df)
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
      Variable_Name = i,
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
  
  # Create a new worksheet with a specific name
  addWorksheet(wb, sheetName = j)
  
  # Write data to the worksheet
  writeData(wb, sheet = j, x = do.call(rbind, expvar))
}

# Save the workbook
saveWorkbook(wb, "mediation_results.xlsx")

```


```{r}
# Save the Excel workbook
saveWorkbook(wb, "mediation_results.xlsx", overwrite = TRUE)
results_df <- do.call(rbind, expvar)
class(results_df$ACME_P_Value)

correlation_df <- correlation_df[order(correlation_df$P_Value), ]
results_df <-results_df[order(results_df$ACME_P_Value), ]
significant_df <- subset(results_df, ACME_P_Value < 0.2)
significant_df 


subset_vars <- as.character(significant_df$Variable_Name)
exp_var <- df[,subset_vars]

library(openxlsx)

file_name <- "Exposure_variables.xlsx"

# Write the dataframe to an Excel file
write.xlsx(exp_var, file_name)
```

```{r}
# Install and load the hilama package
#install.packages("hilama")
n <- 47
p <- 16
q <- 12
M_multi=cdf[,3:18]
is.parallel <- TRUE  # Define is.parallel and set it to TRUE for parallel computation

score_all_ls <- esti_score_all(X = cdf[, 19:30], rho = 0.5, is.parallel = is.parallel, core_num = 4)
result_Theta_ls <- list()

for (j in 1:q) {
  Mj <- M_multi[, j]
  coef_hat_df <- infer_ddlasso(X = X, Mj, rho = 0.5, is.parallel = is.parallel, core_num = core_num,
                               score_ls = score_all_ls, is.adap = FALSE)
  result_Theta_ls[[j]] <- coef_hat_df
}


# Generate random data for X, M, and Y
X <- cdf[,19:30]
M <- cdf[,3:18]
Y <- cdf[4]
dim(Y)
# Run mediation analysis
result <- hilama(X = X, M_multi = M, Y = Y)

# Print the results
print(result)

```