
## Data Processing -------------------------------------------------------------
# Load the data
og_data <- read.csv("brca_tcga_pub2015_clinical_data_original.csv")

colnames(org_data)

# Subset the original data for the study
data <- og_data[, c("Patient.ID", "Sample.ID", "Sex", "Ethnicity.Category", 
                    "Race.Category", "Diagnosis.Age", "Year.Cancer.Initial.Diagnosis", 
                    "Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code", 
                    "Cancer.Type.Detailed", "ICD.10.Classification", 
                    "Days.to.Sample.Collection.", "Days.to.Last.Followup", 
                    "Overall.Survival..Months.", "Overall.Survival.Status")]

# Rename columns
names(data) <- c("Patient ID", "Sample ID", "Gender", "Ethnicity", "Race", 
                 "Age at Diagnosis", "Year of Initial Diagnosis", "Cancer Stage", 
                 "Cancer Type", "ICD-10 Classification", "Days to Sample Collection", 
                 "Days to Last Follow-up", "Overall Survival (Months)", 
                 "Overall Survival Status")

# Remove the variables not used in the study
data <- data[, !(names(data) %in% c("Patient ID", "Sample ID", "ICD-10 Classification", 
                                    "Days to Sample Collection", "Year of Initial Diagnosis", 
                                    "Days to Last Follow-up"))]

# Convert "LIVING" to 0 and "DECEASED" to 1
data$`Overall Survival Status` <- ifelse(data$`Overall Survival Status` == "0:LIVING", 0, 1)

# Remove rows with missing data, including empty strings
cdata <- data[rowSums(data == "" | is.na(data)) == 0, ]



## Exploratory data analysis (EDA) ---------------------------------------------
library(ggplot2)
library(dplyr)
library(gridExtra)

# Summary statistics
summary(cdata)

# Histogram of Survival Time (Months)
hist_survival <- ggplot(cdata, aes(x = `Overall Survival (Months)`)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Survival Time",
       x = "Survival Time (Months)",
       y = "Frequency")

# Histogram of Age_at_Diagnosis
hist_age <- ggplot(cdata, aes(x = `Age at Diagnosis`)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Age at Diagnosis",
       x = "Age at Diagnosis",
       y = "Frequency")

# Display both histograms in the same frame
grid.arrange(hist_survival, hist_age, ncol = 2)

# Box plot of Overall_Survival_Months
boxplot_survival <- ggplot(cdata, aes(y = `Overall Survival (Months)`)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Boxplot of Overall Survival Months",
       x = "",
       y = "Overall Survival Months")

# Box plot of Age_at_Diagnosis
boxplot_age <- ggplot(cdata, aes(y = `Age at Diagnosis`)) +
  geom_boxplot(fill = "skyblue", color = "black") +
  labs(title = "Boxplot of Age at Diagnosis",
       x = "",
       y = "Age at Diagnosis")

# Display both box plots in the same frame
grid.arrange(boxplot_survival, boxplot_age, ncol = 2)

# Bar plot of Cancer_Stage
bar_cancer_stage <- ggplot(cdata, aes(x = `Cancer Stage`)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Cancer Stage",
       x = "Cancer Stage",
       y = "Count")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Bar plot for Gender
bar_gender <- ggplot(cdata, aes(x = Gender)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Gender",
       x = "Gender",
       y = "Count")

# Bar plot for Ethnicity
bar_ethnicity <- ggplot(cdata, aes(x = Ethnicity)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Ethnicity",
       x = "Ethnicity",
       y = "Count")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1,size = 8))

# Bar plot for Race
bar_race <- ggplot(cdata, aes(x = Race)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Distribution of Race",
       x = "Race",
       y = "Count")+
  theme(axis.text.x = element_text(angle = 20, hjust = 1,size = 8))

grid.arrange(bar_cancer_stage, bar_gender, ncol = 2)
grid.arrange(bar_ethnicity, bar_race, ncol = 2)

# Scatter plot of Age at Diagnosis vs. Overall Survival (Months)
plot(cdata$`Age at Diagnosis`, cdata$`Overall Survival (Months)`, 
     main = "Age at Diagnosis vs. Overall Survival Time", 
     xlab = "Age at Diagnosis", ylab = "Overall Survival (Months)")


# Correlation matrix
correlation_matrix <- cor(cdata[, c("Age at Diagnosis", "Overall Survival (Months)")])
correlation_matrix


## Kaplan-Meier survival curves -----------------------------------------------
library(survival)

# Create the survival object
surv_object <- Surv(time = cdata$`Overall Survival (Months)`, event = cdata$`Overall Survival Status`)
# Fit the Kaplan-Meier survival curve
km_fit <- survfit(surv_object ~ 1)
# Table of Kaplan-Meier estimates
summary(km_fit)
# Plot the Kaplan-Meier curve
plot(km_fit, main = "Kaplan-Meier Curve", xlab = "Time (Months)", ylab = "Survival Probability")

# Survival curves for different Cancer Types
surv_fit <- survfit(surv_object ~ cdata$`Cancer Type`)
summary(surv_fit)

# Plot Kaplan-Meier Survival Curves by Cancer Type with colors
# Get unique levels of Cancer Type
cancer_types <- unique(cdata$`Cancer Type`)
# Plot Kaplan-Meier Survival Curves
plot(surv_fit, main = "Kaplan-Meier Survival Curves by Cancer Type", xlab = "Time (Months)", 
     ylab = "Survival Probability", col = c("slateblue", "salmon", "forestgreen", "goldenrod"))
legend("bottomright", fill = c("slateblue", "salmon", "forestgreen", "goldenrod"), 
       legend = cancer_types, ncol = 1, cex = 0.5, x.intersp = 0.1, text.width = 60)


## Log rank test and Cox hazard models------------------------------------------
# logrank test to compare the survival between different cancer types
surv_diff <- survdiff(surv_object ~ cdata$`Cancer Type`)
surv_diff

# Cox proportional hazards model to assess the relationship between survival time and demographics
cox_model_1 <- coxph(Surv(`Overall Survival (Months)`, `Overall Survival Status`) ~ Gender + Ethnicity + Race, data = cdata)
summary(cox_model_1)

# Create survival curve
surv_curve_1 <- survfit(cox_model_1)
# Plot the survival curve
plot(surv_curve_1, main = "Cox Proportional Hazards Survival Curve", xlab = "Time (Months)", ylab = "Survival Probability")

# Cox proportional hazards model to assess the relationship between survival time, age at diagnosis, and cancer stage
cox_model_2 <- coxph(Surv(`Overall Survival (Months)`, `Overall Survival Status`) ~ `Age at Diagnosis` + `Cancer Stage`, data = cdata)
summary(cox_model_2)

# Create survival curve
surv_curve_2 <- survfit(cox_model_2)
# Plot the survival curve
plot(surv_curve_2, main = "Cox Proportional Hazards Survival Curve", xlab = "Time (Months)", ylab = "Survival Probability")


# Cox proportional hazards model to assess the relationship between survival time and cancer TYpes
cox_model_3 <- coxph(Surv(`Overall Survival (Months)`, `Overall Survival Status`) ~ `Cancer Type`, data = cdata)
summary(cox_model_3)

# Create survival curve
surv_curve_3 <- survfit(cox_model_3)
# Plot the survival curve
plot(surv_curve_3, main = "Cox Proportional Hazards Survival Curve", xlab = "Time (Months)", ylab = "Survival Probability")




## Encoding variables to build a model -----------------------------------------
# Gender
cdata$Gender_Label <- ifelse(cdata$Gender == "Female", 0, 1)
# Ethnicity
cdata$Ethnicity_Label <- ifelse(cdata$Ethnicity == "HISPANIC OR LATINO", 1, 0)
# Race
cdata$Race_Label <- ifelse(cdata$Race == "AMERICAN INDIAN OR ALASKA NATIVE", 1,
                           ifelse(cdata$Race == "ASIAN", 2,
                                  ifelse(cdata$Race == "BLACK OR AFRICAN AMERICAN", 3,
                                         ifelse(cdata$Race == "WHITE", 4, 0))))
# Cancer Stage
cdata$Cancer_Stage_Label <- ifelse(cdata$`Cancer Stage` == "Stage I", 1,
                                   ifelse(cdata$`Cancer Stage` == "Stage IA", 2,
                                          ifelse(cdata$`Cancer Stage` == "Stage IB", 3,
                                                 ifelse(cdata$`Cancer Stage` == "Stage II", 4,
                                                        ifelse(cdata$`Cancer Stage` == "Stage IIA", 5,
                                                               ifelse(cdata$`Cancer Stage` == "Stage IIB", 6,
                                                                      ifelse(cdata$`Cancer Stage` == "Stage III", 7,
                                                                             ifelse(cdata$`Cancer Stage` == "Stage IIIA", 8,
                                                                                    ifelse(cdata$`Cancer Stage` == "Stage IIIB", 9,
                                                                                           ifelse(cdata$`Cancer Stage` == "Stage IIIC", 10,
                                                                                                  ifelse(cdata$`Cancer Stage` == "Stage IV", 11,
                                                                                                         ifelse(cdata$`Cancer Stage` == "Stage X", 12, 0))))))))))))
# Cancer Type
cdata$Cancer_Type_Label <- ifelse(cdata$`Cancer Type` == "Breast Invasive Ductal Carcinoma", 1,
                                  ifelse(cdata$`Cancer Type` == "Breast Invasive Lobular Carcinoma", 2,
                                         ifelse(cdata$`Cancer Type` == "Breast Mixed Ductal and Lobular Carcinoma", 3,
                                                ifelse(data$`Cancer Type` == "Invasive Breast Carcinoma", 4, 0))))


# Create a subset of data including predictor variables and outcome variable
model_data <- subset(cdata, select = c("Gender_Label", "Ethnicity_Label", "Race_Label", 
                                       "Cancer_Stage_Label", "Cancer_Type_Label", 
                                       "Age at Diagnosis", "Overall Survival (Months)", 
                                       "Overall Survival Status"))

## Correlation analyis ---------------------------------------------------------
correlation_matrix <- cor(model_data)
correlation_matrix

# Heatmap of correlation matrix
library(corrplot)
corrplot(correlation_matrix, method = "color", type = "upper", 
         order = "hclust", tl.cex = 0.7, 
         col = colorRampPalette(c("red", "white", "blue"))(100))



# Model Building ---------------------------------------------------------------
# Train-test split
set.seed(2024)
train_index <- sample(1:nrow(model_data), 0.8 * nrow(model_data))
train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]




## Model considering all the variables -----------------------------------------
# Model training
model <- glm(`Overall Survival Status` ~ ., data = train_data, family = binomial)

# Model Evaluation
# Predictions on test set
predictions <- predict(model, newdata = test_data, type = "response")
# Convert predicted probabilities to class labels
predicted_classes <- ifelse(predictions > 0.5, 1, 0)

# Model performance metrics
accuracy <- mean(predicted_classes == test_data$`Overall Survival Status`)
print(paste("Accuracy:", accuracy))
# AUC-ROC curve
library(pROC)
roc_curve <- roc(test_data$`Overall Survival Status`, predictions)
auc <- auc(roc_curve)
print(paste("AUC:", auc))

# Model Interpretation
summary(model)




## Model based on significance -------------------------------------------------
# Model training
model_sig <- glm(`Overall Survival Status` ~ Cancer_Stage_Label + `Age at Diagnosis` + `Overall Survival (Months)`, 
                 data = train_data, family = binomial)

# Model Evaluation
# Predictions on test set
predictions_sig <- predict(model_sig, newdata = test_data, type = "response")
# Convert predicted probabilities to class labels
predicted_classes_sig <- ifelse(predictions_sig > 0.5, 1, 0)


# Model performance metrics
accuracy_sig <- mean(predicted_classes_sig == test_data$`Overall Survival Status`)
print(paste("Accuracy:", accuracy_sig))
# AUC-ROC curve
roc_curve_sig <- roc(test_data$`Overall Survival Status`, predictions_sig)
auc_sig <- auc(roc_curve_sig)
print(paste("AUC:", auc_sig))

# Model Interpretation
summary(model_sig)



