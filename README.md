# Survival Analysis in Breast Cancer: Insights from TCGA Data

## Overview
This project analyzes breast cancer patient data from the TCGA (The Cancer Genome Atlas) to study factors affecting overall survival. We performed exploratory data analysis, survival analysis (Kaplan-Meier and Cox regression), and built predictive logistic regression models to classify patient outcomes.

The dataset was sourced from [cBioPortal - BRCA TCGA Public 2015](https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pub2015).

---

## Data Description
The dataset contains clinical and demographic data for breast cancer patients, including:

- Patient ID, Sample ID, Sex, Ethnicity, Race
- Age at diagnosis, Year of initial diagnosis
- Cancer staging and type details
- Survival times and status
- Days to sample collection and follow-up

Columns were selected and renamed for clarity during preprocessing.

---

## Methodology
1. Data cleaning and selection of relevant features.
2. Exploratory Data Analysis (EDA) including correlation heatmaps and distribution plots.
3. Survival analysis using Kaplan-Meier curves and Cox proportional hazards models.
4. Predictive modeling with logistic regression to classify survival status.
5. Evaluation of models using ROC curves and AUC scores.

---

## Results
The data showed varied survival outcomes across cancer stages and types. Age at diagnosis and cancer stage were key predictors of survival. Logistic regression models achieved moderate accuracy (~77.6%) and AUC (~0.63), showing reasonable prediction of patient outcomes.

---
## Requirements

Make sure you have R installed (version 4.3 or higher recommended).

This project requires the following R packages:

ggplot2

dplyr

gridExtra

survival

corrplot

pROC

You can install them using:

```r
install.packages(c("ggplot2", "dplyr", "gridExtra", "survival", "corrplot", "pROC"))
```

---

## Usage
1. Clone the repository.
2. Update the data file path inside the script to point to your local copy of the dataset.
3. Run the script in R
         
          survival_modeling_tcga_breast_cancer.R

4. Check the console and plots generated during the run for results.

---

## Team Contributions

**Poojitha K:** Data cleaning, analysis, modeling, and visualization  
**Mahima Mahabaleshwar S:** Background research and dataset sourcing  
**Yesasvi Sai N:** Project methodology design  
**Muni Manasa V:** Results interpretation and conclusions  
**Tejaswini R:** Presentation and reporting

