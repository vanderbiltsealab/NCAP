library(readxl)
library(nlme)
library(lme4)
library(ggplot2)
library(sjmisc)
library(sjstats)
library(parameters)
library(r2glmm)
library(tibble)
library(dplyr)
library(tidyr)
library(rstatix)

########################## 1) csf change ##########################
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/csfmap.csv")

# remove 10004_scan1 which did not have qalas data
dat_total_subset <- dat_total[!(dat_total$fsid=='10004_scan1'),]
dat_total_subset$sub_id <- as.factor(dat_total_subset$sub_id)

# control for TICV
dat_total_subset$ges_time_z <- scale(dat_total_subset$ges_time , center = TRUE, scale = TRUE)
dat_total_subset$age_z <- scale(dat_total_subset$age , center = TRUE, scale = TRUE)
dat_total_subset$EstimatedTICV_z <- scale(dat_total_subset$EstimatedTICV , center = TRUE, scale = TRUE)

# control for TICV and age
df_result_ges <- data.frame(ROI=character(),
                            p=double(),
                            b=double(),
                            interc=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(),stringsAsFactors=FALSE)
for (i in names(dat_total_subset)[3:10]) {
  print(i)
  fit_func <- formula(paste(i, " ~ ges_time + age + EstimatedTICV"))
  fit_tmp <- lme(fit_func, 
                 random = ~ 1 | sub_id, 
                 data=dat_total_subset, 
                 na.action=na.omit)
  print(summary(fit_tmp))
  sum_tmp <- summary(fit_tmp)
  tmp_p <- sum_tmp$tTable["ges_time", "p-value"]
  tmp_interc <- sum_tmp$tTable["(Intercept)", "Value"]
  b_tmp <- unname(sum_tmp$coefficients$fixed['ges_time'])
  
  # standardized beta and intervals
  idx_tmp <- dat_total_subset[,i]
  y_z_tmp <- scale(idx_tmp, center = TRUE, scale = TRUE)
  fit_z_tmp <- lme(y_z_tmp ~ ges_time_z+age_z+EstimatedTICV_z, 
                   random = ~ 1 | sub_id, 
                   data=dat_total_subset, 
                   na.action=na.omit)
  beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
  beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
  beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
  print(intervals(fit_z_tmp, which=c("fixed")))
  
  # append
  df_result_ges[nrow(df_result_ges) + 1,] = c(i, unname(tmp_p), b_tmp, tmp_interc,
                                              unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
}


########################## 2) leave one out ##########################
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/csfmap.csv")

# remove 10004_scan1 which did not have qalas data
dat_total_subset <- dat_total[!(dat_total$fsid=='10004_scan1'),]
dat_total_subset$sub_id <- as.factor(dat_total_subset$sub_id)

unique_subjects <- unique(dat_total_subset$sub_id)

df_result_ges <- data.frame(ROI1=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(),  stringsAsFactors=FALSE)
# Loop through each unique subject
for (subject in unique_subjects) {
  # Subset the data to exclude the current subject
  subset_data <- dat_total_subset[dat_total_subset$sub_id != subject, ]
  
  subset_data$ges_time_z <- scale(subset_data$ges_time, center = T, scale = T)
  subset_data$EstimatedTICV_z <- scale(subset_data$EstimatedTICV, center = T, scale = T)
  subset_data$age_z <- scale(subset_data$age, center = T, scale = T)
  subset_data$csf_80_z <- scale(subset_data$csf_80, center = T, scale = T)
  
  fit_func <- formula(paste("csf_80", " ~ ges_time + age + EstimatedTICV"))
  fit_tmp <- lme(fit_func, 
                 random = ~ 1 | sub_id, 
                 data=subset_data, 
                 na.action=na.omit)
  print(summary(fit_tmp))
  sum_tmp <- summary(fit_tmp)
  tmp_p <- sum_tmp$tTable["ges_time", "p-value"]
  tmp_interc <- sum_tmp$tTable["(Intercept)", "Value"]
  b_tmp <- unname(sum_tmp$coefficients$fixed['ges_time'])
  
  # standardized beta and intervals
  fit_z_tmp <- lme(csf_80_z ~ ges_time_z + age_z + EstimatedTICV_z, 
                   random = ~ 1 | sub_id, 
                   data=subset_data, 
                   na.action=na.omit)
  beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
  beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
  beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
  print(intervals(fit_z_tmp, which=c("fixed")))
  
  # append
  df_result_ges[nrow(df_result_ges) + 1,] = c(subject, unname(tmp_p), b_tmp, 
                                              unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
}
