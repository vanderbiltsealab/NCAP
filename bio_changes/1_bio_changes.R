library(readxl)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjmisc)
library(sjstats)
library(parameters)
library(r2glmm)
library(rstatix)


########################## 1) hormone ~ ges_time ##########################
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/t1_all.csv")
dat_total_subset_hor <- dat_total[!(dat_total$scan_id=="10001_scan1"),]
dat_total_subset_hor$sub_id <- as.factor(dat_total_subset_hor$sub_id)

dat_total_subset_hor$ges_time_z <- scale(dat_total_subset_hor$ges_time, center = T, scale = T)
dat_total_subset_hor$age_z <- scale(dat_total_subset_hor$age, center = T, scale = T)
dat_total_subset_hor$progesterone_z <- scale(dat_total_subset_hor$progesterone, center = T, scale = T)
dat_total_subset_hor$estradiol_z <- scale(dat_total_subset_hor$estradiol, center = T, scale = T)
dat_total_subset_hor$testosterone_z <- scale(dat_total_subset_hor$testosterone, center = T, scale = T)
dat_total_subset_hor$cortisol_z <- scale(dat_total_subset_hor$cortisol, center = T, scale = T)

# y variables of interest
hormone_interest <- c('progesterone', 'estradiol', 'testosterone', 'cortisol')

# covariates: eTIV_x and age
df_result_ges <- data.frame(ROI1=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(), stringsAsFactors=FALSE)
for (i in hormone_interest) {
  f <- formula(paste(i, " ~ ges_time + age"))
  fit_tmp <- lme(f, random = ~ 1 | sub_id, 
                 data=dat_total_subset_hor, na.action=na.omit)
  print(summary(fit_tmp))
  sum_tmp <- summary(fit_tmp)
  tmp_p <- sum_tmp$tTable["ges_time", "p-value"]
  b_tmp <- unname(sum_tmp$coefficients$fixed["ges_time"])
  # standardized beta and intervals
  dat_total_subset_hor$y_z_tmpi <- scale(dat_total_subset_hor[i], center = TRUE, scale = TRUE)
  fit_z_tmp <- lme(y_z_tmpi ~ ges_time_z + age_z, 
                   random = ~ 1 | sub_id, 
                   data=dat_total_subset_hor, na.action=na.omit)
  print(summary(fit_z_tmp))
  sum_tmp_z <- summary(fit_z_tmp)
  beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
  beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
  beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
  
  # append
  df_result_ges[nrow(df_result_ges) + 1,] = c(i, unname(tmp_p), b_tmp, 
                                              unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
}

write.csv(df_result_ges, "df_result_ges_hormone.csv", 
          row.names=F)


########################## 2) hair hormone ~ ges_time ##########################
dat_total_subset_hor <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/t1_all.csv")
dat_total_subset_hor$sub_id <- as.factor(dat_total_subset_hor$sub_id)

dat_total_subset_hor$ges_time_z <- scale(dat_total_subset_hor$ges_time, center = T, scale = T)
dat_total_subset_hor$age_z <- scale(dat_total_subset_hor$age, center = T, scale = T)
dat_total_subset_hor$cortisol_win_centered_z <- scale(dat_total_subset_hor$cortisol_win_centered, center = T, scale = T)
dat_total_subset_hor$cortisone_win_centered_z <- scale(dat_total_subset_hor$cortisone_win_centered, center = T, scale = T)
dat_total_subset_hor$progesterone_win_centered_z <- scale(dat_total_subset_hor$progesterone_win_centered, center = T, scale = T)
dat_total_subset_hor$distance <- dat_total_subset_hor$segment_sample_use - 1
dat_total_subset_hor$distance_z <- scale(dat_total_subset_hor$distance, center = T, scale = T)

# y variables of interest
hormone_interest <- c('cortisol_win_centered', 'cortisone_win_centered', 'progesterone_win_centered')

# covariate: age and distance 
# note that random distance will not converge
df_result_ges <- data.frame(ROI1=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(), stringsAsFactors=FALSE)
for (i in hormone_interest) {
  f <- formula(paste(i, " ~ ges_time + distance + age"))
  fit_tmp <- lme(f, random = ~ 1 | sub_id, 
                 data=dat_total_subset_hor, na.action=na.omit)
  print(summary(fit_tmp))
  sum_tmp <- summary(fit_tmp)
  tmp_p <- sum_tmp$tTable["ges_time", "p-value"]
  b_tmp <- unname(sum_tmp$coefficients$fixed["ges_time"])
  # standardized beta and intervals
  dat_total_subset_hor$y_z_tmpi <- scale(dat_total_subset_hor[i], center = TRUE, scale = TRUE)
  fit_z_tmp <- lme(y_z_tmpi ~ ges_time_z + distance_z + age_z, 
                   random = ~ 1 | sub_id, 
                   data=dat_total_subset_hor, na.action=na.omit)
  print(summary(fit_z_tmp))
  sum_tmp_z <- summary(fit_z_tmp)
  beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
  beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
  beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
  
  # append
  df_result_ges[nrow(df_result_ges) + 1,] = c(i, unname(tmp_p), b_tmp, 
                                              unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
}

write.csv(df_result_ges, "df_result_ges_hormone_hair.csv", 
          row.names=F)


########################## 3) immune ~ ges_time ##########################
library(psych)

data_t1 <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/t1_all.csv")

# log transform
data_t1$CRP_log <- log(data_t1$CRP)

# subset
data_t1$sub_id <- as.factor(data_t1$sub_id)
data_t1_sub <- subset(data_t1, ges_time>0)
# removing 10012_scan2 as an outlier.
data_t1_sub <- data_t1_sub[!((data_t1_sub$scan_id=="10012_scan2")),]
dat_total_subset <- data_t1_sub

# std
dat_total_subset$ges_time_z <- scale(dat_total_subset$ges_time, center = TRUE, scale = TRUE)
dat_total_subset$age_z <- scale(dat_total_subset$age, center = TRUE, scale = TRUE)

y_interest <- c("CRP_log")

df_result_ges <- data.frame(ROI1=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(), stringsAsFactors=FALSE)
for (i in y_interest) {
  f <- formula(paste(i, "~ ges_time + age"))
  fit_tmp <- lme(f, random = ~ 1 | sub_id, 
                 data=dat_total_subset, 
                 na.action=na.omit)
  print(summary(fit_tmp))
  sum_tmp <- summary(fit_tmp)
  tmp_p <- sum_tmp$tTable["ges_time", "p-value"]
  b_tmp <- unname(sum_tmp$coefficients$fixed['ges_time'])
  # standardized beta and intervals
  dat_total_subset$y_z_tmpi <- scale(dat_total_subset[i], center = TRUE, scale = TRUE)
  fit_z_tmp <- lme(y_z_tmpi ~ ges_time_z + age_z, 
                   random = ~ 1 | sub_id, 
                   data=dat_total_subset, 
                   na.action=na.omit)
  sum_tmp_z <- summary(fit_z_tmp)
  beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
  beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
  beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
  # append
  df_result_ges[nrow(df_result_ges) + 1,] = c(i, unname(tmp_p), b_tmp, 
                                              unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
}

write.csv(df_result_ges, "immune_change.csv", 
          row.names = FALSE)


