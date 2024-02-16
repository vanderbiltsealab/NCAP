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


########################## hormone ~ ges_time ##########################
# Note that hormone has 24 observations after removing 10001_scan1
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/structure_changes/data/t1_all.csv")
dat_total_subset_hor <- dat_total[!(dat_total$sub_id=="10014" | dat_total$sub_id=="10016"
                                    | dat_total$sub_id=="10018" | dat_total$scan_id=="10001_scan1"),]
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


########################## hair hormone ~ ges_time ##########################
# Note that hair hormone does not need to remove 10001_scan1
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/structure_changes/data/t1_all.csv")
dat_total_subset_hor <- dat_total[!(dat_total$sub_id=="10014" | dat_total$sub_id=="10016"
                                    | dat_total$sub_id=="10018" | dat_total$sub_id=="10001"),]
dat_total_subset_hor$sub_id <- as.factor(dat_total_subset_hor$sub_id)

dat_total_subset_hor$ges_time_z <- scale(dat_total_subset_hor$ges_time, center = T, scale = T)
dat_total_subset_hor$age_z <- scale(dat_total_subset_hor$age, center = T, scale = T)
dat_total_subset_hor$cortisol_win_centered_z <- scale(dat_total_subset_hor$cortisol_win_centered, center = T, scale = T)
dat_total_subset_hor$cortisone_win_centered_z <- scale(dat_total_subset_hor$cortisone_win_centered, center = T, scale = T)
dat_total_subset_hor$progesterone_win_centered_z <- scale(dat_total_subset_hor$progesterone_win_centered, center = T, scale = T)
dat_total_subset_hor$segment_sample_use_z <- scale(dat_total_subset_hor$segment_sample_use, center = T, scale = T)

# y variables of interest
hormone_interest <- c('cortisol_win_centered', 'cortisone_win_centered', 'progesterone_win_centered')

# covariates: age
# note that controlling age will not converge
df_result_ges <- data.frame(ROI1=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(), stringsAsFactors=FALSE)
for (i in hormone_interest) {
  f <- formula(paste(i, " ~ ges_time + segment_sample_use"))
  fit_tmp <- lme(f, random = ~ segment_sample_use | sub_id, 
                 data=dat_total_subset_hor, na.action=na.omit)
  print(summary(fit_tmp))
  sum_tmp <- summary(fit_tmp)
  tmp_p <- sum_tmp$tTable["ges_time", "p-value"]
  b_tmp <- unname(sum_tmp$coefficients$fixed["ges_time"])
  # standardized beta and intervals
  dat_total_subset_hor$y_z_tmpi <- scale(dat_total_subset_hor[i], center = TRUE, scale = TRUE)
  fit_z_tmp <- lme(y_z_tmpi ~ ges_time_z + segment_sample_use_z, 
                   random = ~ segment_sample_use_z | sub_id, 
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

write.csv(df_result_ges, "df_result_ges_hormone_hair_stdCov.csv", 
          row.names=F)

