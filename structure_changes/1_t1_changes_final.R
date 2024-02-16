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



########################## 1) t1 ~ ges_time ##########################
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/structure_changes/data/t1_all.csv")
# exclude non-pregnant participants 10014, 10016, and 10018
# exclude 10004_scan1 as 10004_scan1 was collected by t1/t2 not QUALAS
dat_total_subset <- dat_total[!(dat_total$sub_id=="10014" | dat_total$sub_id=="10016" 
                                | dat_total$sub_id=="10018" | dat_total$scan_id=="10004_scan1"),]
dat_total_subset$sub_id <- as.factor(dat_total_subset$sub_id)

dat_total_subset$ges_time_z <- scale(dat_total_subset$ges_time, center = T, scale = T)
dat_total_subset$eTIV_z <- scale(dat_total_subset$eTIV_x, center = T, scale = T)
dat_total_subset$age_z <- scale(dat_total_subset$age, center = T, scale = T)
dat_total_subset$num_pri_preg_z <- scale(dat_total_subset$num_pri_preg, center = T, scale = T)
dat_total_subset$pial_sum <- rowSums(subset(dat_total_subset, 
                                            select = c(lh_pial, rh_pial)), 
                                     na.rm = TRUE)
dat_total_subset$thickness_mean <- rowMeans(subset(dat_total_subset, 
                                                   select = c(lh_MeanThickness_thickness, 
                                                              rh_MeanThickness_thickness)), 
                                            na.rm = TRUE)

# y variables of interest
y_interest <- c('BrainSegVolNotVent_x', 'TotalGrayVol','SubCortGrayVol', 
                'CerebralWhiteMatterVol', 'CortexVol', 'pial_sum', 'thickness_mean')

# covariates: eTIV_x and age
df_result_ges <- data.frame(ROI1=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(), 
                            prec_change=double(), 
                            prec_change_rate=double(), stringsAsFactors=FALSE)
for (i in y_interest) {
  f <- formula(paste(i, " ~ ges_time + eTIV_x + age"))
  fit_tmp <- lme(f, random = ~ 1 | sub_id, 
                 data=dat_total_subset, na.action=na.omit)
  print(summary(fit_tmp))
  sum_tmp <- summary(fit_tmp)
  tmp_p <- sum_tmp$tTable["ges_time", "p-value"]
  b_tmp <- unname(sum_tmp$coefficients$fixed["ges_time"])
  # standardized beta and intervals
  dat_total_subset$y_z_tmpi <- scale(dat_total_subset[i], center = TRUE, scale = TRUE)
  fit_z_tmp <- lme(y_z_tmpi ~ ges_time_z + eTIV_z + age_z, 
                   random = ~ 1 | sub_id, 
                   data=dat_total_subset, na.action=na.omit)
  print(summary(fit_z_tmp))
  sum_tmp_z <- summary(fit_z_tmp)
  beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
  beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
  beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
  
  # append
  df_result_ges[nrow(df_result_ges) + 1,] = c(i, unname(tmp_p), b_tmp, 
                                              unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp),
                                              perc_c, perc_c_rate)
}

write.csv(df_result_ges, "df_result_ges_t1.csv", 
          row.names=F)


########################## 2) leave one out ##########################
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/structure_changes/data/t1_all.csv")
dat_total_subset <- dat_total[!(dat_total$sub_id=="10014" | dat_total$sub_id=="10016" 
                                | dat_total$sub_id=="10018" | dat_total$scan_id=="10004_scan1"),]
dat_total_subset$sub_id <- as.factor(dat_total_subset$sub_id)
dat_total_subset$pial_sum <- rowSums(subset(dat_total_subset, 
                                            select = c(lh_pial, rh_pial)), 
                                     na.rm = TRUE)
dat_total_subset$thickness_mean <- rowMeans(subset(dat_total_subset, 
                                                   select = c(lh_MeanThickness_thickness, 
                                                              rh_MeanThickness_thickness)), 
                                            na.rm = TRUE)

unique_subjects <- unique(dat_total_subset$sub_id)

# Loop through each unique subject
for (subject in unique_subjects) {
  # Subset the data to exclude the current subject
  subset_data <- dat_total_subset[dat_total_subset$sub_id != subject, ]
  
  subset_data$ges_time_z <- scale(subset_data$ges_time, center = T, scale = T)
  subset_data$eTIV_z <- scale(subset_data$eTIV_x, center = T, scale = T)
  subset_data$age_z <- scale(subset_data$age, center = T, scale = T)
  
  # y variables of interest
  y_interest <- c('BrainSegVolNotVent_x', 'TotalGrayVol','SubCortGrayVol', 
                  'CerebralWhiteMatterVol', 'CortexVol', 'pial_sum', 'thickness_mean')
  
  # covariate: eTIV_x and age
  df_result_ges <- data.frame(ROI1=character(),
                              p=double(),
                              b=double(),
                              beta=double(),
                              beta_lower=double(),
                              beta_upper=double(), 
                              prec_change=double(), 
                              prec_change_rate=double(), stringsAsFactors=FALSE)
  for (i in y_interest) {
    f <- formula(paste(i, " ~ ges_time + eTIV_x + age"))
    fit_tmp <- lme(f, random = ~ 1 | sub_id, 
                   data=subset_data, na.action=na.omit)
    print(summary(fit_tmp))
    sum_tmp <- summary(fit_tmp)
    tmp_p <- sum_tmp$tTable["ges_time", "p-value"]
    b_tmp <- unname(sum_tmp$coefficients$fixed["ges_time"])
    # standardized beta and intervals
    subset_data$y_z_tmpi <- scale(subset_data[i], center = TRUE, scale = TRUE)
    fit_z_tmp <- lme(y_z_tmpi ~ ges_time_z + eTIV_z + age_z, 
                     random = ~ 1 | sub_id, 
                     data=subset_data, na.action=na.omit)
    print(summary(fit_z_tmp))
    sum_tmp_z <- summary(fit_z_tmp)
    beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
    beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
    beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
    # percentage change
    ges_min <- min(subset_data$ges_time, na.rm = TRUE)
    ges_max <- max(subset_data$ges_time, na.rm = TRUE)
    new_data_min <- data.frame(ges_time = ges_min, 
                               eTIV_x = mean(subset_data$eTIV_x, na.rm = TRUE), 
                               age = mean(subset_data$age, na.rm = TRUE))
    pred_min <- predict(fit_tmp, new_data_min, level=0)
    new_data_max <- data.frame(ges_time = ges_max, 
                               eTIV_x = mean(subset_data$eTIV_x, na.rm = TRUE), 
                               age = mean(subset_data$age, na.rm = TRUE))
    pred_max <- predict(fit_tmp, new_data_max, level=0)
    perc_c <- (pred_max - pred_min) / pred_min
    perc_c_rate <- (pred_max - pred_min) / (ges_max - ges_min)
    
    # append
    df_result_ges[nrow(df_result_ges) + 1,] = c(i, unname(tmp_p), b_tmp, 
                                                unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp),
                                                perc_c, perc_c_rate)
  }
  
  filename <- paste0("df_result_ges_t1_loo_rm_", subject, ".csv")
  write.csv(df_result_ges, 
            filename, 
            row.names=F)
}


########################## 3) aseg volumes ##########################
library(tidyr)
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/structure_changes/data/t1_all.csv")

dat_total_subset <- dat_total[!(dat_total$sub_id=="10014" | dat_total$sub_id=="10016" 
                                | dat_total$sub_id=="10018" | dat_total$scan_id=="10004_scan1"
                                | dat_total$scan_id=="10003_scan3"),]
dat_total_subset$sub_id <- as.factor(dat_total_subset$sub_id)

dat_total_subset$ges_time_z <- scale(dat_total_subset$ges_time, center = T, scale = T)
dat_total_subset$eTIV_z <- scale(dat_total_subset$eTIV_x, center = T, scale = T)
dat_total_subset$age_z <- scale(dat_total_subset$age, center = T, scale = T)
dat_total_subset$num_pri_preg_z <- scale(dat_total_subset$num_pri_preg, center = T, scale = T)

dat_total_subset$scan <- unlist(lapply(dat_total_subset$scan_id, function(x) strsplit(x, '_')[[1]][2]))
list_subregions <- list()
list_subregions_long <- list()
# from "bankssts_volume" to "insula_volume
for (sub in 1:34) {
  print(sub)
  sub_lh <- colnames(dat_total_subset)[(136 + sub)]
  sub_rh <- gsub("lh", "rh", sub_lh)
  aa <- dat_total_subset[, c(sub_lh, sub_rh, "eTIV_x", "eTIV_z", 
                             "ges_time", "ges_time_z", "age", "age_z", "sub_id", 
                             "num_pri_preg", "num_pri_preg_z"), drop=F]
  aa_long <- aa %>% pivot_longer(
    cols = c(1, 2),
    names_to = c('hemi', '.value'),
    names_pattern = "(.*)_(.*)_volume") %>%
    convert_as_factor(sub_id, hemi)
  list_subregions[[sub]] <- aa
  list_subregions_long[[sub]] <- na.omit(aa_long)
}

# covariate with eTIV and age
df_result_ges <- data.frame(ROI=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(), stringsAsFactors=FALSE)
for (i in 1:length(list_subregions_long)) {
  data <- list_subregions_long[[i]]
  data_dv <- names(data)[ncol(data)]
  print(data_dv)
  fit_func <- formula(paste(data_dv, "~ hemi + eTIV_x + age + ges_time"))
  fit_tmp <- lme(fit_func, random = ~ 1 | sub_id, data=data, na.action=na.omit)
  print(summary(fit_tmp))
  sum_tmp <- summary(fit_tmp)
  tmp_p <- sum_tmp$tTable["ges_time", "p-value"]
  b_tmp <- unname(sum_tmp$coefficients$fixed['ges_time'])
  # standardized beta and intervals
  data$y_z_tmpi <- scale(data[[data_dv]], center = TRUE, scale = TRUE)
  fit_z_tmp <- lme(y_z_tmpi ~ hemi + eTIV_z + age_z + ges_time_z, random = ~ 1 | sub_id,
                   data=data, na.action=na.omit)
  print(summary(fit_z_tmp))
  beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
  beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
  beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
  print(intervals(fit_z_tmp, which=c("fixed")))
  df_result_ges[nrow(df_result_ges) + 1,] = c(data_dv, unname(tmp_p), b_tmp,
                                              unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
}

write.csv(df_result_ges, "/Users/yanbinniu/Projects/NCAP/scripts/analysis_t1/results/df_result_ges_t1_aseg_stdCov.csv", row.names=F)
