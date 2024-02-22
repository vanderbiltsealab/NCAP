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


############## 1) csf and hormone ##############
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/csfmap.csv")
# exclude 10001_scan1 and 10004_scan1
dat_total_subset <- dat_total[!(dat_total$fsid=='10004_scan1'
                                | dat_total$fsid=='10001_scan1'),]
dat_total_subset$sub_id <- as.factor(dat_total_subset$sub_id)

dat_total_subset$age_z <- scale(dat_total_subset$age , center = TRUE, scale = TRUE)
dat_total_subset$EstimatedTICV_z <- scale(dat_total_subset$EstimatedTICV , center = TRUE, scale = TRUE)

hormone_interest <- c('progesterone', 'estradiol', 'testosterone', 'cortisol')

for (j in hormone_interest) {
  
  # loop through all thresholds
  df_result_ges <- data.frame(ROI=character(),
                              p=double(),
                              b=double(),
                              interc=double(),
                              beta=double(),
                              beta_lower=double(),
                              beta_upper=double(), stringsAsFactors=FALSE)
  
  for (i in names(dat_total_subset)[3:9]){
    print(i)
    fit_func <- formula(paste(i, " ~ age + EstimatedTICV +", j))
    fit_tmp <- lme(fit_func, 
                   random = ~ 1 | sub_id, 
                   data=dat_total_subset, 
                   na.action=na.omit)
    print(summary(fit_tmp))
    sum_tmp <- summary(fit_tmp)
    tmp_p <- sum_tmp$tTable[j, "p-value"]
    tmp_interc <- sum_tmp$tTable["(Intercept)", "Value"]
    b_tmp <- unname(sum_tmp$coefficients$fixed[j])
    
    # standardized beta and intervals
    dat_total_subset$y_z_tmpi <- scale(dat_total_subset[i], center = TRUE, scale = TRUE)
    dat_total_subset$y_z_tmpj <- scale(dat_total_subset[j], center = TRUE, scale = TRUE)
    fit_z_tmp <- lme(y_z_tmpi ~ y_z_tmpj+age_z+EstimatedTICV_z, 
                     random = ~ 1 | sub_id, data=dat_total_subset, na.action=na.omit)
    beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpj", "lower"]
    beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpj", "est."]
    beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpj", "upper"]
    print(intervals(fit_z_tmp, which=c("fixed")))
    
    # append
    df_result_ges[nrow(df_result_ges) + 1,] = c(i, unname(tmp_p), b_tmp, tmp_interc,
                                                unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
  }
  
  
  write.csv(df_result_ges, 
            file = paste0("df_result_csf_", j, ".csv"), 
            row.names = FALSE)
}


################## 2) csf and t1 ##################
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/csfmap.csv")
# exclude 10004_scan1
dat_total_subset <- dat_total[!(dat_total$fsid=='10004_scan1'),]
dat_total_subset$sub_id <- as.factor(dat_total_subset$sub_id)
dat_total_subset$pial_sum <- rowSums(subset(dat_total_subset, 
                                              select = c(lh_pial, rh_pial)), 
                                       na.rm = TRUE)
dat_total_subset$thickness_mean <- rowMeans(subset(dat_total_subset, 
                                                   select = c(lh_MeanThickness_thickness, 
                                                              rh_MeanThickness_thickness)), 
                                            na.rm = TRUE)

v1_interest <- c("csf_60", "csf_65", "csf_70", "csf_75",
                    "csf_80", "csf_85", "csf_90", "csf_95")
v2_interest <- c("BrainSegVolNotVent_x", "TotalGrayVol", "SubCortGrayVol",
                 "CortexVol", "pial_sum", "thickness_mean")

dat_total_subset$age_z <- scale(dat_total_subset$age , center = TRUE, scale = TRUE)
dat_total_subset$eTIV_z <- scale(dat_total_subset$eTIV_x , center = TRUE, scale = TRUE)
df_result_ges <- data.frame(ROI1=character(),
                            ROI2=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(), stringsAsFactors=FALSE)
for (i in v2_interest) {
  for (j in v1_interest) {
    
    f <- formula(paste(i, "~ age + eTIV_x + ", j))
    fit_tmp <- lme(f, random = ~ 1 | sub_id, 
                   data=dat_total_subset, na.action=na.omit)
    print(summary(fit_tmp))
    sum_tmp <- summary(fit_tmp)
    tmp_p <- sum_tmp$tTable[j, "p-value"]
    b_tmp <- unname(sum_tmp$coefficients$fixed[j])
    # standardized beta and intervals
    dat_total_subset$y_z_tmpi <- scale(dat_total_subset[i], center = TRUE, scale = TRUE)
    dat_total_subset$y_z_tmpj <- scale(dat_total_subset[j], center = TRUE, scale = TRUE)
    fit_z_tmp <- lme(y_z_tmpi ~ y_z_tmpj + age_z + eTIV_z, 
                     random = ~ 1 | sub_id, 
                     data=dat_total_subset, 
                     na.action=na.omit)
    sum_tmp_z <- summary(fit_z_tmp)
    beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpj", "lower"]
    beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpj", "est."]
    beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpj", "upper"]
    # append
    df_result_ges[nrow(df_result_ges) + 1,] = c(i, j, unname(tmp_p), b_tmp, 
                                                unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
  }
}

write.csv(df_result_ges, "csf_t1_ControlTICV_AGE.csv", 
          row.names = FALSE)


################## 3) csf and diffusion ##################
file_list <- c('/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/matlab_icvf_all.csv')

for (metric in file_list) {
  print(metric)
  
  # get the metric name
  basename <- basename(metric)
  split_string <- strsplit(basename, split = ".csv")
  split_elements <- split_string[[1]]
  split_string2 <- strsplit(split_elements, split = "_")
  metric_name <- split_string2[[1]][2]
  
  # import data
  importDat <- read.csv(metric, header = TRUE)
  importDat$sub <- as.factor(importDat$sub)
  importDat_sub <- importDat[!((importDat$fsid=="10003_scan3") 
                                   | (importDat$fsid=="10004_scan1")
                                   | (importDat$fsid=="10012_scan2")
                                   | (importDat$fsid=="10004_scan3")
                                   | (importDat$fsid=="10011_scan2")),]
  dat_total_subset <- importDat_sub
  
  # convert to the long format
  dat_total_subset$scan <- unlist(lapply(dat_total_subset$fsid, function(x) strsplit(x, '_')[[1]][2]))
  dat_total_subset$csf_80_z <- scale(dat_total_subset$csf_80, center = TRUE, scale = TRUE)
  dat_total_subset$age_z <- scale(dat_total_subset$age, center = TRUE, scale = TRUE)
  dat_total_subset$motion_rel_z <- scale(dat_total_subset$motion_rel, center = TRUE, scale = TRUE)
  list_subregions <- list()
  list_subregions_long <- list()
  for (trac in names(dat_total_subset)) {
    if (grepl('_left', trac)) {
      print(trac)
      # find the left tract index
      idx_left <- which(trac == names(dat_total_subset))
      
      # find the right tract index
      trac_right <- sub("_left", "_right", trac)
      idx_right <- which(trac_right == names(dat_total_subset))
      sub_col_name <- c(trac, trac_right, "sub", "csf_80", 
                        "age", "scan", 
                        "csf_80_z", "age_z", "motion_rel", "motion_rel_z")
      aa <- subset(dat_total_subset, select=sub_col_name)
      aa_long <- aa %>% pivot_longer(
        cols = c(1, 2),
        names_to = c('.value', 'hemi'),
        names_pattern = "(.*)_(.*)") %>%
        convert_as_factor(sub, hemi, scan)
      list_subregions[[length(list_subregions)+1]] <- aa
      list_subregions_long[[length(list_subregions_long)+1]] <- na.omit(aa_long)
    }
  }
  
  # run the models without the interaction term
  df_result_ges <- data.frame(ROI=character(),
                              p=double(),
                              b=double(),
                              beta=double(),
                              beta_lower=double(),
                              beta_upper=double(), stringsAsFactors=FALSE)
  for (i in 1:length(list_subregions_long)) {
    print(i)
    data <- list_subregions_long[[i]] 
    data_dv <- names(data)[length(data)]
    print(data_dv)
    fit_func <- formula(paste(data_dv, "~ hemi + age + motion_rel + csf_80"))
    fit_tmp <- lme(fit_func, random = ~ 1 | sub, 
                   data=data, na.action=na.omit)
    print(summary(fit_tmp))
    sum_tmp <- summary(fit_tmp)
    tmp_p <- sum_tmp$tTable["csf_80", "p-value"]
    b_tmp <- unname(sum_tmp$coefficients$fixed['csf_80'])
    # standardized beta and intervals
    idx_tmp <- length(data)
    y_z_tmp <- scale(data[,idx_tmp], center = TRUE, scale = TRUE)
    fit_z_tmp <- lme(y_z_tmp ~ hemi + age_z + csf_80_z + motion_rel_z, 
                     random = ~ 1 | sub, 
                     data=data, na.action=na.omit)
    beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["csf_80_z", "lower"]
    beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["csf_80_z", "est."]
    beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["csf_80_z", "upper"]
    print(intervals(fit_z_tmp, which=c("fixed")))
    
    # append
    df_result_ges[nrow(df_result_ges) + 1,] = c(data_dv, unname(tmp_p), b_tmp, 
                                                unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
  }
  
  # the above part is not sig, so did not edit the below part yet.
  # group CC into CC_RG(CC_Rostrum(CC_1), CC_Genu(CC_2)), CC_body(CC_345), CC_IS(CC_Isthmus(CC_6), CC_Splenium(CC_7))
  dat_total_subset$CC_RG <- rowMeans(dat_total_subset[, c('CC_1', 'CC_2')])
  dat_total_subset$CC_body <- rowMeans(dat_total_subset[, c('CC_3', 'CC_4', 'CC_5')])
  dat_total_subset$CC_IS <- rowMeans(dat_total_subset[, c('CC_6', 'CC_7')])
  
  for (i in tail(names(dat_total_subset),3)) {
    fit_func <- formula(paste(i, "~ age + csf_80 + motion_rel"))
    fit_tmp <- lme(fit_func, 
                   random = ~ 1 | sub, 
                   data=dat_total_subset, 
                   na.action=na.omit)
    print(summary(fit_tmp))
    sum_tmp <- summary(fit_tmp)
    tmp_p <- sum_tmp$tTable["csf_80", "p-value"]
    b_tmp <- unname(sum_tmp$coefficients$fixed['csf_80'])
    # standardized beta and intervals
    idx_tmp <- which(i == names(dat_total_subset))
    y_z_tmp <- scale(dat_total_subset[,idx_tmp], center = TRUE, scale = TRUE)
    fit_z_tmp <- lme(y_z_tmp ~ age_z + csf_80_z + motion_rel_z, 
                     random = ~ 1 | sub, 
                     data=dat_total_subset, na.action=na.omit)
    beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["csf_80_z", "lower"]
    beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["csf_80_z", "est."]
    beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["csf_80_z", "upper"]
    print(intervals(fit_z_tmp, which=c("fixed")))
    
    # append
    df_result_ges[nrow(df_result_ges) + 1,] = c(i, unname(tmp_p), b_tmp, 
                                                unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
  }
  
  # rank p value
  write.csv(df_result_ges, 
            file = paste0("df_csfmap80_", metric_name, ".csv"), 
            row.names = FALSE)
}

