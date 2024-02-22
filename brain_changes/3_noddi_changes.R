library(nlme)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(MuMIn)
library(r2glmm)
library(sjstats)
library(rstatix)
library(reghelper)
library(tidyr)

# switch to the desired workng directory
setwd('/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data')

file_list <- c("matlab_icvf_all.csv", "matlab_od_all.csv",
               "fa_all.csv", "md_all.csv", "ad_all.csv", "rd_all.csv")

################## matlab noddi and gestational age ##################
# df_final as the final merged dataframe for heatmap plotting
ROI_interest <- c("AF", "ATR", "CG", "CST", "FPT", "ICP", "IFO", "ILF", "OR", "POPT", "SCP", "SLF_I",
                  "SLF_II", "SLF_III", "STR", "UF", "T_PREM", "T_PAR", "T_OCC", "ST_FO", "ST_PREM", "CC_RG",
                  "CC_body", "CC_IS")
df_final <- data.frame(ROI=ROI_interest, stringsAsFactors=FALSE)

for (metric in file_list) {
  print(metric)
  
  importDat <- read.csv(metric, header = TRUE)
  importDat$sub <- as.factor(importDat$sub)
  importDat_sub <- subset(importDat, ges_time>0)
  # remove nan
  importDat_sub <- importDat_sub[!((importDat_sub$fsid=="10003_scan3") 
                                   | (importDat_sub$fsid=="10012_scan2")
                                   | (importDat_sub$fsid=="10004_scan3")
                                   | (importDat_sub$fsid=="10011_scan2")),]
  dat_total_subset <- importDat_sub
  
  # convert to the long format
  dat_total_subset$scan <- unlist(lapply(dat_total_subset$fsid, function(x) strsplit(x, '_')[[1]][2]))
  dat_total_subset$ges_time_z <- scale(dat_total_subset$ges_time, center = TRUE, scale = TRUE)
  dat_total_subset$age_z <- scale(dat_total_subset$age, center = TRUE, scale = TRUE)
  dat_total_subset$motion_rel_z <- scale(dat_total_subset$motion_rel, center = TRUE, scale = TRUE)
  dat_total_subset$progesterone_z <- scale(dat_total_subset$progesterone, center = TRUE, scale = TRUE)
  dat_total_subset$estradiol_z <- scale(dat_total_subset$estradiol, center = TRUE, scale = TRUE)
  dat_total_subset$testosterone_z <- scale(dat_total_subset$testosterone, center = TRUE, scale = TRUE)
  dat_total_subset$cortisol_z <- scale(dat_total_subset$cortisol, center = TRUE, scale = TRUE)
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
      sub_col_name <- c(trac, trac_right, "cortisol", "testosterone", "estradiol", "progesterone",
                        "ges_time", "age", "sub", "motion_rel", "progesterone_z", "estradiol_z",
                        "testosterone_z", "cortisol_z", "ges_time_z", "scan", "age_z", "motion_rel_z")
      
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
  
  #  ~ ges_time
  plot_list <- list()
  df_result_ges <- data.frame(ROI=character(),
                              p=double(),
                              b=double(),
                              beta=double(),
                              beta_lower=double(),
                              beta_upper=double(),stringsAsFactors=FALSE)
  for (i in 1:length(list_subregions_long)) {
    data <- list_subregions_long[[i]] 
    data_dv <- names(data)[length(data)]
    print(data_dv)
    fit_func <- formula(paste(data_dv, "~ hemi + age + ges_time + motion_rel"))
    fit_tmp <- lme(fit_func, random = ~ 1 | sub, data=data, na.action=na.omit)
    print(summary(fit_tmp))
    sum_tmp <- summary(fit_tmp)
    b_tmp <- unname(sum_tmp$coefficients$fixed['ges_time'])
    # standardized beta and intervals
    print("standardized beta and intervals")
    idx_tmp <- length(data)
    data$y_z_tmp <- scale(data[,idx_tmp], center = TRUE, scale = TRUE)
    fit_z_tmp <- lme(y_z_tmp ~ hemi + age_z + ges_time_z + motion_rel_z, 
                     random = ~ 1 | sub, data=data, na.action=na.omit)
    sum_tmp_z <- summary(fit_z_tmp)
    tmp_p <- sum_tmp_z$tTable["ges_time_z", "p-value"]
    beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
    beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
    beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
    print(intervals(fit_z_tmp, which=c("fixed")))
    
    # append
    df_result_ges[nrow(df_result_ges) + 1,] = c(data_dv, unname(tmp_p), b_tmp, 
                                                unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
  }
  
  # group CC into CC_RG(CC_Rostrum(CC_1), CC_Genu(CC_2)), CC_body(CC_345), CC_IS(CC_Isthmus(CC_6), CC_Splenium(CC_7))
  dat_total_subset$CC_RG <- rowMeans(dat_total_subset[, c('CC_1', 'CC_2')])
  dat_total_subset$CC_body <- rowMeans(dat_total_subset[, c('CC_3', 'CC_4', 'CC_5')])
  dat_total_subset$CC_IS <- rowMeans(dat_total_subset[, c('CC_6', 'CC_7')])
  
  for (i in tail(names(dat_total_subset),3)) {
    fit_func <- formula(paste(i, "~ age + ges_time + motion_rel"))
    fit_tmp <- lme(fit_func, random = ~ 1 | sub, 
                   data=dat_total_subset, na.action=na.omit)
    print(summary(fit_tmp))
    sum_tmp <- summary(fit_tmp)
    b_tmp <- unname(sum_tmp$coefficients$fixed['ges_time'])
    # standardized beta and intervals
    idx_tmp <- which(i == names(dat_total_subset))
    y_z_tmp <- scale(dat_total_subset[,idx_tmp], center = TRUE, scale = TRUE)
    fit_z_tmp <- lme(y_z_tmp ~ age_z + ges_time_z + motion_rel_z, 
                     random = ~ 1 | sub, data=dat_total_subset, na.action=na.omit)
    sum_tmp_z <- summary(fit_z_tmp)
    tmp_p <- sum_tmp_z$tTable["ges_time_z", "p-value"]
    beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "lower"]
    beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "est."]
    beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["ges_time_z", "upper"]
    print(intervals(fit_z_tmp, which=c("fixed")))
    
    # append
    df_result_ges[nrow(df_result_ges) + 1,] = c(i, unname(tmp_p), b_tmp, 
                                                unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))
  }
  
  current_metric <- strsplit(metric, "_")[[1]][length(strsplit(metric, "_")[[1]])-1]
  new_column_name_beta <- paste0(current_metric, "_beta")
  df_final[[new_column_name_beta]] <- df_result_ges$beta
  new_column_name_p <- paste0(current_metric, "_p")
  df_final[[new_column_name_p]] <- df_result_ges$p
  
  write.csv(df_result_ges, 
            file = paste0("df_result_ges_noddi_matlab_", current_metric, ".csv"), 
            row.names = FALSE)
}

write.csv(df_final, "all_dti_beta_ctrl_motion.csv", row.names=FALSE)


