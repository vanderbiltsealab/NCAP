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


########################## 1) t1 ~ hormone ##########################
# saliva sample needs to remove 10001_scan1
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/t1_all.csv")
dat_total_subset_25 <- dat_total[!(dat_total$scan_id=="10001_scan1" 
                                   | dat_total$scan_id=="10004_scan1"),]
dat_total_subset_25$sub_id <- as.factor(dat_total_subset_25$sub_id)

dat_total_subset_25$age_z <- scale(dat_total_subset_25$age, center = T, scale = T)
dat_total_subset_25$num_pri_preg_z <- scale(dat_total_subset_25$num_pri_preg, center = T, scale = T)
dat_total_subset_25$progesterone_z <- scale(dat_total_subset_25$progesterone, center = T, scale = T)
dat_total_subset_25$estradiol_z <- scale(dat_total_subset_25$estradiol, center = T, scale = T)
dat_total_subset_25$testosterone_z <- scale(dat_total_subset_25$testosterone, center = T, scale = T)
dat_total_subset_25$cortisol_z <- scale(dat_total_subset_25$cortisol, center = T, scale = T)
dat_total_subset_25$eTIV_z <- scale(dat_total_subset_25$eTIV_x, center = T, scale = T)
dat_total_subset_25$pial_sum <- rowSums(subset(dat_total_subset_25, 
                                            select = c(lh_pial, rh_pial)), 
                                     na.rm = TRUE)
dat_total_subset_25$thickness_mean <- rowMeans(subset(dat_total_subset_25, 
                                                   select = c(lh_MeanThickness_thickness, 
                                                              rh_MeanThickness_thickness)), 
                                            na.rm = TRUE)

# y variables of interest
y_interest <- c('BrainSegVolNotVent_x', 'TotalGrayVol', 'CerebralWhiteMatterVol',
                'SubCortGrayVol', 'CortexVol', 'pial_sum', 'thickness_mean')
hormone_interest <- c('progesterone', 'estradiol', 'testosterone', 'cortisol')

# covariate: eTIV_x and age
df_result_ges <- data.frame(ROI1=character(),
                            ROI2=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(), stringsAsFactors=FALSE)
for (i in hormone_interest) {
  for (j in y_interest) {
    f <- formula(paste(j, " ~ eTIV_x + age +", i))
    
    fit_tmp <- lme(f, random = ~ 1 | sub_id, 
                   data=dat_total_subset_25, na.action=na.omit)
    print(summary(fit_tmp))
    sum_tmp <- summary(fit_tmp)
    tmp_p <- sum_tmp$tTable[i, "p-value"]
    b_tmp <- unname(sum_tmp$coefficients$fixed[i])
    # standardized beta and intervals
    dat_total_subset_25$y_z_tmpi <- scale(dat_total_subset_25[i], center = TRUE, scale = TRUE)
    dat_total_subset_25$y_z_tmpj <- scale(dat_total_subset_25[j], center = TRUE, scale = TRUE)
    fit_z_tmp <- lme(y_z_tmpj ~ eTIV_z + age_z + y_z_tmpi, 
                     random = ~ 1 | sub_id, 
                     data=dat_total_subset_25, na.action=na.omit)
    print(summary(fit_z_tmp))
    sum_tmp_z <- summary(fit_z_tmp)
    beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpi", "lower"]
    beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpi", "est."]
    beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpi", "upper"]
    
    # append
    df_result_ges[nrow(df_result_ges) + 1,] = c(i, j, unname(tmp_p), b_tmp, 
                                                unname(beta_tmp), unname(beta_lower_tmp), unname(beta_upper_tmp))

  }
}
write.csv(df_result_ges, "df_result_ges_t1_and_hormone.csv", 
          row.names=F)


########################## 2) t1 ~ hair hormone ##########################
# hair sample does not need to remove 10001_scan1
# hair sample does not converge if including segment as a random slope
dat_total <- read.csv("/Users/yanbinniu/Projects/NCAP/scripts/github/NCAP/data/t1_all.csv")
dat_total_subset_25 <- dat_total[!(dat_total$scan_id=="10004_scan1"),]
dat_total_subset_25$sub_id <- as.factor(dat_total_subset_25$sub_id)

dat_total_subset_25$age_z <- scale(dat_total_subset_25$age, center = T, scale = T)
dat_total_subset_25$cortisol_win_centered_z <- scale(dat_total_subset_25$cortisol_win_centered, center = T, scale = T)
dat_total_subset_25$cortisone_win_centered_z <- scale(dat_total_subset_25$cortisone_win_centered, center = T, scale = T)
dat_total_subset_25$progesterone_win_centered_z <- scale(dat_total_subset_25$progesterone_win_centered, center = T, scale = T)
dat_total_subset_25$segment_sample_use_z <- scale(dat_total_subset_25$segment_sample_use, center = T, scale = T)
dat_total_subset_25$distance_hair_sample <- dat_total_subset_25$segment_sample_use - 1
dat_total_subset_25$distance_hair_sample_z <- scale(dat_total_subset_25$distance_hair_sample, center = T, scale = T)
dat_total_subset_25$eTIV_z <- scale(dat_total_subset_25$eTIV_x, center = T, scale = T)
dat_total_subset_25$pial_sum <- rowSums(subset(dat_total_subset_25,
                                               select = c(lh_pial, rh_pial)),
                                        na.rm = TRUE)
dat_total_subset_25$thickness_mean <- rowMeans(subset(dat_total_subset_25,
                                                      select = c(lh_MeanThickness_thickness,
                                                                 rh_MeanThickness_thickness)),
                                               na.rm = TRUE)

# y variables of interest
y_interest <- c('BrainSegVolNotVent_x', 'TotalGrayVol', 'CerebralWhiteMatterVol',
                'SubCortGrayVol', 'CortexVol', 'pial_sum', 'thickness_mean')
hormone_interest <- c('cortisol_win_centered', 'cortisone_win_centered')

# covariate: eTIV_x and age
df_result_ges <- data.frame(ROI1=character(),
                            ROI2=character(),
                            p=double(),
                            b=double(),
                            beta=double(),
                            beta_lower=double(),
                            beta_upper=double(), stringsAsFactors=FALSE)
for (i in hormone_interest) {
  for (j in y_interest) {
    print(i)
    print(j)
    
    # standardized beta
    dat_total_subset_25$y_z_tmpi <- scale(dat_total_subset_25[i], center = TRUE, scale = TRUE)
    dat_total_subset_25$y_z_tmpj <- scale(dat_total_subset_25[j], center = TRUE, scale = TRUE)
    
    f <- formula(paste(j, " ~ distance_hair_sample + eTIV_x + age + ", i))
    fit_tmp <- lme(f, random = ~ 1 | sub_id,
                   data=dat_total_subset_25, 
                   na.action=na.omit)
    print(summary(fit_tmp))
    sum_tmp <- summary(fit_tmp)
    tmp_p <- sum_tmp$tTable[i, "p-value"]
    b_tmp <- unname(sum_tmp$coefficients$fixed[i])
    
    fit_z_tmp <- lme(y_z_tmpj ~ distance_hair_sample_z + eTIV_z + age_z + y_z_tmpi,
                     random = ~ 1 | sub_id,
                     data=dat_total_subset_25, 
                     na.action=na.omit)
    print(summary(fit_z_tmp))
    
    sum_tmp_z <- summary(fit_z_tmp)
    tmp_p <- sum_tmp_z$tTable["y_z_tmpi", "p-value"]
    beta_lower_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpi", "lower"]
    beta_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpi", "est."]
    beta_upper_tmp <- intervals(fit_z_tmp, which=c("fixed"))['fixed']$fixed["y_z_tmpi", "upper"]

    # append
    df_result_ges[nrow(df_result_ges) + 1,] = c(i, j, unname(tmp_p), b_tmp, 
                                                unname(beta_tmp), 
                                                unname(beta_lower_tmp), unname(beta_upper_tmp))

  }
}

write.csv(df_result_ges, "df_result_ges_t1_and_hormone_hair.csv",
          row.names=F)

