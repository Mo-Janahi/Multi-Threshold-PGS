
library(TTR)
library(smoother)
library(ADNIMERGE)
library(laGP)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(glmnet)
library(qqman)
library(pROC)
library(plotly)
library(readxl)
library(scales)
library(matrixStats)
library(neuroCombat)
library(psychometric)
library(patchwork)
library(stringr)
library(dplyr)
library(lubridate)
library(lme4)
library(lmerTest)
library(tidyr)
library(viridis)
library(cowplot)
library(ggforce)
library(reshape2)
library(shadowtext)
library(poolr)

source("~/OneDrive - University College London/ANDRE STAGE/Project 1 PGS Nomogram/Scripts/Basic Functions.R")
source("~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/Plotting Functions.R")
source("~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/Preprocessing Functions.R")
source("~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/Extra Functions.R")

plt_folder = "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Paper/Figures/"


ukb_table <- read.csv("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_features_20200830.csv")
adni_table <- merge( ADNIMERGE::adnimerge , ADNIMERGE::ucsffsx51 , by = c("IMAGEUID") )
ucl_table <- read.csv("~/OneDrive - University College London/ANDRE STAGE/Datasets/ucl46_data/UCL_46_HV_subset.csv")
epad_table = read.csv("OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/COMBINE_EPAD.csv")

## DATASET PRE-PROCESSING
x = PRERPOCESS_ALL( ukb_table , adni_table , ucl_table , epad_table , GIF = FALSE)

ukb_all <- x$ukb_all
adni_all<- x$adni_all
ucl_all <- x$ucl_all
epad_all <- x$epad_all

x = HARMONIZE_ALL(ukb_all , adni_all , ucl_all , epad_all)

cols = grep("_nc", names(x$ukb_all_nc) , value = TRUE) 
#c(prs_columns_nc,prs_columns_nc_gpc,"clean_hv_bilateral_nc" ,
#         "clean_hv_bilateral_nc_gpc","clean_hv_bilateral_nc_gpc_icv",
#         "ICV_nc","ICV_nc_gpc" )

ukb_all = merge( ukb_all , x$ukb_all_nc[ , c("eid",cols) ] , by="eid" , all.x = TRUE)
adni_all = merge( adni_all , x$adni_all_nc[ , c("PTID","VISCODE2",cols) ] , by=c("PTID","VISCODE2") , all.x = TRUE)
ucl_all = merge( ucl_all , x$ucl_all_nc[ , c("nshdid",cols) ] , by="nshdid" , all.x = TRUE)
epad_all = merge( epad_all , x$epad_all_nc[ , c("GWAS_ID","visit",cols) ] , by=c("GWAS_ID","visit") , all.x = TRUE)

nc_cols = c("clean_hv_bilateral" , "clean_hv_bilateral_nc" , "ICV" , "ICV_nc" , unlist(lapply(  grep("^INTERSECT_PRS_TH.*_nc$" , names(ukb_all) , value = TRUE)  , function(x)  c( str_replace(x, "_nc" , "" ) , x) )))
harmonized_table = rbind( 
  cbind( ukb_all[  ,nc_cols] ,table = "ukb") , cbind( adni_all[ ,nc_cols],table = "adni"),
  cbind( epad_all[ ,nc_cols],table = "epad"),  cbind( ucl_all[  ,nc_cols],table = "ucl")
)

myplots = lapply( nc_cols ,  function(col ) ggplot(harmonized_table, aes(x = get(col), fill = table)) +
           geom_density( aes(y = after_stat(density)) ,position = "identity", alpha = 0.5)  +
           xlab(col) + scale_fill_brewer(palette = "Set1") )


pdf(file = paste0(plt_folder,"Supplementary Figure - Harmonization Results - regressed_alone",".pdf") , width = 30 , height = 15)

plot_grid(plotlist=myplots)

dev.off()


write.csv(ukb_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/with_imaging_40k/PREPROCESSED_UKB_27_01_2024 Regressed Alone.csv")
write.csv(adni_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data//PREPROCESSED_ADNI_27_01_2024 Regressed Alone.csv")
write.csv(ucl_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/ucl46_data//PREPROCESSED_I46_27_01_2024 Regressed Alone.csv")
write.csv(epad_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/PREPROCESSED_EPAD_27_01_2024 Regressed Alone.csv")


ukb_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/with_imaging_40k/PREPROCESSED_UKB_27_01_2024.csv")
adni_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data//PREPROCESSED_ADNI_27_01_2024.csv")
ucl_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/ucl46_data//PREPROCESSED_I46_27_01_2024.csv")
epad_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/PREPROCESSED_EPAD_27_01_2024.csv")

prs_columns <- grep("^INTERSECT_PRS_TH.*[^_nc]$" , names(ucl_all) , value = TRUE)
prs_columns_nc <- paste0(prs_columns,"_nc")
prs_columns_nc_gpc <- paste0(prs_columns,"_nc_gpc")
prs_columns_nc_gpc_icv <- paste0(prs_columns,"_nc_gpc_icv")


ukb_all$id = ukb_all$eid
adni_all$id = adni_all$PTID
ucl_all$id = ucl_all$nshdid
epad_all$id = epad_all$GWAS_ID

adni_all = adni_all[ adni_all$age < max(ukb_all$age) , ]
ucl_all = ucl_all[ ucl_all$age < max(ukb_all$age) , ]
epad_all = epad_all[ epad_all$age < max(ukb_all$age) , ]



#### ABSTRACTION EXPERIMENTS ####
SEEDS = 4560:4590

ABSTRACTION_RESULTS = data.frame(matrix(ncol = 25, nrow = 0))
ABSTRACTION_RESULTS = setNames(ABSTRACTION_RESULTS , c("run_ind", "exp_ind" , "experiment", "dataset" , "seed" , "RMSE" , "ASSOC", "RMSE_t" , "ASSOC_t","RMSE_t2" , "ASSOC_t2" ,
                                                       "RMSE_ADNI" , "ASSOC_ADNI", "RMSE_UCL" , "ASSOC_UCL", "RMSE_EPAD" , "ASSOC_EPAD" ,
                                                       "RMSE_ADNI_CN" , "ASSOC_ADNI_CN","RMSE_ADNI_MCI" , "ASSOC_ADNI_MCI","RMSE_ADNI_AD" , "ASSOC_ADNI_AD",
                                                       "exp_grp"))

DETAILED_RESULTS = list()

exp_ind = 1
run_ind = 1
mod_num = 1
pred_col = "clean_hv_bilateral_nc_gpc_icv"

COVS = append( append( list( c("age"),
                             c("age","sex"),
                             c("age","sex","ICV_nc_gpc") ) ,
                       lapply(prs_columns_nc_gpc, function(x) c("age", "sex","ICV_nc_gpc",x)) ) ,
               list(c("age","sex","ICV_nc_gpc","PRS_LASSO")))

for( covs in COVS){
  for ( seed in SEEDS ) {
    #  lasso_run = ("PRS_LASSO" %in% covs)
    #  covs_tmp = covs
    #  if (lasso_run)
    #    covs_tmp = c( covs[-length(covs)] , prs_columns_nc_gpc)
    # 
    #  grab_cols = c("id" , covs_tmp , pred_col)
    # 
    #  ukb_rows = complete.cases(ukb_all[,grab_cols])
    #  adni_rows = complete.cases(adni_all[,grab_cols])
    #  ucl_rows = complete.cases(ucl_all[,grab_cols])
    #  epad_rows = complete.cases(epad_all[,grab_cols])
    # 
    #  if (lasso_run){
    # 
    #    lasso = LASSO_SCORE(ukb_all[ukb_rows , ] , ukb_all[ukb_rows,pred_col] , prs_columns_nc_gpc , print_weights = TRUE , seed = seed)
    #    prs_lasso_col = "PRS_LASSO"
    # 
    #    ukb_all[ukb_rows,prs_lasso_col] = lasso$score
    #    adni_all[adni_rows, prs_lasso_col] = rowSums(sweep( adni_all[adni_rows , prs_columns_nc_gpc] ,2 , lasso$weights , "*"))
    #    ucl_all[ucl_rows,prs_lasso_col] = rowSums(sweep( ucl_all[ucl_rows , prs_columns_nc_gpc] ,2 , lasso$weights , "*"))
    #    epad_all[epad_rows, prs_lasso_col] = rowSums(sweep( epad_all[epad_rows , prs_columns_nc_gpc] ,2 , lasso$weights , "*"))
    # 
    #    grab_cols = c(grab_cols , prs_lasso_col)
    # 
    #  }
    # 
    # set.seed(seed)
    # n = sum(ukb_rows)
    # train_ind = sample( 1:n , .5*n)
    # train_table = ukb_all[ukb_rows , ][train_ind, ]
    # test_table = ukb_all[ukb_rows , ][-train_ind, ]
    # 
    # epad_all$DX = epad_all$apoe_result
    # ucl_all$DX = ucl_all$PHENOTYPE
    # adni_temp = adni_all[adni_rows,c(grab_cols,"DX")]
    # ucl_temp = ucl_all[ucl_rows,c(grab_cols,"DX")]
    # epad_temp = epad_all[epad_rows,c(grab_cols,"DX")]
    # 
    # out_sample_test_table = rbind(adni_temp , ucl_temp , epad_temp )

     covs_lin = paste(covs,collapse = ' + ')
     print(paste( pred_col," ~ ",covs_lin ,"WITH:",seed))

    # table = train_table
    # lin_m = lm( as.formula( paste0(pred_col," ~ " , covs_lin) ) , data = table )
    # sqrt(mean((lin_m$fitted.values - table[ , pred_col])^2))


 # GPR MODELS ##
 # gpr = GPR_PREDICT(train_set = train_table , test_set = test_table , test_set_2 = out_sample_test_table , x_cols=covs , y_col = pred_col , out_sample_sizes = c( nrow(adni_temp) , nrow(ucl_temp), nrow(epad_temp)  ) )

 hv_pred_mean_col = paste(pred_col,mod_num,seed,"predicted_mean",sep = "_")
 hv_pred_sd_col = paste(pred_col,mod_num,seed,"predicted_sd",sep = "_")
 z_score_col = paste(pred_col,mod_num,seed,sep = "_")

# ukb_all[ukb_rows , paste0(z_score_col,"_train") ] = NA
# ukb_all[ukb_rows , paste0(z_score_col,"_test") ] = NA
# 
# ukb_all[ukb_rows , paste0(hv_pred_mean_col,"_train") ] = NA
# ukb_all[ukb_rows , paste0(hv_pred_mean_col,"_test") ] = NA
# 
# ukb_all[ukb_rows , paste0(hv_pred_sd_col,"_train") ] = NA
# ukb_all[ukb_rows , paste0(hv_pred_sd_col,"_test") ] = NA
# 
# ukb_all[ukb_rows,][ train_ind , paste0(hv_pred_mean_col,"_train") ] =  gpr$training_set$predicted_means
# ukb_all[ukb_rows,][ -train_ind , paste0(hv_pred_mean_col,"_test") ] = gpr$test_set$predicted_means
# 
# ukb_all[ukb_rows,][ train_ind , paste0(hv_pred_sd_col,"_train") ] = gpr$training_set$predicted_sds
# ukb_all[ukb_rows,][ -train_ind , paste0(hv_pred_sd_col,"_test") ] = gpr$test_set$predicted_sds
# 
# ukb_all[ukb_rows,][ train_ind , paste0(z_score_col,"_train") ] = gpr$training_set$zscores
# ukb_all[ukb_rows,][ -train_ind , paste0(z_score_col,"_test") ] = gpr$test_set$zscores
# 
# adni_sub = gpr$test_2_set[  1:sum(adni_rows) , ]
# ucl_sub  = gpr$test_2_set[ (1:sum(ucl_rows) + sum(adni_rows) ) , ]
# epad_sub = gpr$test_2_set[ (1:sum(epad_rows) + sum(adni_rows) + sum(ucl_rows) ) , ]
# 
# adni_all[adni_rows,hv_pred_mean_col ] = adni_sub$predicted_means
# ucl_all[ ucl_rows, hv_pred_mean_col ] = ucl_sub$predicted_means
# epad_all[epad_rows,hv_pred_mean_col ] = epad_sub$predicted_means
# 
# adni_all[adni_rows,hv_pred_sd_col ] = adni_sub$predicted_sds
# ucl_all[ ucl_rows, hv_pred_sd_col ] = ucl_sub$predicted_sds
# epad_all[epad_rows,hv_pred_sd_col ] = epad_sub$predicted_sds
# 
# adni_all[adni_rows,z_score_col ] = adni_sub$zscores
# ucl_all[ ucl_rows, z_score_col ] = ucl_sub$zscores
# epad_all[epad_rows,z_score_col ] = epad_sub$zscores


rmse = sqrt(mean( (ukb_all[, pred_col] - ukb_all[ , paste0(hv_pred_mean_col,"_train")] )^2 , na.rm = TRUE))
corr = cor(ukb_all[, pred_col] , ukb_all[ , paste0(hv_pred_mean_col,"_train")] , use = "complete.obs")

rmse_t = sqrt(mean((ukb_all[, pred_col] - ukb_all[ , paste0(hv_pred_mean_col,"_test")] )^2, na.rm = TRUE))
corr_t = cor(ukb_all[, pred_col] , ukb_all[ , paste0(hv_pred_mean_col,"_test")] , use = "complete.obs")


rmse_t_2_adni = sqrt(mean( (adni_all[, pred_col] - adni_all[ , hv_pred_mean_col] )^2 , na.rm = TRUE))
corr_t_2_adni = cor(adni_all[, pred_col] , adni_all[ , hv_pred_mean_col] , use = "complete.obs")
rmse_t_2_ucl = sqrt(mean( (ucl_all[, pred_col] - ucl_all[ , hv_pred_mean_col] )^2 , na.rm = TRUE))
corr_t_2_ucl = cor(ucl_all[, pred_col] ,  ucl_all[ , hv_pred_mean_col] , use = "complete.obs")
rmse_t_2_epad = sqrt( mean((epad_all[, pred_col] - epad_all[ , hv_pred_mean_col] )^2, na.rm = TRUE))
corr_t_2_epad = cor(epad_all[, pred_col] , epad_all[ , hv_pred_mean_col] , use = "complete.obs")

rmse_t_2_adni_cn = sqrt(mean( (adni_all[adni_all$DX %in% c("CN"), pred_col] - adni_all[adni_all$DX %in% c("CN") , hv_pred_mean_col] )^2 , na.rm = TRUE))
corr_t_2_adni_cn = cor(adni_all[ adni_all$DX %in% c("CN") , pred_col] , adni_all[ adni_all$DX %in% c("CN") , hv_pred_mean_col] , use = "complete.obs")

rmse_t_2_adni_mci = sqrt(mean( (adni_all[adni_all$DX %in% c("MCI"), pred_col] - adni_all[adni_all$DX %in% c("MCI") , hv_pred_mean_col] )^2 , na.rm = TRUE))
corr_t_2_adni_mci = cor(adni_all[adni_all$DX %in% c("MCI"), pred_col] , adni_all[adni_all$DX %in% c("MCI") , hv_pred_mean_col] , use = "complete.obs")

rmse_t_2_adni_ad = sqrt(mean( (adni_all[ adni_all$DX %in% c("Dementia") , pred_col] - adni_all[adni_all$DX %in% c("Dementia") , hv_pred_mean_col] )^2 , na.rm = TRUE))
corr_t_2_adni_ad = cor(adni_all[adni_all$DX %in% c("Dementia"), pred_col] , adni_all[adni_all$DX %in% c("Dementia") , hv_pred_mean_col] , use = "complete.obs")

rmse_t_2 = mean( c(rmse_t_2_adni , rmse_t_2_epad) )
corr_t_2 = mean( c(corr_t_2_adni , corr_t_2_epad) )

    ABSTRACTION_RESULTS <- rbind(ABSTRACTION_RESULTS , 
                                 list(run_ind , exp_ind , paste0("gpr_",paste(covs,collapse = '_')) , "ukbb", seed ,
                                      rmse , corr, rmse_t , corr_t, rmse_t_2 , corr_t_2 ,
                                      rmse_t_2_adni , corr_t_2_adni , rmse_t_2_ucl , corr_t_2_ucl , rmse_t_2_epad , corr_t_2_epad ,
                                      rmse_t_2_adni_cn , corr_t_2_adni_cn , rmse_t_2_adni_mci , corr_t_2_adni_mci , rmse_t_2_adni_ad, 
                                      corr_t_2_adni_ad , exp_ind))
    run_ind = run_ind + 1

  }
  
  cols = grep( paste0("clean_hv_bilateral_nc_gpc_icv_",mod_num,"_[0-9][0-9][0-9][0-9]$") , names(epad_all) , value = TRUE)

  mean = rowMeans(ukb_all[ , paste0(cols,"_train")])
  sd   = rowSds(as.matrix(ukb_all[ , paste0(cols,"_train")]))
  ukb_all[ , paste0("train_clean_hv_bilateral_nc_gpc_icv_",mod_num,"_mean")] = mean
  ukb_all[ , paste0("train_clean_hv_bilateral_nc_gpc_icv_",mod_num,"_sd")] = sd

  mean = rowMeans(ukb_all[ , paste0(cols,"_test")])
  sd   = rowSds(as.matrix(ukb_all[ , paste0(cols,"_test")]))
  ukb_all[ , paste0("test_clean_hv_bilateral_nc_gpc_icv_",mod_num,"_mean")] = mean
  ukb_all[ , paste0("test_clean_hv_bilateral_nc_gpc_icv_",mod_num,"_sd")] = sd

  mean = rowMeans(adni_all[ , cols])
  sd   = rowSds(as.matrix(adni_all[ , cols]))
  adni_all[ , paste0("clean_hv_bilateral_nc_gpc_icv_",mod_num,"_mean")] = mean
  adni_all[ , paste0("clean_hv_bilateral_nc_gpc_icv_",mod_num,"_sd")] = sd

  mean = rowMeans(ucl_all[ , cols])
  sd   = rowSds(as.matrix(ucl_all[ , cols]))
  ucl_all[ , paste0("clean_hv_bilateral_nc_gpc_icv_",mod_num,"_mean")] = mean
  ucl_all[ , paste0("clean_hv_bilateral_nc_gpc_icv_",mod_num,"_sd")] = sd

  mean = rowMeans(epad_all[ , cols])
  sd   = rowSds(as.matrix(epad_all[ , cols]))
  epad_all[ , paste0("clean_hv_bilateral_nc_gpc_icv_",mod_num,"_mean")] = mean
  epad_all[ , paste0("clean_hv_bilateral_nc_gpc_icv_",mod_num,"_sd")] = sd

  exp_ind = exp_ind + 1
  mod_num = mod_num + 1
}

ABSTRACTION_RESULTS = setNames(ABSTRACTION_RESULTS , c("run_ind", "exp_ind" , "experiment", "dataset" , "seed" , "RMSE" , "ASSOC", "RMSE_t" , "ASSOC_t","RMSE_t2" , "ASSOC_t2" ,
                                                       "RMSE_ADNI" , "ASSOC_ADNI", "RMSE_UCL" , "ASSOC_UCL", "RMSE_EPAD" , "ASSOC_EPAD" ,
                                                       "RMSE_ADNI_CN" , "ASSOC_ADNI_CN","RMSE_ADNI_MCI" , "ASSOC_ADNI_MCI","RMSE_ADNI_AD" , "ASSOC_ADNI_AD",
                                                       "exp_grp"))

ABSTRACTION_RESULTS$exp_grp = ABSTRACTION_RESULTS$exp_ind 
ABSTRACTION_RESULTS[ grep( "^gpr.*PRS_TH.*" , ABSTRACTION_RESULTS$experiment) , "exp_grp"] = floor(mean(ABSTRACTION_RESULTS[ grep( "^gpr.*PRS_TH.*" , ABSTRACTION_RESULTS$experiment) , "exp_ind"]))


means = colMeans(as.matrix(ukb_all[ , prs_columns_nc_gpc]))
sds =  colSds(as.matrix(ukb_all[ , prs_columns_nc_gpc]))                 
ukb_all[ , paste0(prs_columns_nc_gpc,"_zscore")] =  sweep( sweep(ukb_all[ , prs_columns_nc_gpc]  , 2 , means) , 2 , sds , "/")

pred_col = "clean_hv_bilateral_nc_gpc_icv"
COVS =list(c("age","sex","ICV_nc_gpc","PRS_LASSO"))
LASSO_WEIGHTS = data.frame(matrix(ncol = length(prs_columns) + 1 , nrow = 0))
for(covs in COVS)
  for ( seed in SEEDS ) {
    covs_tmp = covs
    covs_tmp = c( covs[-length(covs)] , paste0(prs_columns_nc_gpc,"_zscore"))
    grab_cols = c("id" , covs_tmp , pred_col)
    ukb_rows = complete.cases(ukb_all[,grab_cols])
    
    lasso = LASSO_SCORE(ukb_all[ukb_rows , ] , ukb_all[ukb_rows,pred_col] , paste0(prs_columns_nc_gpc,"_zscore") , print_weights = TRUE , seed = seed)
    print(lasso$weights)
    LASSO_WEIGHTS <- rbind(LASSO_WEIGHTS , c(seed,lasso$weights))
  }
names(LASSO_WEIGHTS) = c("seed", paste0(prs_columns_nc_gpc,"_zscore"))
LASSO_WEIGHTS = cbind( experiment = "ukb_intersect_pgs", LASSO_WEIGHTS)

p = ggplot(data = melt(LASSO_WEIGHTS, c("experiment","seed")) , aes(x=variable, y=value )) + geom_boxplot(aes(fill=experiment)) + 
  scale_fill_discrete(guide = "none" ) +
  scale_x_discrete(labels= do.call(rbind, strsplit(prs_columns_nc_gpc , "_"))[,4] ) + 
  xlab("PGS THRESHOLD") +
  ylab("LASSO WEIGHT") + 
  ggtitle("ukb_intersect_pgs")

pdf(file = paste0(plt_folder,"Supplementary Figure - Lasso Weights",".pdf") , width = 10 , height = 10)
p
dev.off()



save(ABSTRACTION_RESULTS , file = "~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/ABSTRACTION_RESULTS 03_02_2024.RData")
#save(DETAILED_RESULTS , file = "~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/LIN_DETAILED_RESULTS 27_01_2024.RData")

write.csv(ukb_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/with_imaging_40k/PREPROCESSED_UKB_03_02_2024_RESULTS.csv")
write.csv(adni_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data//PREPROCESSED_ADNI_03_02_2024_RESULTS.csv")
write.csv(ucl_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/ucl46_data//PREPROCESSED_I46_03_02_2024_RESULTS.csv")
write.csv(epad_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/PREPROCESSED_EPAD_03_02_2024_RESULTS.csv")

load(file = "~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/ABSTRACTION_RESULTS 03_02_2024.RData")
#load(file = "~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/LIN_DETAILED_RESULTS 27_01_2024.RData")
#load(file = "~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/GPS_RESULTS 27_01_2024.RData")
 
ukb_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/with_imaging_40k/PREPROCESSED_UKB_03_02_2024_RESULTS.csv")
ucl_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/ucl46_data//PREPROCESSED_I46_03_02_2024_RESULTS.csv")
adni_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data//PREPROCESSED_ADNI_03_02_2024_RESULTS.csv")
epad_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/PREPROCESSED_EPAD_03_02_2024_RESULTS.csv")


# A few figures to show how the GPCs compare across datasets
par(mfrow = c(2,2))
make_GPC_plot("snp_weights_gpc_EA_1" , "snp_weights_gpc_EA_2" , ylim = c(-0.02 , 0.06))
make_GPC_plot("snp_weights_gpc_NA_1" , "snp_weights_gpc_NA_2" , legend_loc = "topright")
make_GPC_plot("snp_weights_gpc_NA_1" , "snp_weights_gpc_NA_3", legend_loc = "topright")
make_GPC_plot("snp_weights_gpc_NA_2" , "snp_weights_gpc_NA_3")

make_GPC_plot = function( gpc_1 = "snp_weights_gpc_EA_1" , gpc_2 = "snp_weights_gpc_EA_2" , ylim = NULL , legend_loc = "topleft"){
  
  plot(ukb_all[,gpc_1] , ukb_all[,gpc_2] , col = scales::alpha("black",0.5) , pch = 1 , xlab = gpc_1 , ylab = gpc_2 , ylim = ylim)
  points(adni_all[ adni_all$VISCODE.x == 'bl' , gpc_1] , adni_all[ adni_all$VISCODE.x == 'bl' , gpc_2] , col = scales::alpha("blue",0.5) , pch = 15)
  points(epad_all[,gpc_1] , epad_all[,gpc_2] , col = scales::alpha("red",0.5) , pch = 17)
  points(ucl_all[,gpc_1] , ucl_all[,gpc_2] , col = scales::alpha("yellow",0.5) , pch = 16)
  
  points(adni_all[ adni_all$VISCODE.x == 'bl' & adni_all$prob_EU < 0.8 , gpc_1] ,
         adni_all[ adni_all$VISCODE.x == 'bl' & adni_all$prob_EU < 0.8, gpc_2]  , pch = 0)
  
  points(epad_all[ epad_all$prob_EU < 0.8 , gpc_1] , epad_all[ epad_all$prob_EU < 0.8, gpc_2] , pch = 2)
  
  points(ucl_all[ ucl_all$prob_EU < 0.8 , gpc_1] , ucl_all[ ucl_all$prob_EU < 0.8, gpc_2]  , pch = 1)
  

  legend(legend_loc, legend=c("UKB", "ADNI" , "EPAD" , "I46"),
         col=c( scales::alpha("black",0.5) ,  scales::alpha("blue",0.5) , scales::alpha("red",0.5) , scales::alpha("yellow",0.5) ), pch=c(1,15,17,16), cex=0.8)
}


# Model comparison figure
par(mfrow = c(2,3) )
lins = grep("^lin" , ABSTRACTION_RESULTS$experiment )



t = ABSTRACTION_RESULTS[ ABSTRACTION_RESULTS$exp_grp %in% c(3,10,17) , ] 
labels = c("ASI" ,"ASIP-single" , "ASIP-lasso")

p1 = make_boxplot_2( t , c("RMSE", "RMSE_t" , "RMSE_t2" ) ,
                     labs = labels , title = "(A) RMSE Performance" ,
                     xlab = c("TRAINING" , "TESTING" , "OUT OF SAMPLE") ,
                     ylab = "Root Mean Square Error (mm3)" ,
                     rm_legend = TRUE)

p2 = make_boxplot_2( t , c("ASSOC", "ASSOC_t" , "ASSOC_t2" ) , 
                     labs = labels , title = "(B) Correlation Performance" ,
                     xlab = c("TRAINING" , "TESTING" , "OUT OF SAMPLE") ,
                     ylab = "Correlation to HV (%)" ,
                     legend_name = "zscore group" ,  rm_legend = FALSE)


pdf(file = paste0(plt_folder,"HV_Figure 4 - 1",".pdf") , width = 10 , height = 7)

p1 + p2 

dev.off()

p1 = make_boxplot_2( t , c("RMSE", "RMSE_t" ,"RMSE_ADNI" , "RMSE_UCL" , "RMSE_EPAD" ) ,
                     labs = labels , title = "(A) RMSE Performance" ,
                     xlab = c("TRAINING" , "TESTING" , "ADNI" , "I46" , "EPAD") ,
                     ylab = "Root Mean Square Error (mm3)" ,
                     rm_legend = TRUE)

p2 = make_boxplot_2( t , c("ASSOC", "ASSOC_t" , "ASSOC_ADNI" , "ASSOC_UCL" , "ASSOC_EPAD" ) , 
                     labs = labels , title = "(B) Correlation Performance" ,
                     xlab = c("TRAINING" , "TESTING" , "ADNI" , "I46" , "EPAD") ,
                     ylab = "Correlation to HV (%)" ,
                     legend_name = "zscore group" ,  rm_legend = FALSE)


pdf(file = paste0(plt_folder,"HV_Figure 4 - 2",".pdf") , width = 10 , height = 7)

p1 + p2 

dev.off()


t = ABSTRACTION_RESULTS 
labels = c("A" ,"AS" , "ASI" , paste0( "ASIP_" , gsub( "INTERSECT_PRS_TH_" , "" , prs_columns) ) , "ASIP_lasso")
p1 = make_boxplot_2( t , grp_col =  "exp_ind" , c("RMSE", "RMSE_t" , "RMSE_t2" ) ,
                     labs = labels , title = "(A) RMSE Performance" ,
                     xlab = c("TRAINING" , "TESTING" , "OUT OF SAMPLE") ,
                     ylab = "Root Mean Square Error (mm3)" ,
                     rm_legend = TRUE)

p2 = make_boxplot_2( t , grp_col =  "exp_ind" , c("ASSOC", "ASSOC_t" , "ASSOC_t2" ) , 
                     labs = labels , title = "(B) Correlation Performance" ,
                     xlab = c("TRAINING" , "TESTING" , "OUT OF SAMPLE" ) ,
                     ylab = "Correlation to HV (%)" ,
                     legend_name = "zscore group" ,  rm_legend = FALSE)

pdf(file = paste0(plt_folder,"HV_Figure 4 - 3",".pdf") , width = 15 , height = 11)

p1 + p2 

dev.off()



p1 = make_boxplot_2( t , grp_col =  "exp_ind" , c("RMSE_ADNI" , "RMSE_UCL" , "RMSE_EPAD" ) ,
                     labs = labels , title = "(A) RMSE Performance" ,
                     xlab = c("ADNI" , "I46" , "EPAD") ,
                     ylab = "Root Mean Square Error (mm3)" ,
                     rm_legend = TRUE)

p2 = make_boxplot_2( t , grp_col =  "exp_ind" , c( "ASSOC_ADNI" , "ASSOC_UCL" , "ASSOC_EPAD" ) , 
                     labs = labels , title = "(B) Correlation Performance" ,
                     xlab = c("ADNI" , "I46" , "EPAD") ,
                     ylab = "Correlation to HV (%)" ,
                     legend_name = "zscore group" ,  rm_legend = FALSE)

pdf(file = paste0(plt_folder,"HV_Figure 4 - 4",".pdf") , width = 15 , height = 11)

p1 + p2 

dev.off()


p1 = make_boxplot_2( t , grp_col =  "exp_ind" , c("RMSE" , "RMSE_t" , "RMSE_ADNI_CN" , "RMSE_ADNI_MCI" , "RMSE_ADNI_AD" ) ,
                     labs = labels , title = "(A) RMSE Performance" ,
                     xlab = c("UKB_TRAIN" , "UKB_TEST" ,"ADNI_CN" , "ADNI_MCI" , "ADNI_AD") ,
                     ylab = "Root Mean Square Error (mm3)" ,
                     rm_legend = TRUE)

p2 = make_boxplot_2( t , grp_col =  "exp_ind" , c("ASSOC" , "ASSOC_t" , "ASSOC_ADNI_CN" , "ASSOC_ADNI_MCI" , "ASSOC_ADNI_AD" ) , 
                     labs = labels , title = "(B) Correlation Performance" ,
                     xlab = c("UKB_TRAIN","UKB_TEST", "ADNI_CN" , "ADNI_MCI" , "ADNI_AD") ,
                     ylab = "Correlation to HV (%)" ,
                     legend_name = "zscore group" ,  rm_legend = FALSE)

pdf(file = paste0(plt_folder,"HV_Figure 4 - 5",".pdf") , width = 15 , height = 11)

p1 + p2 

dev.off()


adni_eval_cols = c("MMSE" , "CDRSB" ,"EcogSPTotal" , "ADAS13" ,  "MOCA" , "GDTOTAL" , "PTAU" , "TAU" , "ABETA")


adni_eval_cols_memory = c("EcogSPMem" , "MMSE_RECALL" , "ADAS_MEMORY" , "MOCA_MEMORY" , "CDMEMORY")
adni_eval_cols_language = c("EcogSPLang" , "MMSE_LANGUAGE" , "ADAS_LANGUAGE" , "MOCA_LANGUAGE")
adni_eval_cols_orientation = c("MMSE_ORIENTATION" , "ADAS_ORIENTATION" , "MOCA_ORIENTATION" , "CDORIENT")
adni_eval_cols_visuospatial = c("EcogSPVisspat" , "ADAS_VISUOSPATIAL" , "MOCA_VISUOSPATIAL")
adni_eval_cols_execfunc = c("EcogSPEF","MMSE_ATTCALC" , "MOCA_EF" , "CDJUDGE")


epad_eval_cols = c( "mmse_total" , "cdr_sum_of_box" , "rbans_total_scale" ,  "flanker_score" , "dotcount_total", 
                    "fav_delay_total_correct", "fav_learn_r2_correct", "fav_learn_r1_correct" , "fav_rec_total_correct" ,
                    "gds_total" , "ptau_result" , "ttau_result" , "abeta_1_42_result" , "st_trial_total" , "fms_med_total")


epad_eval_cols_memory = c("rbans_delayed_memory_index" , "rbans_immediate_memory_index" , "mmse_recall_total" ,  "cdr_memory") 
epad_eval_cols_language = c("rbans_language_index" ,"mmse_naming_total")
epad_eval_cols_orientation = c("mmse_orientation_total" , "cdr_orientation")
epad_eval_cols_visuospatial = c("flanker_score" , "fms_med_total" , "st_trial_total" , "rbans_visuo_constructional_index")
epad_eval_cols_execfunc = c("mmse_calc_total" , "cdr_judgement")



# rename the years since bl and id's in epad and adni to match.
adni_all$visit_year = adni_all$Years.bl
adni_all$visit = adni_all$VISCODE2
epad_all$id = epad_all$patient_id

# calculate the baseline values for the performance metrics in epads. 
fill_cols = paste0("clean_hv_bilateral_nc_gpc_icv_", 1:17,"_mean" )

for(eval_col in c(epad_eval_cols,epad_eval_cols_memory , epad_eval_cols_language , epad_eval_cols_orientation , epad_eval_cols_visuospatial , epad_eval_cols_execfunc)){
  eval_col_bl = paste0(eval_col,".bl")
  epad_all[,eval_col_bl] = NA
  epad_all[ epad_all$visit_year == 0, eval_col_bl] = epad_all[ epad_all$visit_year == 0 , eval_col]
  epad_all = as.data.frame(epad_all %>% group_by(id) %>% mutate(across( all_of(eval_col_bl)  , ~ifelse(!is.na(.), ., first(na.omit(.))))) )
}

for(eval_col in c(adni_eval_cols , fill_cols)){
  eval_col_bl = paste0(eval_col,".bl")
  adni_all[,eval_col_bl] = NA
  adni_all[ adni_all$visit_year == 0 & !is.na(adni_all[,eval_col]) , eval_col_bl] = adni_all[ adni_all$visit_year == 0 & !is.na(adni_all[,eval_col]), eval_col]
  adni_all = as.data.frame(adni_all %>% group_by(id) %>% mutate(across( all_of(eval_col_bl)  , ~ifelse(!is.na(.), ., first(na.omit(.))))) )
}

adni_all[ !is.na(adni_all$ICV_nc_gpc) , "ICV_nc_gpc.bl"] = adni_all[ !is.na(adni_all$ICV_nc_gpc) , "ICV_nc_gpc"]
epad_all[ !is.na(epad_all$ICV_nc_gpc) , "ICV_nc_gpc.bl"] = epad_all[ !is.na(epad_all$ICV_nc_gpc) , "ICV_nc_gpc"]

adni_all <- data.frame( adni_all %>% arrange(id, visit) %>% group_by(id) %>% fill(ICV_nc_gpc.bl, .direction = "downup") %>% ungroup() )
epad_all <- data.frame( epad_all %>% arrange(id, visit) %>% group_by(id) %>% fill(ICV_nc_gpc.bl, .direction = "downup") %>% ungroup() )

adni_all[ adni_all$visit_year == 0 , "ICV_nc_gpc"] = adni_all[ adni_all$visit_year == 0 , "ICV_nc_gpc.bl"]
epad_all[ epad_all$visit_year == 0 , "ICV_nc_gpc"] = epad_all[ epad_all$visit_year == 0 , "ICV_nc_gpc.bl"]

for( eval_col in c(adni_eval_cols , adni_eval_cols_memory , adni_eval_cols_language , adni_eval_cols_visuospatial , adni_eval_cols_orientation , adni_eval_cols_execfunc) ){
  eval_col_bl =  paste0(eval_col,".bl")
  adni_all[ (adni_all$VISCODE2 %in% c("bl")) & !is.na(adni_all[, eval_col]) , eval_col_bl] = adni_all[ (adni_all$VISCODE2 %in% c("bl")) & !is.na(adni_all[, eval_col]),  eval_col]
  adni_all = as.data.frame(adni_all %>% group_by(id) %>% mutate(across( all_of( eval_col_bl)  , ~ifelse(!is.na(.), ., first(na.omit(.))))) )
  adni_all[ (adni_all$VISCODE2 %in% c("sc")) & !is.na(adni_all[, eval_col]) , eval_col_bl] = adni_all[ (adni_all$VISCODE2 %in% c("sc")) & !is.na(adni_all[, eval_col]),  eval_col]

  adni_all = data.frame(adni_all %>% group_by(id) %>% mutate( eval_col_bl = if_else(VISCODE2 == "bl" & is.na(eval_col_bl),
                                                                                    first(eval_col_bl[VISCODE2 == "sc" & !is.na(eval_col_bl)]),
                                                                                    eval_col_bl) ) %>% ungroup() )
  adni_all[ (adni_all$visit_year %in% c(0)) & is.na(adni_all[ , eval_col]) , eval_col] =  adni_all[ (adni_all$visit_year %in% c(0)) & is.na(adni_all[ , eval_col]) , eval_col_bl]
  
}


adni_all[ adni_all$VISCODE2 %in% c("bl") , "visit_year"] = 0

adni_ADAS_2 = read.csv("Downloads/Neuropsychological (2)/ADAS_ADNIGO23_27Feb2024.csv")
adni_ADAS_2$temp_MEMORY = rowSums(adni_ADAS_2[ , c("Q1SCORE","Q4SCORE","Q8SCORE")])
adni_ADAS_2$temp_LANGUAGE = rowSums(adni_ADAS_2[ , c("Q10SCORE","Q11SCORE","Q12SCORE")])
adni_ADAS_2$temp_VISUOSPATIAL = adni_ADAS_2$Q3SCORE
adni_ADAS_2$temp_ORIENTATION = adni_ADAS_2$Q7SCORE
adni_ADAS_2$temp_TOTAL11 = adni_ADAS_2$TOTSCORE
adni_ADAS_2$temp_TOTAL13 = adni_ADAS_2$TOTAL13
adni_ADAS_2$temp_VISDATE = adni_ADAS_2$VISDATE

adni_all = merge(adni_all, adni_ADAS_2[ adni_ADAS_2$PTID %in% adni_all$PTID , c("PTID","VISCODE2","temp_MEMORY" , "temp_LANGUAGE" , "temp_VISUOSPATIAL" , "temp_ORIENTATION" , "temp_TOTAL11" , "temp_TOTAL13" , "temp_VISDATE")] , by = c("PTID" , "VISCODE2") , all = TRUE) 

adni_all[ is.na(adni_all$ADAS11) , "ADAS11" ] = adni_all[ is.na(adni_all$ADAS11) , "temp_TOTAL11" ]
adni_all[ is.na(adni_all$ADAS13) , "ADAS13" ] = adni_all[ is.na(adni_all$ADAS13) , "temp_TOTAL13" ]

adni_all[ is.na(adni_all$ADAS_MEMORY) , "ADAS_MEMORY" ] = adni_all[ is.na(adni_all$ADAS_MEMORY) , "temp_MEMORY" ]
adni_all[ is.na(adni_all$ADAS_LANGUAGE) , "ADAS_LANGUAGE" ] = adni_all[ is.na(adni_all$ADAS_LANGUAGE) , "temp_LANGUAGE" ]
adni_all[ is.na(adni_all$ADAS_VISUOSPATIAL) , "ADAS_VISUOSPATIAL" ] = adni_all[ is.na(adni_all$ADAS_VISUOSPATIAL) , "temp_VISUOSPATIAL" ]
adni_all[ is.na(adni_all$ADAS_ORIENTATION) , "ADAS_ORIENTATION" ] = adni_all[ is.na(adni_all$ADAS_ORIENTATION) , "temp_ORIENTATION" ]
adni_all[ is.na(adni_all$EXAMDATE) , "EXAMDATE" ] = adni_all[ is.na(adni_all$EXAMDATE) , "temp_VISDATE" ]

# all baseline and characteristic columns should remain the same, so make sure they are filled in with the same values for every participant
colunms_that_should_be_the_same = c( "id" , "RID" , "sex" , "PTGENDER" , "PTEDUCAT" , "PTETHCAT" , 
                                    "PTRACCAT" , "PTMARRY" , "APOE4" , "SITE" , 
                                    "COLPROT" , "ORIGPROT" , "AGE" , 
                                    grep( "snp_weig" , names(adni_all) , value = TRUE) ,
                                    grep( "INTERSECT" , names(adni_all) , value = TRUE) , 
                                    grep( "ICV" , names(adni_all) , value = TRUE) ,
                                    grep( ".bl$" , names(adni_all) , value = TRUE) ,
                                    grep( "prob_" , names(adni_all) , value = TRUE))

adni_all = as.data.frame(adni_all %>% group_by(PTID) %>% mutate(across(colunms_that_should_be_the_same, ~ { if (is.Date(.)) {
  coalesce(., na.omit(.)[1])
} else {
  ifelse(!is.na(.), ., first(na.omit(.)))
}})))

# we can now recalculate the age and replace the Years.bl and Month.bl with their intended values.
adni_all$age <- adni_all$AGE + time_length( difftime(as.Date(adni_all$EXAMDATE) , as.Date(adni_all$EXAMDATE.bl) ) , "years")
adni_all$Years.bl <- time_length( difftime(as.Date(adni_all$EXAMDATE) , as.Date(adni_all$EXAMDATE.bl) ) , "years")
adni_all$Month.bl <- time_length( difftime(as.Date(adni_all$EXAMDATE) , as.Date(adni_all$EXAMDATE.bl) ) , "months")
adni_all$VISCODE_NUM = as.numeric( str_replace( str_replace( str_replace( adni_all$VISCODE2 , "bl" , "0")  , "sc" , "0" ) , "m" , "") )


#write.csv(adni_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data//PREPROCESSED_ADNI_29_02_2024_RESULTS_EVAL.csv")
#write.csv(epad_all , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/PREPROCESSED_EPAD_29_02_2024_RESULTS_EVAL.csv")

adni_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data//PREPROCESSED_ADNI_29_02_2024_RESULTS_EVAL.csv")
epad_all = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/PREPROCESSED_EPAD_29_02_2024_RESULTS_EVAL.csv")

for( no_dementia in c(FALSE , TRUE)){
  P_TABLE = data.frame(matrix(ncol = 8, nrow = 0))
  P_TABLE = setNames(P_TABLE , c("model-shift" , "table" , "test_cat" , "test_name", "slope" , "p-val" , "slope-no-extreme","p-val-no-extreme"))
  
  for(j in c(1:16))
    for(table_pref in c("adni" , "epad")){
      table = get(paste0(table_pref,"_all"))
      if(no_dementia) table = table[ table$DX.bl != 5 , ]
      table = table[ table$visit_year == 0 , ]
      for( category in c("",paste0("_",c("memory", "language" , "orientation" , "visuospatial" , "execfunc"))) ){
        #for( category in c("") ){
        lm_results = data.frame(matrix(ncol = 4, nrow = 0))
        for(eval_col in get( paste0(table_pref,"_eval_cols",category)) ){
          eval_col = paste0(eval_col,".bl")
          if (sum(complete.cases(table[ , c(eval_col ,
                                            paste0("clean_hv_bilateral_nc_gpc_icv_17_mean.bl"),
                                            paste0("clean_hv_bilateral_nc_gpc_icv_",j,"_mean.bl"))])) > 0 ){
            
            X =  shift_model( table , eval_col , "clean_hv_bilateral_nc_gpc_icv_17_mean.bl" , paste0("clean_hv_bilateral_nc_gpc_icv_",j,"_mean.bl") , plt = FALSE)
            lm_results = rbind( lm_results , X )
            P_TABLE = rbind(P_TABLE , list( paste("model 17 vs",j) , table_pref , category , eval_col , X[[1]] , X[[2]] , X[[3]] , X[[4]]) )
          }
        }
        if(nrow(lm_results) > 0){
          print(paste( "17 vs" , j , "combining:", paste0(table_pref,"_eval_cols",category)))
          setNames(lm_results , c("slope" , "p-val" , "slope-no-extreme","p-val-no-extreme"))
          P_TABLE = setNames(P_TABLE , c("model-shift" , "table" , "test_cat" , "test_name", "slope" , "p-val" , "slope-no-extreme","p-val-no-extreme"))
          rs = cor(get(paste0(table_pref,"_all"))[ get(paste0(table_pref,"_all"))[ , "visit_year"] == 0 , get( paste0(table_pref,"_eval_cols",category)) ] , use="complete.obs")
          ps = lm_results[ , 2]
          ps_no_ex = lm_results[ , 4]
          com_adj = stouffer( ps , adjust = "gao" , R = rs) 
          com_ps = com_adj$p
          com_m = com_adj$m
          com_adj_no_ex = stouffer( ps_no_ex , adjust = "gao" , R = rs)
          com_ps_no_ex = com_adj_no_ex$p
          com_m_no_ex = com_adj_no_ex$m
          means = c(com_m , com_ps , com_m_no_ex , com_ps_no_ex)
          P_TABLE = rbind(P_TABLE , c( paste("model 17 vs" , j) , table_pref , category , "COMBINED PVALUES" , means ) )
          
        }
      }
    }
  
  if(no_dementia){
    write.csv(P_TABLE , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/P_TABLE_NO_DEMENTIA_04_03_2024.csv")
  }else{
    write.csv(P_TABLE , "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/P_TABLE_04_03_2024.csv")
  }
}

P_TABLE = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/P_TABLE_04_03_2024.csv")
P_TABLE_NO_DEMENTIA = read.csv("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/P_TABLE_NO_DEMENTIA_04_03_2024.csv")


other_models = c(3)
num_other_models = length(other_models)
results_table =  data.frame(matrix(ncol = (5 + 8) , nrow = 0))

outlier_threshold = 11
for( j in other_models )
  for( table_pref in c( "epad")){
    
    table = get(paste0(table_pref,"_all"))
    #table = table[ table$DX.bl != 5 , ]
    
    for( category in c("",paste0("_",c("memory", "language" , "orientation" , "visuospatial" , "execfunc"))) ){
      model_results = data.frame(matrix(ncol = 8, nrow = 0))
      for( eval_col in get( paste0(table_pref,"_eval_cols",category)) ){
        print(paste0("checking ",eval_col," in ",table_pref))
        results_row = c()
        
        table$shift = ( pnorm(unlist(table[ , "clean_hv_bilateral_nc_gpc_icv_17_mean.bl"])) -
                          pnorm(unlist(table[ , paste0("clean_hv_bilateral_nc_gpc_icv_",j,"_mean.bl") ])) )* 100
        
        table$shift = table$shift + 100
        
        # rescale the ICV for the mixed effects model.
        table$ICV_new = (table$ICV_nc_gpc.bl - mean(table$ICV_nc_gpc.bl , na.rm = TRUE)) / sd(table$ICV_nc_gpc.bl , na.rm = TRUE)
        
        
        # exclude all the subjects where we only have one time point for this evaluation metric.
        # to do this: create a pivot table where each unique ID has a row of which visits they have attended, 
        # if that list if <= 1 then we have too few data points for this subject. 
        df = data.frame( table[ !is.na(table[ , eval_col]) , c("id","visit_year")] %>% group_by(id , visit_year) %>% tally() %>%
                           pivot_wider(names_from=visit_year, values_from=n) )
        d = data.frame(df[,1],rowSums(df[,-1] , na.rm = TRUE))
        
        #mean_num_visits = mean( d[,2] , na.rm = TRUE)
        mean_num_visits = mean( d[ d[,2] > 1 , 2] , na.rm = TRUE)
        
        ids_to_exclude = d[ d[,2] <= 1 , 1]
        table_temp = table[ !(table$id %in% ids_to_exclude) & !is.na(table[,eval_col]), ]
        table_temp_no_extreme = table_temp[ abs(table_temp$shift -100) < outlier_threshold , ]
        
        print(nrow(table_temp))
        print(nrow(table_temp_no_extreme))
        for( table_temp_temp in list(table_temp , table_temp_no_extreme)){
          
          LMM_effect = NA
          LMM_pvalue = NA
          anova_pval_s = NA
          anova_pval_st = NA
          
          if( sum(complete.cases(table_temp_temp[ , c(eval_col , "age" , "sex", "ICV_new", "clean_hv_bilateral_nc_gpc_icv_3_mean.bl" , "shift" , "visit_year")])) > 0 ){
            
            base_model = lmer( as.formula( paste0(eval_col," ~ age + sex + ICV_new + (clean_hv_bilateral_nc_gpc_icv_3_mean.bl + shift)*visit_year + (1 + visit_year | id )"))   , data = table_temp_temp ,
                               control = lmerControl(optimizer ="nmkbw") , REML = FALSE)
            temp = summary(base_model)
            LMM_effect = temp$coefficients[9,1] * sd(table_temp_temp[,eval_col] , na.rm = TRUE)
            LMM_pvalue = temp$coefficients[9,5]
            
            
            no_shift_model = lmer( as.formula(paste0(eval_col," ~  age + sex + ICV_new + (clean_hv_bilateral_nc_gpc_icv_3_mean.bl * visit_year) +  (1 + visit_year | id )"))  , data = table_temp_temp,
                                   control = lmerControl(optimizer ="nmkbw") , REML = FALSE)
            no_time_shift_model = lmer( as.formula(paste0(eval_col," ~ age + sex + ICV_new + (clean_hv_bilateral_nc_gpc_icv_3_mean.bl * visit_year) + shift +  (1 + visit_year | id)"))  , data = table_temp_temp,
                                        control = lmerControl(optimizer ="nmkbw") , REML = FALSE)
            
            anova_mod_s = anova(no_shift_model , base_model )
            anova_pval_s = scientific(anova_mod_s$`Pr(>Chisq)`[2] ) 
            
            anova_mod_st = anova(no_time_shift_model , base_model )
            anova_pval_st =  scientific(anova_mod_st$`Pr(>Chisq)`[2] ) 
            
            
            
          } else {print("not enough visits")}
          results_row = c(results_row , c( LMM_effect , LMM_pvalue , anova_pval_s , anova_pval_st) )
        }
        results_table = rbind( results_table , c( paste("model 17 vs" , j) , table_pref, category, eval_col , mean_num_visits , results_row) )
        model_results = rbind(model_results , unlist(lapply(results_row, as.numeric)))
      }
    
    if(nrow(model_results) > 0){
      print(paste( "17 vs" , j , "combining:", paste0(table_pref,"_eval_cols",category)))
      rs = cor(get(paste0(table_pref,"_all"))[ get(paste0(table_pref,"_all"))[ , "visit_year"] == 0 , get( paste0(table_pref,"_eval_cols",category)) ] , use="complete.obs")
      ps = model_results[ , 4]
      ps_no_ex = model_results[ , 8]
      com_adj = stouffer( ps , adjust = "gao" , R = rs) 
      com_ps = com_adj$p
      com_m = com_adj$m
      com_adj_no_ex = stouffer( ps_no_ex , adjust = "gao" , R = rs)
      com_ps_no_ex = com_adj_no_ex$p
      com_m_no_ex = com_adj_no_ex$m
      means = c(com_m , 0 , 0, com_ps , com_m_no_ex , 0 , 0, com_ps_no_ex)
      results_table = rbind(results_table , c( paste("model 17 vs" , j) , table_pref , category , "COMBINED PVALUES" , 0 , means ) )
    }
    }
  }


results_table = setNames(results_table , c("Models Compared" , "table", "eval_col" , "test_category" , "num_visits" ,
                                            "Effect Size" , "LMM P-value" , "no-shift vs shift" , "shift-time vs no-shift-time",
                                            "Effect Size No-Extreme" , "LMM P-value No-Extreme" ,
                                            "no-shift vs shift No-Extreme" , "shift-time vs no-shift-time No-Extreme" ))


#write.csv(results_table , "OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/LONGIDUTINAL_MODELS_RESULTS_05_03_2024.csv")

write.csv(results_table , "OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/LONGIDUTINAL_MODELS_NO_DEMENTIA_RESULTS_05_04_2024.csv")


results_table = read.csv("OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/LONGIDUTINAL_MODELS_RESULTS_29_02_2024.csv")
results_table_no_dementia = read.csv("OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/LONGIDUTINAL_MODELS_NO_DEMENTIA_RESULTS_29_02_2024.csv")


sig_eval_cols = results_table[ !is.na(results_table[,6])  , c("table","eval_col") ]

sig_eval_cols = results_table[ !is.na(results_table[ , 17]) & (results_table[ , 17] < 0.05) & results_table[ , 4] %in% c("") , 2:3 ]

# sig_eval_cols = data.frame( table = rep("adni",1) , eval_col = c("MMSE") )

for(i in 1:nrow(sig_eval_cols)){
  for(j in c(3,4,9,16)){
    eval_table = get(paste0(sig_eval_cols[i,1],"_all")) 
    eval_col = sig_eval_cols[i,2]
    test = tryCatch( {
    assign( paste0("p_long_" ,sig_eval_cols[i,1] , "_" , eval_col,"_17_",j), plot_longitudinal( eval_table, eval_col, "clean_hv_bilateral_nc_gpc_icv_17_mean.bl" ,paste0("clean_hv_bilateral_nc_gpc_icv_",j,"_mean.bl") , ylab = "") )
    assign( paste0("p_cross_",sig_eval_cols[i,1] , "_" , eval_col,"_17_",j), shift_model( eval_table[ eval_table$visit_year %in% c(0) , ] , paste0(eval_col,".bl") , "clean_hv_bilateral_nc_gpc_icv_17_mean.bl" , paste0("clean_hv_bilateral_nc_gpc_icv_",j,"_mean.bl") , plt = TRUE , legend = TRUE , 
                                                                                          ylab = paste("residual",eval_col) ))
                
    }, error = function(e){
      print(e)
    })
      
  }
}


p_cross_adni_MMSE  = shift_model(adni_all , "MMSE.bl" , "clean_hv_bilateral_nc_gpc_icv_17_mean.bl" , "clean_hv_bilateral_nc_gpc_icv_3_mean.bl" , plt = TRUE , title = "(A)" , ylab = "Residual MMSE" , legend = TRUE , ylimit = c(0,40))
p_cross_adni_CDRSB = shift_model(adni_all , "CDRSB.bl" , "clean_hv_bilateral_nc_gpc_icv_17_mean.bl" , "clean_hv_bilateral_nc_gpc_icv_3_mean.bl" , plt = TRUE , title = "(B)" , ylab = "Residual CDRSB" , legend = TRUE, ylimit = c(-5,20))
p_long_adni_MMSE   = plot_longitudinal( adni_all, "MMSE" , "clean_hv_bilateral_nc_gpc_icv_17_mean.bl" , "clean_hv_bilateral_nc_gpc_icv_3_mean.bl",  ylab = "" , title = "(C)" , ylimit = c(0,40) )
p_long_adni_CDRSB  = plot_longitudinal( adni_all, "CDRSB" , "clean_hv_bilateral_nc_gpc_icv_17_mean.bl" , "clean_hv_bilateral_nc_gpc_icv_3_mean.bl", ylab = "" , title = "(D)" , ylimit = c(-5,20) )
  
pdf(file = paste0(plt_folder,"HV_Figure 5 - 1",".pdf") , width = 10 , height = 10 )

p_cross_adni_MMSE + p_long_adni_MMSE + p_cross_adni_CDRSB + p_long_adni_CDRSB 

dev.off()


pdf(file = paste0(plt_folder , "HV_Figure 5 - 2",".pdf") , width = 25 , height = 20)

p_cross_adni_MMSE_17_3 + ggtitle("\ncross sectional analysis") + p_long_adni_MMSE_17_3 + ggtitle("LASSO VS NO GENETICS\nlongitudinal analysis") + 
  p_cross_adni_CDRSB_17_3 + ggtitle("\ncross sectional analysis") + p_long_adni_CDRSB_17_3 + ggtitle("\nlongitudinal analysis") + 
  p_cross_adni_ADAS13_17_3 + p_long_adni_ADAS13_17_3 + 
  p_cross_adni_MOCA_17_3 + p_long_adni_MOCA_17_3 + 
  p_cross_adni_EcogSPTotal_17_3 + p_long_adni_EcogSPTotal_17_3

dev.off()

pdf(file = paste0(plt_folder , "HV_Figure 5 - 3",".pdf") , width = 25 , height = 20)

p_cross_adni_MMSE_17_4 + ggtitle("\ncross sectional analysis") + p_long_adni_MMSE_17_4 + ggtitle("LASSO VS PGS_TH_1e-8\nlongitudinal analysis") + 
  p_cross_adni_CDRSB_17_4 + ggtitle("\ncross sectional analysis") + p_long_adni_CDRSB_17_4 + ggtitle("\nlongitudinal analysis") + 
  p_cross_adni_ADAS13_17_4 + p_long_adni_ADAS13_17_4 + 
  p_cross_adni_MOCA_17_4 + p_long_adni_MOCA_17_4 +
  p_cross_adni_EcogSPTotal_17_4 + p_long_adni_EcogSPTotal_17_4 

dev.off()

pdf(file = paste0(plt_folder , "HV_Figure 5 - 4",".pdf") , width = 25 , height = 20)

p_cross_adni_MMSE_17_9 + ggtitle("\ncross sectional analysis") + p_long_adni_MMSE_17_9 + ggtitle("LASSO VS PGS_TH_0.05\nlongitudinal analysis") + 
  p_cross_adni_CDRSB_17_9 + ggtitle("\ncross sectional analysis") + p_long_adni_CDRSB_17_9 + ggtitle("\nlongitudinal analysis") + 
  p_cross_adni_ADAS13_17_9 + p_long_adni_ADAS13_17_9 +  
  p_cross_adni_MOCA_17_9 + p_long_adni_MOCA_17_9 + 
  p_cross_adni_EcogSPTotal_17_9 + p_long_adni_EcogSPTotal_17_9

dev.off()

pdf(file = paste0(plt_folder , "HV_Figure 5 - 5",".pdf") , width = 25 , height = 20)

p_cross_adni_MMSE_17_16 + ggtitle("\ncross sectional analysis") + p_long_adni_MMSE_17_16 + ggtitle("LASSO VS PGS_TH_1\nlongitudinal analysis") + 
  p_cross_adni_CDRSB_17_16 + ggtitle("\ncross sectional analysis") + p_long_adni_CDRSB_17_16 + ggtitle("\nlongitudinal analysis") + 
  p_cross_adni_ADAS13_17_16 + p_long_adni_ADAS13_17_16 +
  p_cross_adni_MOCA_17_16 + p_long_adni_MOCA_17_16 +
  p_cross_adni_EcogSPTotal_17_16 + p_long_adni_EcogSPTotal_17_16 

dev.off()


excluded_table =  data.frame(matrix(ncol = 9 , nrow = 0))

for(table_pref in list("adni" , "epad") ){
  table = get(paste0(table_pref,"_all"))
  for( i in 1:16 ){
    table[, "shift"] = (pnorm(unlist(table[ , "clean_hv_bilateral_nc_gpc_icv_17_mean.bl"])) -
                           pnorm(unlist(table[ , paste0("clean_hv_bilateral_nc_gpc_icv_",i,"_mean.bl") ])) )* 100 
    
    for( category in c("",paste0("_",c("memory", "language" , "orientation" , "visuospatial" , "execfunc"))) )
    for( eval_col in get(paste0(table_pref,"_eval_cols",category)) ){
    
      nrow_cross = length( table[ !is.na(table$shift) & table$visit_year %in% c(0) & !is.na(table[ , paste0(eval_col,".bl")]) , "shift"] )
      nrow_cross_no_extreme = length( table[ !is.na(table$shift) & table$visit_year %in% c(0) & !is.na(table[ , paste0(eval_col,".bl")]) & abs(table$shift) < 11 , "shift"] )
      if(table_pref == "adni")
        nrow_cross_no_dementia = length( table[ !is.na(table$shift) & (table$DX.bl != 5) & table$visit_year %in% c(0) & !is.na(table[ , paste0(eval_col,".bl")]) , "shift"] )
      else
        nrow_cross_no_dementia = NA
      
      df = data.frame( table[ !is.na(table[ , eval_col]) , c("id","visit_year")] %>% group_by(id , visit_year) %>% tally() %>%
                         pivot_wider(names_from=visit_year, values_from=n) )
      d = data.frame(df[,1],rowSums(df[,-1] , na.rm = TRUE))
      
      mean_num_visits = mean( d[ d[,2] > 1 , 2] , na.rm = TRUE)
      
      ids_to_exclude = d[ d[,2] <= 1 , 1]
      
      table_temp = table[ !(table$id %in% ids_to_exclude) & !is.na(table[,eval_col ]), ]
      
      nrow_long = length( table_temp[ !is.na(table_temp$shift), "shift"] )
      nrow_long_no_extreme = length( table_temp[ !is.na(table_temp$shift) & abs(table_temp$shift) < 11 , "shift"] )
      if(table_pref == "adni")
        nrow_long_no_dementia = length( table_temp[ !is.na(table_temp$shift) & (table_temp$DX.bl != 5) , "shift"] )
      else
        nrow_long_no_dementia = NA
      
      excluded_table = rbind( excluded_table , list( table_pref , paste0("17_",i) , nrow_cross , nrow_cross_no_extreme , nrow_cross_no_dementia , eval_col , nrow_long , nrow_long_no_extreme , nrow_long_no_dementia) )
    }
  }
}

names(excluded_table) = c("table" , "models compared" , "all baseline available" , "without outliers" , "without dementia" , "evaluation column" , "all available" , "no outliers", "no dementia")

write.csv(excluded_table , "OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/OUTLIERS_01_03_2024.csv")


df =  read.csv( c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_HV_FREESURFER_40k.prsice") , sep = "\t")
df = df[ nrow(df):1 , ]

df$Num_SNP_scaled = (df$Num_SNP)/max(df$Num_SNP)
max_size = max(df$Num_SNP)

ggplot(df, aes(x = Num_SNP )) + geom_point( y = 0 ) + geom_point( y = 1) + 
  lapply( 1:14, function(i) {
    geom_curve(aes(x = df[ i , "Num_SNP"], y = 0, xend = df[ i , "Num_SNP"], yend = 1), curvature = (1 - df[ i , "Num_SNP_scaled"])/2 , ncp = 100) }
  )





for( i in 1:14){
  assign( paste0("p",i) , ggplot() + geom_arc(aes(x0 = 0, y0 = 0, r = df[i,"Num_SNP_scaled"], start = pi/2 + atan2( x = df[i,"Num_SNP"] , y= .01) , end = pi/2 - atan2( x = df[i,"Num_SNP"] , y = .01 ) )) )
  min_y = min(ggplot_build( get(paste0("p",i)) )$data[[1]]$y)
  max_y = max(ggplot_build( get(paste0("p",i)) )$data[[1]]$y)
  assign( paste0("df_poly",i) , rbind(c(0,max_y) , c(0,min_y) ,data.frame(x = ggplot_build( get(paste0("p",i)) )$data[[1]]$x,y = ggplot_build( get(paste0("p",i)) )$data[[1]]$y) , c(0,max_y)) )
}

ggplot() + lapply( 1:14, function(i) {
  geom_polygon(data = get(paste0("df_poly",i)) , aes(x,y) , alpha = 0.1 )
})

