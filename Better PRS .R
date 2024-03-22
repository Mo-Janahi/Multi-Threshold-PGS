

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
library(patchwork)
library(reshape2)
library(GGally)
library(MASS)
library(stringdist)

source("~/OneDrive - University College London/ANDRE STAGE/Project 1 PGS Nomogram/Scripts/Basic Functions.R")
source("~/OneDrive - University College London/ANDRE STAGE/Project 1 PGS Nomogram/Scripts/Plotting Functions.R")
source("~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/Preprocessing Functions.R")
source("~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/Extra Functions.R")

#####################
### PREPROCESSING ###
#####################

# UKB Table is saved as CSV file 
ukb_table <- read.csv("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_features_20200830.csv")
ukb = PREPROCESS_UKBB(ukb_table)
ukb_male <- ukb$ukb_img_ex_outliers_male
ukb_female <- ukb$ukb_img_ex_outliers_female
ukb_all <- rbind(ukb_male , ukb_female)

#ukb_no_hv = PREPROCESS_UKBB_NO_HV(ukb_table)
#ukb_male_no_hv <- ukb_no_hv$ukb_ex_outliers_male
#ukb_female_no_hv <- ukb_no_hv$ukb_ex_outliers_female

# ADNI Table is read from downloaded library
adni_table <- merge( ADNIMERGE::adnimerge , ADNIMERGE::ucsffsx51 , by = c("IMAGEUID") )
x <- PREPROCESS_ADNI(adni_table)
adni_male <- x$ADNI_male
adni_female <- x$ADNI_female
adni_all <- rbind(adni_male , adni_female)

# GREP all columns that start with "PRS_TH"
# this should be all 14 or so PGS thresholds generated
# I have standardized the column names across datasets (UKBB / ADNI / UCL)
PRS_columns <- grep("^PRS_TH" , names(ukb_male) , value = TRUE) 
# match all column names that start with "PRS_AD"
PRS_AD_columns <- grep("^AD_PRS" , names(ukb_male) , value = TRUE)

ukb_genetic_pcs = paste0("Genetic.PC.",1:20)
adni_genetic_pcs = paste0("PCA",1:10)

for(column in c( PRS_columns , PRS_AD_columns) ){
  
  ukb_male <- CORRECT_PRS( ukb_male , column , ukb_genetic_pcs )
  ukb_female <- CORRECT_PRS( ukb_female , column , ukb_genetic_pcs )
  
  adni_male <- CORRECT_PRS( adni_male , column , adni_genetic_pcs )
  adni_female <- CORRECT_PRS( adni_female , column , adni_genetic_pcs )
  
}

PRS_columns = grep("^PRS_TH.*NO_GPC_ZSCORE$" , names(ukb_male) , value = TRUE) 
PRS_AD_columns <- grep("^AD_PRS.*NO_GPC_ZSCORE$" , names(ukb_male) , value = TRUE)

PRS_columns_sup_all =  grep("^PRS_TH" , names(ukb_all) , value = TRUE) 
PRS_AD_columns_sup_all = grep("^AD_PRS_TH" , names(ukb_all) , value = TRUE) 

PRS_columns_sup = PRS_columns_sup_all[c(1,4,7,9,11,14)]
PRS_AD_columns_sup = PRS_AD_columns_sup_all[c(1,4,7,9,11,14)]

adni_all$id = paste0(adni_all$PTID , "_", adni_all$VISCODE.x)



SEEDS = 4560:4589


######################
### HV EXPERIMENTS ###
######################




# next I will show the nomogram separation that can be achieved by using only one PGS

# we down sample to half the dataset in order for results to be comparable the the lasso experiments later. 
# This is necessary since the lasso required training, and so the dataset is split in half for training/testing


# Each experiment results table will have: 
#   - Run_Index: unique for every row in the table
#   - Exp_Index: unique for every experiment that is repeated 30 times. helps with generating summary statistics
#   - Grp_Index: same as Exp_Index but additionally grouped by single/double/lasso category   
#   - Experiment Label: HV or AD experiment
#   - Dataset: UKBB, ADNI, UCL
#   - Seed: the seed that was used to randomly split the dataset in half.
#   - Sex: Male, Female
#   - Prs: either one prs, two (prs1&prs2), or lasso
#   - Diff: a measure of difference between the separated groups. 

DIFFS_TABLE = data.frame(matrix(ncol = 10, nrow = 0))
run_ind = 1
exp_ind = 1 
grp_ind = 1

for( dataset in c("UKB" , "ADNI" ) ){
  
  if(dataset == "UKB"){
    male_table  =  ukb_male[!is.na(ukb_male$AD_PRS_TH_1) & !is.na(ukb_male$clean_hv_bilateral) , ]
    female_table = ukb_female[!is.na(ukb_female$AD_PRS_TH_1) & !is.na(ukb_female$clean_hv_bilateral) , ]
    AGE_COL = "age"
    ID_COL = "eid"  
  }else{
    if(dataset == "ADNI"){
      male_table <- filter(adni_male ,  VISCODE.y=="scmri" & DX=="CN" & !is.na(clean_hv_bilateral) )
      female_table <- filter(adni_female ,  VISCODE.y=="scmri" & DX=="CN" & !is.na(clean_hv_bilateral) )
      AGE_COL = "age"
      ID_COL = "PTID"
    }
  }
  
  
  #######################
  # HV PGS CORRELATIONS #
  #######################
  
  # showing the correlation between each of the genetic Scores and HV
  
  CORRS = setNames(data.frame(matrix(ncol = 10, nrow = length(PRS_columns) )), 
                   c("PRS_Threshold", "Male_est", "Male_error" , "Male_Pvalue", 
                     "Female_est", "Female_error" , "Female_Pvalue", 
                     "both_est", "both_error" , "both_Pvalue"))
  CORRS$PRS_Threshold = PRS_columns
  
  for( col in PRS_columns){
    # since I want the regression coefficient to represent the correlations in the range of 0-1, I need to normalize
    # the hv and prs to the range of 0-1
    male_hv <- normalize( male_table$clean_hv_bilateral )
    male_prs <- normalize(male_table[ , col])
    x = summary(lm(male_hv ~ male_prs ))
    male_res = x$coefficients[2,c(1,2,4)]
    
    female_hv <- normalize( female_table$clean_hv_bilateral )
    female_prs <- normalize(female_table[ , col])
    x = summary(lm(female_hv ~ female_prs ))
    female_res = x$coefficients[2,c(1,2,4)]
    
    all_hvs <- c( female_hv , male_hv )
    all_prs <- c( female_prs , male_prs)
    x = summary(lm(all_hvs ~ all_prs ))
    both_res = x$coefficients[2,c(1,2,4)]
    
    CORRS[ CORRS$PRS_Threshold == col , ] =  unlist(c(col ,  male_res , female_res , both_res))
  }
  print(CORRS)
  
  for (sex in c("male", "female")) {
    table = get(paste(sex, "table" , sep = "_"))
    n = nrow(table)
    
    #############################
    # HV SINGLE PGS EXPERIEMNTS #
    #############################
    
    for (col in PRS_columns ) {
      for (seed in SEEDS) {
  
        diff = SUBSAMPLE_DIFF_NOMOGRAM(table , c(col) , AGE_COL , ID_COL , seed)
        
        DIFFS_TABLE <- rbind(DIFFS_TABLE , list(run_ind , exp_ind , grp_ind , "HV", dataset, seed , sex , col , diff , 1))
        run_ind <- run_ind + 1
      }
      exp_ind = exp_ind + 1
    }
    grp_ind <- grp_ind + 1
    
    print(paste( "DONE:" , dataset , sex , "HV" , col))

    #############################
    # HV DOUBLE PGS EXPERIEMNTS #
    #############################
    pairs_1 = combn( PRS_columns , 2 , simplify = FALSE)
    
    for (i in 1:length(pairs_1) ) {
      PRS_1 = pairs_1[[i]][1] 
      PRS_2 = pairs_1[[i]][2] 
      
      table$TEMP = table[, PRS_1 ] + table[, PRS_2 ]

      for (seed in SEEDS) {
        diff = SUBSAMPLE_DIFF_NOMOGRAM(table , c("TEMP") , AGE_COL , ID_COL , seed)
        DIFFS_TABLE <- rbind(DIFFS_TABLE , list( run_ind , exp_ind , grp_ind  ,"HV", dataset, seed , sex , paste(PRS_1, PRS_2, sep = "&") , diff , 2))
        run_ind <- run_ind + 1
      }
      exp_ind = exp_ind + 1
    }
    grp_ind <- grp_ind + 1

    print(paste( "DONE:" , dataset , sex , "HV" , paste(PRS_1, PRS_2, sep = "&")  ))
    
    
    ############################
    # HV LASSO PGS EXPERIEMNTS #
    ############################
    
    # Serious analysis, perform a lasso regression on the 14 columns,
    #                   use the resulting weights to produce a new score,
    #                   use that new score to train models with top/bottom
    #                   scoring samples
    
    for (seed in SEEDS) {
      
      XX <- GPR_ANALYSIS_SEPARATED_LASSO( table ,  PRS_columns , 0.3 , seed = seed , age_column = AGE_COL , id_column = ID_COL, 
                                      XX = seq(min(table[, AGE_COL]) , max(table[, AGE_COL]), 0.1))
      diff <- NOMOGRAM_DIFFERENCE(XX$upper , XX$lower , age_range = range(table[, AGE_COL]))
      
      DIFFS_TABLE <- rbind(DIFFS_TABLE ,list(run_ind , exp_ind , grp_ind  , "HV", dataset, seed , sex , "LASSO" , diff, 3))
      run_ind <- run_ind + 1
      
      LASSO_SCORE( table ,  "clean_hv_bilateral" , PRS_columns , print_weights = TRUE , seed = seed )
    }
    exp_ind = exp_ind + 1
    grp_ind <- grp_ind + 1
    
    print(paste( "DONE:" , dataset , sex , "HV" , "LASSO"))
    
    
  } # closing the male/female table loop
} # closing the dataset loop

DIFFS_TABLE = setNames(DIFFS_TABLE , c("run_index","exp_index","grp_index", "experiment", "dataset" , "seed", "sex", "prs" , "diff" , "prs_num"))
save(DIFFS_TABLE , file = "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/MULTI_PGS_RESULTS_2023_09_12.RData")

load( file = "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/MULTI_PGS_RESULTS_2023_09_12.RData")


#############################################################################################################
## A classic application of PGS's. AD PRS correlation with AD potentially improved by using multiple PGSs  ##
#############################################################################################################


######################
### AD EXPERIMENTS ###
######################

DIFFS_TABLE_NO_AGE = data.frame(matrix(ncol = 10, nrow = 0))

for (dataset in c( "ADNI")){
    
  if(dataset == "UKB"){

    ukb_male$AD_SCORE <- ukb_male$mother_AD + ukb_male$father_AD
    ukb_female$AD_SCORE <- ukb_female$mother_AD + ukb_female$father_AD
    
    
    #ukb_male$AD_SCORE <- as.factor(ukb_male$familial_AD_score)
    #ukb_female$AD_SCORE <- as.factor(ukb_female$familial_AD_score)
    
    
    male_table  =  ukb_male[ !is.na(ukb_male$AD_PRS_TH_1) & !is.na(ukb_male$familial_AD_score) , ]
    female_table = ukb_female[ !is.na(ukb_female$AD_PRS_TH_1) & !is.na(ukb_female$familial_AD_score) , ]
    
    
    AGE_COL = "AGE_Latest"
    ID_COL = "eid"
    AD_SCORE_COL = "AD_SCORE"
    
  }else{
    if(dataset == "ADNI"){
      
      adni_male$AD_SCORE <- as.integer(as.factor(adni_male$DX))
      adni_male$AD_SCORE <- (adni_male$AD_SCORE -1) / 2 
      
      adni_female$AD_SCORE <- as.integer(as.factor(adni_female$DX))
      adni_female$AD_SCORE <- (adni_female$AD_SCORE -1) / 2
      
      male_table  =  adni_male[ !is.na(adni_male$AD_SCORE) , ]
      female_table = adni_female[ !is.na(adni_female$AD_SCORE) , ]
      
      male_table$AD_SCORE_no_AGE = REGRESS_OUT_COLUMN( male_table$AD_SCORE , male_table$AGE)
      female_table$AD_SCORE_no_AGE = REGRESS_OUT_COLUMN( female_table$AD_SCORE , female_table$AGE)
      
      AGE_COL = "AGE"
      ID_COL = "id"
      AD_SCORE_COL = "AD_SCORE_no_AGE"
      
    }else{
      if(dataset == "UCL"){
        male_table  =  ucl_male
        female_table = ucl_female
        AGE_COL = "age"
        ID_COL = "nshdid"
      }
    }
  }

  #######################
  # AD PGS CORRELATIONS #
  #######################
  
  CORRS = setNames(data.frame(matrix(ncol = 10, nrow = length(PRS_AD_columns) )), 
                   c("PRS_Threshold", "Male_est", "Male_error" , "Male_Pvalue", 
                                      "Female_est", "Female_error" , "Female_Pvalue", 
                                      "Both_est", "Both_error" , "Both_Pvalue"))
  CORRS$PRS_Threshold = PRS_AD_columns
  
  for( col in PRS_AD_columns){
    
    # since I want the regression coefficient to represent the correlations in the range of 0-1, I need to normalize
    # the familial AD score and prs to the range of 0-1. 
    male_fad <- normalize( male_table[,AD_SCORE_COL] )
    male_prs <- normalize(male_table[ , col])
    x = summary(lm(male_fad ~ male_prs ))
    male_res = x$coefficients[2,c(1,2,4)]
    
    female_fad <- normalize( female_table[,AD_SCORE_COL] )
    female_prs <- normalize(female_table[ , col])
    x = summary(lm(female_fad ~ female_prs ))
    female_res = x$coefficients[2,c(1,2,4)]
    
    all_fads <- c( female_fad , male_fad )
    all_prs <- c( female_prs , male_prs)
    x = summary(lm(all_fads ~ all_prs ))
    both_res = x$coefficients[2,c(1,2,4)]
    
    CORRS[ CORRS$PRS_Threshold == col , ] =  unlist(c(col ,  male_res , female_res , both_res))
  }
  
  
  for (sex in c("male", "female")) {
    table = get(paste(sex, "table" , sep = "_"))
    
    #############################
    # AD SINGLE PGS EXPERIEMNTS #
    #############################
    
    for (col in PRS_AD_columns) {
      for (seed in SEEDS) {
        
        n = nrow(table)
        set.seed(seed)
        train_rows <- sample(1:n , .5*n )
        test_set <- table[ -train_rows, ]
        
        #test_set[ , AD_SCORE_COL] = factor(test_set[ , AD_SCORE_COL])
        #or = polr( formula( paste0(AD_SCORE_COL," ~ ",AGE_COL ,"+", col)) , data = test_set, Hess = TRUE)
        #diff = or$deviance
        #DIFFS_TABLE <- rbind(DIFFS_TABLE , list(run_ind , exp_ind , grp_ind  , "AD", dataset , seed , sex , col , diff , 1))
        
        test_set[ , AD_SCORE_COL] = factor(test_set[ , AD_SCORE_COL])
        or = polr( formula( paste0(AD_SCORE_COL," ~ ",AGE_COL ,"+", col)) , data = test_set, Hess = TRUE)
        diff = or$deviance
        DIFFS_TABLE_NO_AGE <- rbind(DIFFS_TABLE_NO_AGE , list(run_ind , exp_ind , grp_ind  , "AD", dataset , seed , sex , col , diff , 1))
        
        run_ind <- run_ind + 1
      }
      exp_ind = exp_ind + 1
      print(paste( "DONE:" , dataset , sex , "AD" , col , "SEED:",seed))
      
    }
    grp_ind <- grp_ind + 1
    
    print(paste( "DONE:" , dataset , sex , "AD" , col))
    
    
    ##########################
    # TWO AD PGS EXPERIEMNTS #
    ##########################
    
    for (prs_pairs in combn( PRS_AD_columns , 2 , simplify = FALSE)) {
      prs_1 = prs_pairs[1]
      prs_2 = prs_pairs[2]
      for (seed in SEEDS) {
        
        n = nrow(table)
        set.seed(seed)
        train_rows <- sample(1:n , .5*n )
        test_set <- table[ -train_rows, ]
        
        #test_set[ , AD_SCORE_COL] = factor(test_set[ , AD_SCORE_COL])
        #col = test_set[ , prs_1] + test_set[ , prs_2]
        #or = polr( formula( paste0(AD_SCORE_COL," ~ ",AGE_COL,"+ col")) , data = test_set, Hess = TRUE)
        #diff = or$deviance
        
        #DIFFS_TABLE <- rbind(DIFFS_TABLE , list(run_ind , exp_ind , grp_ind , "AD" , dataset , seed , sex , paste(prs_1 , prs_2,sep = "&") , diff , 2) )
        
        test_set[ , AD_SCORE_COL] = factor(test_set[ , AD_SCORE_COL])
        col = test_set[ , prs_1] + test_set[ , prs_2]
        or = polr( formula( paste0(AD_SCORE_COL," ~ ",AGE_COL,"+ col")) , data = test_set, Hess = TRUE)
        diff = or$deviance
        
        DIFFS_TABLE_NO_AGE <- rbind(DIFFS_TABLE_NO_AGE , list(run_ind , exp_ind , grp_ind , "AD" , dataset , seed , sex , paste(prs_1 , prs_2,sep = "&") , diff , 2) )
        
        run_ind <- run_ind + 1
      }
      exp_ind = exp_ind + 1
    }
    grp_ind <- grp_ind + 1
    
    print(paste( "DONE:" , dataset , sex , "AD" , paste(prs_1 , prs_2,sep = "&") ))
    
    
    ############################
    # AD LASSO PGS EXPERIEMNTS #
    ############################
    for (seed in SEEDS) {
      table$AD_LASSO_SCORE = LASSO_SCORE( table , y_column = table[ , AD_SCORE_COL] , x_columns = PRS_AD_columns , seed = seed, print_weights = FALSE)
      
      n = nrow(table)
      set.seed(seed)
      train_rows <- sample(1:n , .5*n )
      test_set <- table[ -train_rows, ]
      test_set[ , AD_SCORE_COL] = factor(test_set[ , AD_SCORE_COL])
      
      or = polr( formula( paste0(AD_SCORE_COL," ~ ",AGE_COL ,"+ AD_LASSO_SCORE")) , data = test_set, Hess = TRUE)
      diff = or$deviance
      
      #DIFFS_TABLE <- rbind(DIFFS_TABLE , list(run_ind , exp_ind , grp_ind , "AD" , dataset , seed , sex , "LASSO" , diff , 3) )
      

      DIFFS_TABLE_NO_AGE <- rbind(DIFFS_TABLE_NO_AGE , list(run_ind , exp_ind , grp_ind , "AD" , dataset , seed , sex , "LASSO" , diff , 3) )
      
      run_ind <- run_ind + 1
      
    }
    exp_ind = exp_ind + 1
    grp_ind <- grp_ind + 1
    
    print(paste( "DONE:" , dataset , sex , "AD" , "LASSO"))
  }
}

DIFFS_TABLE_NO_AGE = setNames(DIFFS_TABLE , c("run_index","exp_index","grp_index", "experiment", "dataset" , "seed", "sex", "prs" , "diff" , "prs_num"))

#DIFFS_TABLE = setNames(DIFFS_TABLE , c("run_index","exp_index","grp_index", "experiment", "dataset" , "seed", "sex", "prs" , "diff" , "prs_num"))
#save(DIFFS_TABLE , file = "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/MULTI_PGS_RESULTS_2023_09_12.RData")

# OR LOAD PREVIUOSLY COMPUTED TABLE:
load("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/MULTI_PGS_RESULTS_2023_09_12.RData")

# load("/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/MULTI_PGS_RESULTS_NO_AGE_2023_03_17.RData")


#plotting
plt_folder = "/Users/mohammedjanahi/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Paper/Figures/"
 
# Figure 3: Summary of Single vs Double vs Lasso PGS Experiments
p1 = ggplot( data = DIFFS_TABLE[DIFFS_TABLE$experiment == "HV" & DIFFS_TABLE$dataset == "UKB" & DIFFS_TABLE$sex == "female", ] ) +
  geom_boxplot( aes(x = factor(prs_num) , y = diff , group = factor(prs_num) ) , outlier.shape = NA ) + 
  scale_x_discrete(name="PGS Model" ,labels=c("Single","Double","LASSO")) + 
  ylab("Hippocampal Volume Nomogram Separation") + ggtitle("UKBB") +
  theme(plot.title = element_text(hjust = 0.5))
p2 = ggplot( data = DIFFS_TABLE[DIFFS_TABLE$experiment == "AD" & DIFFS_TABLE$dataset == "ADNI" & DIFFS_TABLE$sex == "male", ] ) +
  geom_boxplot( aes(x = factor(prs_num) , y = diff , group = factor(prs_num) ) , outlier.shape = NA ) +
  scale_x_discrete(name="PGS Model" ,labels=c("Single","Double","LASSO")) + 
  ylab("Alzhiemer's Model Residual Deviance") + ggtitle("ADNI") +
  theme(plot.title = element_text(hjust = 0.5))


pdf(file = paste0(plt_folder,"HV_AD_Figure 3 - 1 - test",".pdf") , width = 8 , height = 11)
p1 + p2
dev.off()

pdf(file = paste0(plt_folder,"HV_Figure 3 - 1 - Test",".pdf") , width = 10 , height = 3)
p1  + ylim(20, 160) + coord_flip() 
dev.off()

pdf(file = paste0(plt_folder,"AD_Figure 3 - 1",".pdf") , width = 4 , height = 8)
p2
dev.off()



order_prs = c("1e.08","1e.07","1e.06","1e.05","1e.04","1e.03","0.01","0.05","0.1","0.2","0.4","0.5","0.75","1")

ll = expand.grid(as.list(order_prs), as.list(order_prs))[, c("Var2","Var1")]
ll$Var1 = as.character(ll$Var1)
ll$Var2 = as.character(ll$Var2)
ll = ll[ ll$Var2 != ll$Var1 , ]
ll =  ll %>% filter(!duplicated(paste0(pmax(Var1, Var2), pmin(Var1, Var2))))
ll = paste( ll$Var2 , ll$Var1 , sep = " & ")

DIFFS_TABLE$prs_short = unlist(lapply( str_extract_all( (DIFFS_TABLE[ DIFFS_TABLE$dataset %in% c("UKB") & DIFFS_TABLE$experiment %in% c("HV"), "prs"]), "(?<=TH_).{1,5}(?=_NO)|LASSO" ) , function(e) paste( e , collapse = " & ") ))

DIFFS_TABLE$prs_short = factor(DIFFS_TABLE$prs_short , levels = c(order_prs , ll , "LASSO"))

p = ggplot( DIFFS_TABLE[ DIFFS_TABLE$dataset %in% c("UKB") & DIFFS_TABLE$experiment %in% c("HV"), ] , aes( x = prs_short  , y = diff, fill = prs_short)) +
  geom_boxplot( show.legend = FALSE) + ggtitle("NOMOGRAM SEPARATION ACROSS PGSs IN UKBB") + 
  xlab("PGS Threshold") +
  facet_grid( ~ prs_num , scale = "free_x" , space = "free_x") + 
  theme(axis.text.x = element_text( angle=90))

pdf(file = paste0(plt_folder,"HV_Figure 3 - 3",".pdf") , width = 15 , height = 10)
p 
dev.off()

p = ggplot( DIFFS_TABLE[ DIFFS_TABLE$dataset %in% c("ADNI") & DIFFS_TABLE$experiment %in% c("HV"), ] , aes( x = prs_short  , y = diff, fill = prs_short)) +
  geom_boxplot( show.legend = FALSE) + ggtitle("NOMOGRAM SEPARATION ACROSS PGSs IN ADNI") + 
  xlab("PGS Threshold") +
  facet_grid( ~ prs_num , scale = "free_x" , space = "free_x") + 
  theme(axis.text.x = element_text( angle=90))

pdf(file = paste0(plt_folder,"HV_Figure 3 - 4",".pdf") , width = 15 , height = 10)
p 
dev.off()

# Supplementary figure with all substrata of tests.
pdf(file = paste0(plt_folder,"HV_Figure 3 - 2",".pdf") , width = 6 , height = 8)
par(mfrow=c(2,2))
exp = "HV"
for(dataset in c("UKB","ADNI"))
  for(sex in c("male","female"))
    boxplot(DIFFS_TABLE[ DIFFS_TABLE$experiment == exp & DIFFS_TABLE$dataset == dataset & DIFFS_TABLE$sex == sex, "diff"] 
            ~ DIFFS_TABLE[ DIFFS_TABLE$experiment == exp & DIFFS_TABLE$dataset == dataset & DIFFS_TABLE$sex == sex, "grp_index"] ,
            names = c("Single","Double","Lasso") ,
            main = paste(dataset,sex) ,
            xlab = "PGS Model" ,
            ylab = if(exp=="HV") "Nomogram Separation (mm3)" else "Residual Deviance", 
            outline = FALSE)
dev.off()

pdf(file = paste0(plt_folder,"AD_Figure 3 - 2",".pdf") , width = 6 , height = 8)
par(mfrow=c(2,2))
exp = "AD"
for(dataset in c("ADNI","UKB"))
  for(sex in c("male","female"))
    boxplot(DIFFS_TABLE[ DIFFS_TABLE$experiment == exp & DIFFS_TABLE$dataset == dataset & DIFFS_TABLE$sex == sex, "diff"] 
            ~ DIFFS_TABLE[ DIFFS_TABLE$experiment == exp & DIFFS_TABLE$dataset == dataset & DIFFS_TABLE$sex == sex, "grp_index"] ,
            names = c("Single","Double","Lasso") ,
            main = paste(dataset,sex) ,
            xlab = "PGS Model" ,
            ylab = if(exp=="HV") "Nomogram Separation (mm3)" else "Residual Deviance", 
            outline = FALSE)
dev.off()

# Supplementary figure with ADNI AD experiments with/without age adjusted diagnostic scores
pdf(file = paste0(plt_folder,"AD_Figure 3 - 3",".pdf") , width = 6 , height = 8)
par(mfrow=c(2,2))
exp = "AD"
dataset = "ADNI"
for(sex in c("male","female"))
  for(table in list( DIFFS_TABLE , DIFFS_TABLE_NO_AGE) )
      boxplot(DIFFS_TABLE[ DIFFS_TABLE$experiment == exp & DIFFS_TABLE$dataset == dataset & DIFFS_TABLE$sex == sex, "diff"] 
              ~ DIFFS_TABLE[ DIFFS_TABLE$experiment == exp & DIFFS_TABLE$dataset == dataset & DIFFS_TABLE$sex == sex, "grp_index"] ,
              names = c("Single","Double","Lasso") ,
              main = paste( ifelse( nrow(table) < 10000 , "Age Adjusted" , "" ),sex) ,
              xlab = "PdGS Model" ,
              ylab = "Diagnostic Score Separation", 
              outline = FALSE)
dev.off()


p = ggplot( DIFFS_TABLE[ DIFFS_TABLE$dataset %in% c("UKB") & DIFFS_TABLE$experiment %in% c("HV") & DIFFS_TABLE$prs_num == 1, ] , aes( x = run_ind  , y = diff, fill = prs)) +
      geom_boxplot() + ggtitle("NOMOGRAM SEPARATION ACROSS PGSs IN UKBB") + 
  ylab("nomogram difference (mm3)") + xlab("PGS Threshold")

pdf(file = paste0(plt_folder,"HV_Figure 3 - 3",".pdf") , width = 10 , height = 8)
p
dev.off()

PRS = PRS_columns_sup
PRS_AD = PRS_AD_columns_sup
# basic correlations to hv
ukb_male_hv_cor = cor( ukb_male[,PRS] , ukb_male[,"clean_hv_bilateral"] , use="complete.obs")
ukb_female_hv_cor = cor( ukb_female[,PRS] , ukb_female[,"clean_hv_bilateral"]  , use="complete.obs")
ukb_hv_cor = (ukb_male_hv_cor + ukb_female_hv_cor) / 2

adni_male_hv_cor = cor( adni_male[adni_male$DX=="CN",PRS] , adni_male[adni_male$DX=="CN","clean_hv_bilateral"] , use="complete.obs")
adni_female_hv_cor = cor( adni_female[adni_female$DX=="CN",PRS] , adni_female[adni_female$DX=="CN","clean_hv_bilateral"] , use="complete.obs")
adni_hv_cor = (adni_male_hv_cor + adni_female_hv_cor) / 2

ukb_male_ad_cor = cor( ukb_male[,PRS_AD] , ukb_male[,"familial_AD_score"] , use="complete.obs")
ukb_female_ad_cor = cor( ukb_female[,PRS_AD] , ukb_female[,"familial_AD_score"]  , use="complete.obs")
ukb_ad_cor = (ukb_male_ad_cor + ukb_female_ad_cor) / 2

adni_male_ad_cor = cor( adni_male[,PRS_AD] , adni_male[,"AD_SCORE"] , use="complete.obs")
adni_female_ad_cor = cor( adni_female[,PRS_AD] , adni_female[,"AD_SCORE"] , use="complete.obs")
adni_ad_cor = (adni_male_ad_cor + adni_female_ad_cor) / 2

# correlations of polygenic scores to each other
ukb_male_prs_cor_hv = cor( ukb_male[,PRS] , ukb_male[,PRS] , use = "complete.obs")
ukb_female_prs_cor_hv = cor( ukb_female[,PRS] , ukb_female[,PRS] , use = "complete.obs")
ukb_prs_cor_hv = (ukb_male_prs_cor_hv+ukb_female_prs_cor_hv)/2

adni_male_prs_cor_hv = cor( adni_male[,PRS] , adni_male[,PRS] , use = "complete.obs")
adni_female_prs_cor_hv = cor( adni_female[,PRS] , adni_female[,PRS] , use = "complete.obs")
adni_prs_cor_hv = (adni_male_prs_cor_hv+adni_female_prs_cor_hv)/2

ukb_male_prs_cor_ad = cor( ukb_male[,PRS_AD] , ukb_male[,PRS_AD] , use = "complete.obs")
ukb_female_prs_cor_ad = cor( ukb_female[,PRS_AD] , ukb_female[,PRS_AD] , use = "complete.obs")
ukb_prs_cor_ad = (ukb_male_prs_cor_ad+ukb_female_prs_cor_ad)/2

adni_male_prs_cor_ad = cor( adni_male[,PRS_AD] , adni_male[,PRS_AD] , use = "complete.obs")
adni_female_prs_cor_ad = cor( adni_female[,PRS_AD] , adni_female[,PRS_AD] , use = "complete.obs")
adni_prs_cor_ad = (adni_male_prs_cor_ad+adni_female_prs_cor_ad)/2


ukb_sub_hv = ukb_all[!is.na(ukb_all$PRS_TH_1),]
ukb_sub_ad = ukb_all[!is.na(ukb_all$AD_PRS_TH_1),]

adni_sub_hv = adni_all[!is.na(adni_all$AD_PRS_TH_1),]
adni_sub_ad = adni_all[!is.na(adni_all$AD_PRS_TH_1),]

# supplementary figure showing the correlations between the different pgs high/low subtrata
intersections_ukb_hv = find_intersections( ukb_sub_hv , "eid" , PRS)
intersections_ukb_ad = find_intersections( ukb_sub_ad , "eid" , PRS_AD)
intersections_adni_hv = find_intersections( adni_sub_hv , "id" , PRS)
intersections_adni_ad = find_intersections( adni_sub_ad , "id" , PRS_AD)

# Figure 2: Summary of Polygenic Scores Correlations and Intersections experiments
pl_ukb_hv_cor_both = make_hv_heatmap(ukb_hv_cor , FALSE , 3 , 5 , FALSE , title = "(A) PGS/HV Correlation" ) 
pl_ukb_hv_prs_cor_both = make_cor_heatmap(ukb_prs_cor_hv , FALSE , 2 , 3 , title = "(B) PGS/PGS Correlation")
pl_ukb_hl_hv_sub = make_cor_heatmap( intersections_ukb_hv , leg = FALSE , mid = 0.5 , lims = c(0,1) , two_point = TRUE , title = "(D) High/Low Intersection" ,
                                     xlab = "High                                                                                                Low\nPGS Threshold Bin" ,
                                     ylab = "PGS Threshold Bin\n High                                                                        Low") +
                  geom_hline(yintercept = length(PRS)+0.5) + geom_vline(xintercept = length(PRS)+0.5)
pl_para_hv_ukb = make_paracoord_plot(ukb_sub_hv , "eid" , PRS[1] , PRS , highlight_perc = 0.1 , title = "(C) Order of samples by TH = 10e-7")
layout <- "
ADDDDD
BBCCCC
BBCCCC
"
pdf(file = paste0(plt_folder,"HV_Figure 2 - 1",".pdf") , width = 15 , height = 9)
pl_ukb_hv_cor_both + pl_ukb_hv_prs_cor_both + pl_ukb_hl_hv_sub + pl_para_hv_ukb +
  plot_layout(design = layout) + plot_annotation(title = "Polygenic Scores Correlations and Intersections") & theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Figure 2: Summary of Polygenic Scores Correlations and Intersections experiments
pl_adni_ad_cor_both = make_hv_heatmap(adni_ad_cor  , FALSE , 3 , 5 , FALSE , title = "(A) PGS/AD Correlation" ) 
pl_adni_ad_prs_cor_both = make_cor_heatmap( adni_prs_cor_ad , FALSE , 2 , 3 , title = "(B) PGS/PGS Correlation")
pl_adni_hl_ad_sub = make_cor_heatmap( intersections_adni_ad , leg = FALSE , mid = 0.5 , lims = c(0,1) , two_point = TRUE  , title = "(D) High/Low Intersection" , 
                                      xlab = "High                                                                                                Low\nPGS Threshold Bin" ,
                                      ylab = "PGS Threshold Bin\n High                                                                        Low") +
  geom_hline(yintercept = length(PRS_AD)+0.5) + geom_vline(xintercept = length(PRS_AD)+0.5)
pl_para_ad_adni = make_paracoord_plot(adni_sub_ad , "id" , PRS_AD[1] , PRS_AD , highlight_perc = 0.1 , title = "(C) Order of samples by TH = 10e-7")
layout <- "
ADDDDD
BBCCCC
BBCCCC
"

pdf(file = paste0(plt_folder,"AD_Figure 2 - 1",".pdf") , width = 15 , height = 9)
pl_adni_ad_cor_both + pl_adni_ad_prs_cor_both + pl_adni_hl_ad_sub + pl_para_ad_adni +
  plot_layout(design = layout) + plot_annotation(title = "Polygenic Scores Correlations and Intersections") & theme(plot.title = element_text(hjust = 0.5))
dev.off()



# hv correlations heatmap figures
p_ukb_hv_cor_male = make_hv_heatmap(ukb_male_hv_cor , FALSE , 3 , 5 , FALSE , ylab = "HV PGS Threshold") + ggtitle("\n\nMALE")
p_ukb_hv_cor_female = make_hv_heatmap(ukb_female_hv_cor , FALSE , 3 , 5 , TRUE) + ggtitle("UKB\n\nFEMALE")
p_ukb_hv_cor = make_hv_heatmap(ukb_hv_cor , FALSE , 3 , 5 , TRUE) + ggtitle("\n\nMEAN")

p_adni_hv_cor_male = make_hv_heatmap(adni_male_hv_cor , FALSE , 3 , 5 , TRUE) + ggtitle("\n\nMALE")
p_adni_hv_cor_female = make_hv_heatmap(adni_female_hv_cor , FALSE , 3 , 5 , TRUE) + ggtitle("ADNI\n\nFEMALE")
p_adni_hv_cor = make_hv_heatmap(adni_hv_cor , FALSE , 3 , 5 , TRUE) + ggtitle("\n\nMEAN")

p_ukb_ad_cor_male = make_hv_heatmap(ukb_male_ad_cor , FALSE , 3 , 5 , FALSE , ylab = "AD PRS Threshold") + ggtitle("\n\nMALE")
p_ukb_ad_cor_female = make_hv_heatmap(ukb_female_ad_cor , FALSE , 3 , 5 , TRUE) + ggtitle("UKB\n\nFEMALE")
p_ukb_ad_cor = make_hv_heatmap(ukb_ad_cor , FALSE , 3 , 5 , TRUE) + ggtitle("\n\nMEAN")

p_adni_ad_cor_male = make_hv_heatmap(adni_male_ad_cor , FALSE , 3 , 5 , TRUE) + ggtitle("\n\nMALE")
p_adni_ad_cor_female = make_hv_heatmap(adni_female_ad_cor , FALSE , 3 , 5 , TRUE) + ggtitle("ADNI\n\nFEMALE")
p_adni_ad_cor = make_hv_heatmap(adni_ad_cor , FALSE , 3 , 5 , TRUE) + ggtitle("\n\nMEAN")

pdf(file = paste0(plt_folder,"HV_AD_Figure 2 - 2",".pdf") , width = 15 , height = 12)
p_ukb_hv_cor_male + p_ukb_hv_cor_female + p_ukb_hv_cor + 
  p_adni_hv_cor_male + p_adni_hv_cor_female + p_adni_hv_cor +
  p_ukb_ad_cor_male + p_ukb_ad_cor_female + p_ukb_ad_cor +
  p_adni_ad_cor_male + p_adni_ad_cor_female + p_adni_ad_cor + plot_layout(ncol = 6)
dev.off()

pdf(file = paste0(plt_folder,"HV_Figure 2 - 2",".pdf") , width = 12 , height = 6)
p_ukb_hv_cor_male + p_ukb_hv_cor_female + p_ukb_hv_cor + 
  p_adni_hv_cor_male + p_adni_hv_cor_female + p_adni_hv_cor + plot_layout(ncol = 6)
dev.off()

pdf(file = paste0(plt_folder,"AD_Figure 2 - 2",".pdf") , width = 12 , height = 6)
p_ukb_ad_cor_male + p_ukb_ad_cor_female + p_ukb_ad_cor +
  p_adni_ad_cor_male + p_adni_ad_cor_female + p_adni_ad_cor + plot_layout(ncol = 6)
dev.off()

# prs correlations heatmap figures
p_ukb_hv_prs_cor_male = make_cor_heatmap(ukb_male_prs_cor_hv , FALSE , 2 , 3 ,  rem_x = TRUE , title = "MALE" , big_y_label = "UKB\n\n")
p_ukb_hv_prs_cor_female = make_cor_heatmap(ukb_female_prs_cor_hv , FALSE , 2 , 3 ,  rem_x = TRUE , rem_y = TRUE , title = "HV PGS CORRELATIONS\n\nFEMALE")
p_ukb_hv_prs_cor_both = make_cor_heatmap(ukb_prs_cor_hv , FALSE , 2 , 3 ,  rem_x = TRUE , rem_y = TRUE , title = "MEAN")

p_adni_hv_prs_cor_male = make_cor_heatmap(adni_male_prs_cor_hv , FALSE , 2 , 3 , big_y_label = "ADNI\n\n")
p_adni_hv_prs_cor_female = make_cor_heatmap(adni_female_prs_cor_hv , FALSE , 2 , 3 ,  rem_y = TRUE)
p_adni_hv_prs_cor_both = make_cor_heatmap(adni_prs_cor_hv , FALSE , 2 , 3 , rem_y = TRUE)

pdf(file = paste0(plt_folder,"HV_Figure 2 - 3",".pdf") , width = 20 , height = 12)
p_ukb_hv_prs_cor_male + p_ukb_hv_prs_cor_female + p_ukb_hv_prs_cor_both +
  p_adni_hv_prs_cor_male + p_adni_hv_prs_cor_female + p_adni_hv_prs_cor_both 
dev.off()

p_ukb_ad_prs_cor_male = make_cor_heatmap(ukb_male_prs_cor_ad , FALSE , 2 , 3 ,  rem_x = TRUE , title = "MALE" , big_y_label = "UKB\n\n")
p_ukb_ad_prs_cor_female = make_cor_heatmap(ukb_female_prs_cor_ad , FALSE , 2 , 3 ,  rem_x = TRUE , rem_y = TRUE , title = "AD PGS CORRELATIONS\n\nFEMALE")
p_ukb_ad_prs_cor_both = make_cor_heatmap(ukb_prs_cor_ad , FALSE , 2 , 3 ,  rem_x = TRUE , rem_y = TRUE , title = "MEAN")

p_adni_ad_prs_cor_male = make_cor_heatmap(adni_male_prs_cor_ad , FALSE , 2 , 3 , big_y_label = "ADNI\n\n")
p_adni_ad_prs_cor_female = make_cor_heatmap(adni_female_prs_cor_ad , FALSE , 2 , 3 ,  rem_y = TRUE)
p_adni_ad_prs_cor_both = make_cor_heatmap(adni_prs_cor_ad , FALSE , 2 , 3 , rem_y = TRUE)

pdf(file = paste0(plt_folder,"AD_Figure 2 - 3",".pdf") , width = 20 , height = 12)
p_ukb_ad_prs_cor_male + p_ukb_ad_prs_cor_female + p_ukb_ad_prs_cor_both +
  p_adni_ad_prs_cor_male + p_adni_ad_prs_cor_female + p_adni_ad_prs_cor_both 
dev.off()

# prs high/low correlations heatmap figures
p_ukb_hl_hv_sub = make_cor_heatmap( intersections_ukb_hv , leg = FALSE , mid = 0.5 , lims = c(0,1) , two_point = TRUE , title = "UKB") + geom_hline(yintercept = length(PRS)+0.5) + geom_vline(xintercept = length(PRS)+0.5)
p_ukb_hl_ad_sub = make_cor_heatmap( intersections_ukb_ad , leg = FALSE , mid = 0.5 , lims = c(0,1) , two_point = TRUE , title = "UKB") + geom_hline(yintercept = length(PRS)+0.5) + geom_vline(xintercept = length(PRS)+0.5)
p_adni_hl_hv_sub = make_cor_heatmap( intersections_adni_hv , leg = FALSE , mid = 0.5 , lims = c(0,1) , two_point = TRUE , title = "ADNI") + geom_hline(yintercept = length(PRS)+0.5) + geom_vline(xintercept = length(PRS)+0.5)
p_adni_hl_ad_sub = make_cor_heatmap( intersections_adni_ad , leg = FALSE , mid = 0.5 , lims = c(0,1) , two_point = TRUE , title = "ADNI") + geom_hline(yintercept = length(PRS)+0.5) + geom_vline(xintercept = length(PRS)+0.5)

pdf(file = paste0(plt_folder,"HV_Figure 2 - 4",".pdf") , width = 25 , height = 10)
p_ukb_hl_hv_sub +  p_adni_hl_hv_sub + 
  plot_annotation(title = "HIGH/LOW HV PGS INTERSECTION") & theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(file = paste0(plt_folder,"AD_Figure 2 - 4",".pdf") , width = 25 , height = 10)
p_adni_hl_ad_sub + p_ukb_hl_ad_sub + 
  plot_annotation(title = "HIGH/LOW AD PGS INTERSECTION") & theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# prs paracord plots
p_para_hv_ukb = make_paracoord_plot(ukb_sub_hv , "eid" , "PRS_TH_1e.08" , PRS , highlight_perc = 0.1 , title = "UKB")
p_para_ad_ukb = make_paracoord_plot(ukb_sub_ad , "eid" , "AD_PRS_TH_1e.08" , PRS_AD , highlight_perc = 0.1 , title = "UKB")
p_para_hv_adni = make_paracoord_plot(adni_sub_hv , "id" , "PRS_TH_1e.08" , PRS , highlight_perc = 0.1 , title = "ADNI")
p_para_ad_adni = make_paracoord_plot(adni_sub_ad , "id" , "AD_PRS_TH_1e.08" , PRS_AD , highlight_perc = 0.1 , title = "ADNI")

pdf(file = paste0(plt_folder,"HV_Figure 2 - 5",".pdf") , width = 15 , height = 5)
p_para_hv_ukb + p_para_hv_adni +
  plot_annotation(title = "HV PGS Order of samples by TH = 10e-7 ") & theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = paste0(plt_folder,"AD_Figure 2 - 5",".pdf") , width = 15 , height = 5)
p_para_ad_adni + p_para_ad_ukb +
  plot_annotation(title = "AD PGS Order of samples by TH = 10e-7 ") & theme(plot.title = element_text(hjust = 0.5))
dev.off()

for( prs in PRS){
  assign( paste0("p_para_hv_ukb_",prs) , 
          make_paracoord_plot(ukb_sub_hv , "eid" , prs , PRS , highlight_perc = 0.1 , title = paste0("sorted by ",prs))  )
  
  assign( paste0("p_para_ad_adni_",prs) , 
          make_paracoord_plot(adni_sub_ad , "id" , paste0("AD_",prs) , PRS_AD , highlight_perc = 0.1 , title = paste0("sorted by ",prs))  )
}


#p_para_hv_ukb_PRS_TH_1e.08 + p_para_hv_ukb_PRS_TH_1e.07 + p_para_hv_ukb_PRS_TH_1e.06 + p_para_hv_ukb_PRS_TH_1e.05 +
#  p_para_hv_ukb_PRS_TH_1e.04 + p_para_hv_ukb_PRS_TH_1e.03 + p_para_hv_ukb_PRS_TH_0.01 + p_para_hv_ukb_PRS_TH_0.05 +
#  p_para_hv_ukb_PRS_TH_0.1 + p_para_hv_ukb_PRS_TH_0.2 + p_para_hv_ukb_PRS_TH_0.4 + p_para_hv_ukb_PRS_TH_0.5 +
#  p_para_hv_ukb_PRS_TH_0.75 + p_para_hv_ukb_PRS_TH_1 + plot_layout(ncol = 2) + plot_annotation(title = "UKB HV PGS") & theme(plot.title = element_text(hjust = 0.5))

pdf(file = paste0(plt_folder,"HV_Figure 2 - 6",".pdf") , width = 15 , height = 15)
p_para_hv_ukb_PRS_TH_1e.08 + p_para_hv_ukb_PRS_TH_1e.05 + p_para_hv_ukb_PRS_TH_0.01 + 
  p_para_hv_ukb_PRS_TH_0.1 + p_para_hv_ukb_PRS_TH_0.4 + p_para_hv_ukb_PRS_TH_1 + 
  plot_layout(ncol = 2) + plot_annotation(title = "UKB HV PGS") & theme(plot.title = element_text(hjust = 0.5))
dev.off()


pdf(file = paste0(plt_folder,"AD_Figure 2 - 6",".pdf") , width = 15 , height = 15)
p_para_ad_adni_PRS_TH_1e.08 + p_para_ad_adni_PRS_TH_1e.05 + p_para_ad_adni_PRS_TH_0.01 + 
  p_para_ad_adni_PRS_TH_0.1 + p_para_ad_adni_PRS_TH_0.4 + p_para_ad_adni_PRS_TH_1 + 
  plot_layout(ncol = 2) + plot_annotation(title = "ADNI AD PGS") & theme(plot.title = element_text(hjust = 0.5))
dev.off()


# Get upper triangle of the correlation matrix
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)]<- NA
  return(cormat)
}

make_cor_heatmap = function(table , leg=TRUE ,d = 3 , size = 4 , mid = 0.5 , lims = c(0,1) , rem_y = FALSE , rem_x = FALSE , title = NA , big_y_label = "" , two_point=TRUE , ylab = "PGS Threshold" , xlab = "PGS Threshold"){

  p = ggplot(data = melt(get_lower_tri(round(table , digits = d)) , na.rm = TRUE), aes(x=Var1, y=Var2, fill=value)) +
    geom_tile(color="white") +
    geom_text(aes(Var1, Var2, label = value), color = "black", size = size) + 
    scale_x_discrete(name = xlab , labels = gsub( "AD_|_NO_GPC|PRS_TH_|_high|_low" , "" , row.names(table) ) )   +
    scale_y_discrete(name = paste0(big_y_label,ylab) , labels = gsub( "AD_|_NO_GPC|PRS_TH_|_high|_low" , "" , row.names(table) ) )   
  
  if(two_point)
    p = p + scale_fill_gradient(low = "white", high = "red", limit = lims, space = "Lab", name="Pearson\nCorrelation" , guide = ifelse(leg,"colourbar","none"))
  else
    p = p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = mid, limit = lims, space = "Lab", name="Pearson\nCorrelation" , guide = ifelse(leg,"colourbar","none"))
  
  if(rem_y)
    p = p + theme(axis.text.y = element_blank()  , axis.ticks.y = element_blank() , axis.title.y = element_blank() )
  
  if(rem_x)
    p = p + theme(axis.text.x = element_blank()  , axis.ticks.x = element_blank() , axis.title.x = element_blank() )
  
  if(!is.na(title))
    p = p + labs(title = title) + theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

make_hv_heatmap = function(table , leg=TRUE ,d = 3 , size = 4 , rem_y = FALSE , title = NA , ylab = "PGS Threshold"){
  
  p =  ggplot(data.frame(x = factor(row.names(table) , levels = row.names(table)), y = table), aes(x = x, y = 1, fill = y)) +
    geom_tile(color="white" , stat="identity") +
    scale_fill_gradient2(low = "white", high = "red", midpoint = (min(table) + max(table))/2, limit = c( min(table) ,max(table)), space = "Lab", name="Pearson\nCorrelation" , guide = ifelse(leg,"colourbar","none")) +
    geom_text(aes(x, 1, label = round(y,d)), color = "black", size = size) + 
    scale_x_discrete(name = ylab , labels = gsub("AD_","",gsub( "_NO_GPC" , "" , gsub( "PRS_TH_" , "" , row.names(table) ))) )  +
    theme(axis.title.x=element_blank()  , axis.text.x = element_blank()  , axis.ticks.x = element_blank() ) +
    coord_flip() +
    theme(plot.title = element_text(hjust = 0.5))
  
  if(rem_y)
    p = p + theme(axis.text.y = element_blank()  , axis.ticks.y = element_blank() , axis.title.y = element_blank() )
  
  if(!is.na(title))
    p = p + labs(title = title) + theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

make_paracoord_plot = function(table , id_col , sort_col = "PRS_TH_1" , PRS_columns , highlight_perc = 0.1 , title = ""){
  
  table = table[order(table[ , sort_col]) , ]
  
  data = data.frame(matrix( table[,id_col] , nrow = nrow(table) , ncol = length(PRS_columns)+3 ))
  for(i in 3:(length(PRS_columns)+2))
    data[ , i] = order(table[,PRS_columns[i-2]] , decreasing = TRUE)
  
  perc =  highlight_perc * nrow(data)
  data[,1] = 1
  data[ (1:perc), 1] = 2
  data[ ((nrow(data)-perc):nrow(data)) , 1] = 3
  
  data[,1] = as.factor(data[,1])
  
  data[,2] = data[,3]
  data[,ncol(data)] = data[,ncol(data)-1]
  
  names(data) = c("categories" , "  ", PRS_columns , " ")
  
  data = data[order(data[,1]) , ]
  
  p = ggparcoord(data[seq(1,nrow(data) , 50) , ], columns = 2:(length(PRS_columns)+3), groupColumn = 1, alpha = 0 , scale="globalminmax") +
    scale_color_manual(values=c( "grey", "blue", "red") , guide="none") +
    geom_line(size=1.5 , alpha=.3) +
    geom_hline(yintercept=perc , alpha=0.4) +
    geom_hline(yintercept=nrow(data)- perc , alpha=0.4) +
    labs( x = "PGS Threshold" , y  ="RANK" , title = title) +
    scale_x_discrete(labels = gsub("AD_","",gsub( "_NO_GPC_ZSCORE" , "" , gsub( "PRS_TH_" , "" , names(data)[2:ncol(data)] ))) ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

find_intersections = function( table , id_col , prs_cols){

  table = data.frame(table %>% distinct( table[ , id_col] ,.keep_all = TRUE))
  
  data = data.frame(matrix( table[,id_col] , nrow = nrow(table) , ncol = length(prs_cols) ))
  
  for(i in 1:length(prs_cols) )
    data[ , i] = table[order(table[,prs_cols[i]]) , id_col] 
  names(data) = prs_cols
  
  L = floor(.3*nrow(data))
  data_lists = list()
  for(i  in 1:length(prs_cols) ){
    data_lists[[ (2*i)-1 ]] = ( data[ 1:L ,i] )
    data_lists[[ 2*i ]] = ( data[ (nrow(data)-L):nrow(data) ,i] )
    
  }
  n = length(data_lists)
  intersections = matrix(mapply(function(u,v) length(intersect(u,v)), rep(data_lists, n), rep(data_lists, each=n)), ncol=n)
  
  rownames(intersections) = c(rbind( paste0(prs_cols,"_high") , paste0(prs_cols,"_low") ))
  colnames(intersections) = c(rbind( paste0(prs_cols,"_high") , paste0(prs_cols,"_low") ))
  
  intersections = intersections[ c( paste0(prs_cols,"_high") , paste0(prs_cols,"_low")) , c( paste0(prs_cols,"_high") , paste0(prs_cols,"_low")) ] / (L+1)
  
  return( round(intersections , digits = 2))
}

FIND_DIFFS <- function( table , cols , seeds){
  
  AD_SINGLE_DIFFS = data.frame(matrix(ncol = 4, nrow = 0))
  iter_ind = 1 
  
  n = nrow(table)
  
  for (col in cols ){
    
    all_diffs = c()
    for(seed in seeds){
      set.seed(seed)
      train_rows <- sample(1:n, .5*n)
      test_set <- table[-train_rows , ] 
      test_set = test_set[ order(test_set[,col]) , ]
      table_high <- tail( test_set , .3*n )
      table_low <- head( test_set , .3*n )
      
      diff <- (mean(table_high$familial_AD_score,na.rm=TRUE)) - (mean(table_low$familial_AD_score,na.rm=TRUE))
      
      AD_SINGLE_DIFFS <- rbind(AD_SINGLE_DIFFS , list(iter_ind , sex , col , diff) )
      
    }
    iter_ind = iter_ind + 1
  }
  
  AD_SINGLE_DIFFS = setNames( AD_SINGLE_DIFFS , c("index" , "sex", "prs" , "diff") )
  
  return(AD_SINGLE_DIFFS)
}


DIFFS_TABLE_2 = data.frame(matrix(ncol = 16, nrow = 0))
run_ind = 1
exp_ind = 1 
grp_ind = 1

d_s_imp = c()
l_s_imp = c()
l_d_imp = c()
#plotting
for(exp in c("HV","AD")){
  for(dataset in c("UKB","ADNI")){
    for(sex in c("male","female")){
      for(seed in SEEDS){
        diffs = DIFFS_TABLE[ DIFFS_TABLE$experiment == exp & DIFFS_TABLE$dataset == dataset & DIFFS_TABLE$sex == sex & DIFFS_TABLE$seed==seed,]
        
        
        diffs_single = diffs[diffs$prs_num == 1 , "diff"]
        diffs_double = diffs[diffs$prs_num == 2 , "diff"]
        diffs_lasso = diffs[diffs$prs_num == 3  , "diff"]
        
        #double_improvement = 0
        #for( prs_col in diffs[ diffs$prs_num == 1 , "prs"] ){
        # for every prs, find the index of that prs and the index of all double prs's that include that one prs,
        # the paired result is then the difference between the first element of this list, versus the rest of the list. 
        #  diffs_to_use = diffs[grepl(diffs$prs[1] , diffs$prs) , "diff"]
        #  diff_one = diffs_to_use[1]
        #  diffs_two = mean(diffs_to_use[-1])
        
        #  double_improvement = double_improvement + diffs_two / diff_one
        #}
        
        #double_improvement = double_improvement / length(diffs[ diffs$prs_num == 1 , "prs"])
        
        
        double_improvement = mean(diffs_double , na.rm=TRUE) /  mean(diffs_single , na.rm=TRUE)
        lasso_improvement = mean(diffs_lasso , na.rm=TRUE) /  mean(diffs_single , na.rm=TRUE)
        lasso_d_improvement = mean(diffs_lasso , na.rm=TRUE) /  mean(diffs_double , na.rm=TRUE)
        
        d_s_imp = c(d_s_imp , double_improvement)
        l_s_imp = c(l_s_imp , lasso_improvement)
        l_d_imp = c(l_d_imp , lasso_d_improvement)
        
        DIFFS_TABLE_2 = rbind( DIFFS_TABLE_2 , list(run_ind , exp_ind , grp_ind , exp, dataset , sex , seed ,
                                                    mean( diffs_single , na.rm = TRUE) , sd(diffs_single , na.rm = TRUE), 
                                                    mean( diffs_double , na.rm = TRUE) , sd(diffs_double , na.rm = TRUE), 
                                                    mean( diffs_lasso , na.rm = TRUE) , sd(diffs_lasso , na.rm = TRUE),
                                                    double_improvement , lasso_improvement , lasso_d_improvement ) )
        
        run_ind = run_ind + 1
      }
    }
    exp_ind = exp_ind + 1
  }
  run_ind = run_ind + 1
}
DIFFS_TABLE_2 = setNames(DIFFS_TABLE_2 , c("run_index","exp_index","grp_index", "experiment", "dataset", "sex", "seed" ,
                                           "single_prs_diff_mean" , "single_prs_diff_sd" ,
                                           "double_prs_diff_mean" , "double_prs_diff_sd" , 
                                           "lasso_prs_diff_mean" , "lasso_prs_diff_sd" , 
                                           "single_double_diff" , "single_lasso_diff" , "double_lasso_diff"))


SUMMARY_DIFFS = data.frame(matrix(ncol = 6, nrow = 0))
for(exp in c("HV","AD")){
  for(dataset in c("UKB","ADNI")){
    for(sex in c("male","female")){
      diffs = DIFFS_TABLE_2[ DIFFS_TABLE_2$experiment == exp & DIFFS_TABLE_2$dataset == dataset & DIFFS_TABLE_2$sex == sex, ]
      mean_s_d_imp = round(mean(diffs[,"single_double_diff"] ) , digits = 3)
      mean_s_l_imp = round(mean(diffs[,"single_lasso_diff"] ) , digits = 3)
      mean_d_l_imp = round(mean(diffs[,"double_lasso_diff"] ) , digits = 3)
      SUMMARY_DIFFS = rbind( SUMMARY_DIFFS , list(exp , dataset , sex , mean_s_d_imp , mean_s_l_imp , mean_d_l_imp))
    }
  }
}
names(SUMMARY_DIFFS) = c("EXPERIMENT" , "DATASET" , "SEX" , "Single to Double improvment" , "Single to Lasso improvment" , "Double to Lasso improvment")

print(paste("UKB:",colMeans(SUMMARY_DIFFS[c(1,2) , c(4,5,6)])))
print(paste("ADNI:",colMeans(SUMMARY_DIFFS[c(3,4) , c(4,5,6)])))

print(SUMMARY_DIFFS[ c(1,2), c(1,2,3,4,6,5)])





LASSO_WEIGHTS = data.frame(matrix(ncol = length(prs_columns) + 1 , nrow = 0))


for( dataset in c( "ADNI") ){
  
  if(dataset == "UKB"){
    male_table  =  ukb_male[!is.na(ukb_male$AD_PRS_TH_1) & !is.na(ukb_male$clean_hv_bilateral) , ]
    female_table = ukb_female[!is.na(ukb_female$AD_PRS_TH_1) & !is.na(ukb_female$clean_hv_bilateral) , ]

  }else{
    if(dataset == "ADNI"){
      male_table <- filter(adni_male ,  VISCODE=="bl" & DX=="CN" & !is.na(clean_hv_bilateral) )
      female_table <- filter(adni_female ,  VISCODE=="bl" & DX=="CN" & !is.na(clean_hv_bilateral) )

    }
  }
  
  for (sex in c("male", "female")) {
    table = get(paste(sex, "table" , sep = "_"))
    
    for (seed in SEEDS) {
      lasso = LASSO_SCORE( table ,  table$clean_hv_bilateral , PRS_columns , print_weights = TRUE , seed = seed )
      #LASSO_WEIGHTS <- rbind(LASSO_WEIGHTS , c( paste0(dataset , "_" , sex) , seed , lasso$weights))
      print(lasso$weights)
    }
    
    
  } 
} 
