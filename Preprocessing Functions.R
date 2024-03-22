# File Containing all pre-processing functions that deal with Nomograms
# By: Mohammed Janahi
#
#



PRERPOCESS_ALL <- function( ukb_table , adni_table , ucl_table , epad_table , GIF = FALSE){
  
  # call the function to filter table, take care of outliers, and correct hippocampal volume for ICV and scan date
  x <- PREPROCESS_UKBB(ukb_table , GIF = GIF)
  
  # grab the separate gender tables
  ukb_male <- x$ukb_img_ex_outliers_male
  ukb_female <- x$ukb_img_ex_outliers_female
  ukb_male_ex <- x$excluded_males
  ukb_female_ex <- x$excluded_females
  
  ukb_all <- rbind(ukb_male,ukb_female)
  ukb_ex_all <- rbind(ukb_male_ex,ukb_female_ex)
  
  names(ukb_all)[names(ukb_all)=="ICV.x"] = "ICV"
  names(ukb_ex_all)[names(ukb_ex_all)=="ICV.x"] = "ICV"
  
  
  # ADNI Table is read from downloaded library
  x <- PREPROCESS_ADNI(adni_table  ,include_neuro_tests = TRUE)
  adni_male <- x$ADNI_male
  adni_female <- x$ADNI_female
  adni_all <- rbind(adni_male , adni_female)
  
  
  # UCL 46' Data is saved as CSV file
  x <- PREPROCESS_I46( ucl_table )
  ucl_male <- x$ucl_male
  ucl_female <- x$ucl_female
  ucl_all <- rbind(ucl_male , ucl_female)
  
  # EPADS Preprocessing
  # read in the file that has most of the pre-processing done. 
  x <- PREPROCESS_EPAD( epad_table )
  epad_male <- x$epad_male
  epad_female <- x$epad_female
  epad_all <- rbind(epad_male , epad_female)
  
  return(list( ukb_all = ukb_all , ukb_ex = ukb_ex_all , adni_all = adni_all , ucl_all = ucl_all , epad_all = epad_all))
  
}


HARMONIZE_ALL <- function(ukb_all , adni_all , ucl_all , epad_all){
  
  prs_columns <- grep("^INTERSECT_PRS_TH" , names(ucl_all) , value = TRUE)
  prs_columns_nc <- paste0(prs_columns,"_nc")
  prs_columns_nc_gpc <- paste0(prs_columns,"_nc_gpc")

  # Redefine the tables to not include NAs in the relavent columns
  ukb_all_nc = ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all[,prs_columns[1]]) & !is.na(ukb_all$ICV) , ]
  #adni_all_nc = adni_all[!is.na(adni_all$clean_hv_bilateral) & !is.na(adni_all$DX) & (adni_all$VISCODE.x == "bl") & !is.na(adni_all[,prs_columns[1]])& !is.na(adni_all$ICV) , ]
  adni_all_nc = adni_all[!is.na(adni_all$clean_hv_bilateral) & !is.na(adni_all$DX) & !is.na(adni_all[,prs_columns[1]])& !is.na(adni_all$ICV) , ]
  ucl_all_nc <- filter(ucl_all , !is.na(ucl_all$clean_hv_bilateral) & !is.na(ucl_all[,prs_columns[1]])  & !is.na(ucl_all$ICV) )
  epad_all_nc <- filter(epad_all , !is.na(epad_all$clean_hv_bilateral) & !is.na(epad_all[,prs_columns[1]])  & !is.na(epad_all$ICV)  )
  
  
  #constructing the covariate columns
  ukb_ICVs <- ukb_all_nc$ICV
  adni_ICVs <- adni_all_nc$ICV
  ucl_ICVs <- ucl_all_nc$ICV
  epad_ICVs <- epad_all_nc$ICV
  ICVs <- c(ukb_ICVs , adni_ICVs , ucl_ICVs , epad_ICVs)
  
  
  ukb_sex <- as.logical(ukb_all_nc$sex)
  adni_sex <- adni_all_nc$sex
  ucl_sex <- as.logical(ucl_all_nc$sex)
  epad_sex <- epad_all_nc$sex
  sexs <- as.factor(c(ukb_sex , adni_sex , ucl_sex , epad_sex))
  
  ages <- c(ukb_all_nc$age , adni_all_nc$age , ucl_all_nc$age , epad_all_nc$age )
  
  ukb_DX <- rep("1" , nrow(ukb_all_nc))
  adni_DX <- adni_all_nc$DX
  ucl_DX <- rep("1" , nrow(ucl_all_nc))
  epad_DX <- rep("1" , nrow(epad_all_nc))
  
  diagnosis <- as.double(as.factor(c(ukb_DX , adni_DX , ucl_DX  ,epad_DX)))
  
  
  ukb_hvs = ukb_all_nc$clean_hv_bilateral
  adni_hvs = adni_all_nc$clean_hv_bilateral
  ucl_hvs = ucl_all_nc$clean_hv_bilateral
  epad_hvs = epad_all_nc$clean_hv_bilateral
  hvs = c(ukb_hvs , adni_hvs , ucl_hvs , epad_hvs)
  
  snp_wieghts_gpcs = c("snp_weights_gpc_EA_1" , "snp_weights_gpc_EA_2" , "snp_weights_gpc_NA_1", "snp_weights_gpc_NA_2", "snp_weights_gpc_NA_3")
  for(i in 1:5)
    assign( paste0("gpc_",i) , c(ukb_all_nc[,snp_wieghts_gpcs[i]] ,
                                 adni_all_nc[,snp_wieghts_gpcs[i]] , 
                                 ucl_all_nc[,snp_wieghts_gpcs[i]] , 
                                 epad_all_nc[,snp_wieghts_gpcs[i]] ) )
  
  
  ukb_indices  = 1:nrow(ukb_all_nc)
  adni_indices = (nrow(ukb_all_nc)+1):(nrow(ukb_all_nc) + nrow(adni_all_nc))
  ucl_indices  = (nrow(ukb_all_nc) + nrow(adni_all_nc) + 1):(nrow(ukb_all_nc) + nrow(adni_all_nc) + nrow(ucl_all_nc))
  epad_indices = (nrow(ukb_all_nc) + nrow(adni_all_nc) + nrow(ucl_all_nc) + 1):(length(hvs))
  
  ## HARMONIZE PGSs ##
  b = c( rep(0 , nrow(ukb_all_nc)) , rep(1 , nrow(adni_all_nc)) , rep(2 , nrow(ucl_all_nc)), rep(3 , nrow(epad_all_nc)) )
  
  
  d = t(as.matrix(rbind(ukb_all_nc[ , prs_columns] , adni_all_nc[ , prs_columns] , ucl_all_nc[ , prs_columns] , epad_all_nc[ , prs_columns] )))
  covs <- model.matrix(~ gpc_1 + gpc_2 + gpc_3 + gpc_4 + gpc_5)
  nc = neuroCombat(dat=d , batch = b , ref.batch = 0 , mod = covs)
  
  i= 1
  for( col in prs_columns_nc){
    ukb_all_nc[ , col ] = (nc$dat.combat[i, ukb_indices ])
    adni_all_nc[, col ] = (nc$dat.combat[i, adni_indices ])
    ucl_all_nc[ , col ] = (nc$dat.combat[i, ucl_indices ])
    epad_all_nc[, col ] = (nc$dat.combat[i, epad_indices ])
    
    assign( paste0("COV_",col) , nc$dat.combat[i,] )
    
    i = i + 1
  }
  
  ## HARMONIZE ICVS ##
  
  d = t(as.matrix(rbind(ukb_all_nc[, c("ICV" ,"clean_hv_bilateral")] , 
                        adni_all_nc[,c("ICV" ,"clean_hv_bilateral")] , 
                        ucl_all_nc[, c("ICV" ,"clean_hv_bilateral")] , 
                        epad_all_nc[,c("ICV" ,"clean_hv_bilateral")] )))
  
  
  # Keep only age sex and gpc
  covs <- model.matrix( as.formula(paste0("~ sexs + ages + ", paste(paste0("gpc_",1:5) , collapse = " + ") )) )
  nc = neuroCombat(dat=d , batch = b , ref.batch = 0 , mod = covs , )
  
  ukb_all_nc[ , "ICV_nc" ] = (nc$dat.combat[1, ukb_indices ])
  adni_all_nc[, "ICV_nc" ] = (nc$dat.combat[1, adni_indices ])
  ucl_all_nc[, "ICV_nc" ] = (nc$dat.combat[1, ucl_indices ])
  epad_all_nc[, "ICV_nc" ] = (nc$dat.combat[1, epad_indices ])
  
  ICVs_nc <- nc$dat.combat[1,]
  
  
  ## HARMONIZE HVs ##
  d = t(as.matrix(rbind(ukb_all_nc[ , c("clean_hv_left" ,"clean_hv_right","clean_hv_bilateral")] , 
                        adni_all_nc[ ,c("clean_hv_left" ,"clean_hv_right","clean_hv_bilateral")] , 
                        ucl_all_nc[ , c("clean_hv_left" ,"clean_hv_right","clean_hv_bilateral")] ,
                        epad_all_nc[ ,c("clean_hv_left" ,"clean_hv_right","clean_hv_bilateral")])))
  
  
  # age sex icv_nc and prs_nc
  covs <- model.matrix( as.formula(paste0("~ICVs_nc + sexs + ages + " , paste( paste0("COV_",prs_columns_nc) , collapse = " + ") )) )
  nc = neuroCombat(dat=d , batch = b , ref.batch = 0 , mod = covs)
  
  ukb_all_nc[ , "clean_hv_bilateral_nc" ] = (nc$dat.combat[3, ukb_indices ])
  adni_all_nc[, "clean_hv_bilateral_nc" ] = (nc$dat.combat[3, adni_indices ])
  ucl_all_nc[, "clean_hv_bilateral_nc" ] = (nc$dat.combat[3, ucl_indices ])
  epad_all_nc[, "clean_hv_bilateral_nc" ] = (nc$dat.combat[3, epad_indices ])
  
  # now regress the GPCs out of the PRSs , the ICV's, and the HV's
  for(col in c(prs_columns_nc,"ICV_nc","clean_hv_bilateral_nc")){
    
    # x =  REGRESS_OUT_COLUMNS( rbind(rbind(rbind( ukb_all_nc[ , c(col , snp_wieghts_gpcs) ] , 
    #                                              adni_all_nc[ , c(col , snp_wieghts_gpcs) ]), 
    #                                       ucl_all_nc[ , c(col , snp_wieghts_gpcs) ]), 
    #                                 epad_all_nc[ , c(col , snp_wieghts_gpcs) ])
    #                           , col , snp_wieghts_gpcs)
    # 
    # ukb_all_nc[ ,paste0(col,"_gpc") ] = x[  1:nrow(ukb_all_nc) ]
    # adni_all_nc[,paste0(col,"_gpc") ] = x[ (1:nrow(adni_all_nc)) + (nrow(ukb_all_nc)) ]
    # ucl_all_nc[ ,paste0(col,"_gpc") ] = x[ (1:nrow(ucl_all_nc)) + (nrow(ukb_all_nc) + nrow(adni_all_nc)) ]
    # epad_all_nc[,paste0(col,"_gpc") ] = x[ (1:nrow(epad_all_nc)) + (nrow(ukb_all_nc)+nrow(adni_all_nc)+nrow(ucl_all_nc)) ]
    
    
    ukb_all_nc[ ,paste0(col,"_gpc") ] = REGRESS_OUT_COLUMNS( ukb_all_nc[ , c(col , snp_wieghts_gpcs) ] , col , snp_wieghts_gpcs )
    adni_all_nc[,paste0(col,"_gpc") ] = REGRESS_OUT_COLUMNS( adni_all_nc[ , c(col , snp_wieghts_gpcs) ] , col , snp_wieghts_gpcs )
    ucl_all_nc[ ,paste0(col,"_gpc") ] = REGRESS_OUT_COLUMNS( ucl_all_nc[ , c(col , snp_wieghts_gpcs) ] , col , snp_wieghts_gpcs )
    epad_all_nc[,paste0(col,"_gpc") ] = REGRESS_OUT_COLUMNS( epad_all_nc[ , c(col , snp_wieghts_gpcs) ] , col , snp_wieghts_gpcs )
  }
  
  # regress the ICV out of the HV's
  col = "clean_hv_bilateral_nc_gpc"
  
  # x =  REGRESS_OUT_COLUMN( c(ukb_all_nc[ ,col] , adni_all_nc[ ,col]  , ucl_all_nc[ ,col] , epad_all_nc[ ,col] ) ,
  #                          c(ukb_all_nc$ICV_nc_gpc , adni_all_nc$ICV_nc_gpc , ucl_all_nc$ICV_nc_gpc , epad_all_nc$ICV_nc_gpc))
  # 
  # ukb_all_nc[ ,paste0(col,"_icv") ] = x[  1:nrow(ukb_all_nc) ]
  # adni_all_nc[,paste0(col,"_icv") ] = x[ (1:nrow(adni_all_nc)) + (nrow(ukb_all_nc)) ]
  # ucl_all_nc[ ,paste0(col,"_icv") ] = x[ (1:nrow(ucl_all_nc)) + (nrow(ukb_all_nc) + nrow(adni_all_nc)) ]
  # epad_all_nc[,paste0(col,"_icv") ] = x[ (1:nrow(epad_all_nc)) + (nrow(ukb_all_nc)+nrow(adni_all_nc)+nrow(ucl_all_nc)) ] 
  # 
  
  ukb_all_nc[ ,paste0(col,"_icv") ] = REGRESS_OUT_COLUMN( ukb_all_nc[ ,col] , ukb_all_nc$ICV_nc_gpc )
  adni_all_nc[,paste0(col,"_icv") ] = REGRESS_OUT_COLUMN( adni_all_nc[ ,col] , adni_all_nc$ICV_nc_gpc )
  ucl_all_nc[ ,paste0(col,"_icv") ] = REGRESS_OUT_COLUMN( ucl_all_nc[ ,col] , ucl_all_nc$ICV_nc_gpc )
  epad_all_nc[,paste0(col,"_icv") ] = REGRESS_OUT_COLUMN( epad_all_nc[ ,col] , epad_all_nc$ICV_nc_gpc )
  
  return(list(ukb_all_nc = ukb_all_nc , adni_all_nc = adni_all_nc ,
              ucl_all_nc = ucl_all_nc , epad_all_nc = epad_all_nc))
  
  }

PREPROCESS_UKBB <- function( ukb , GIF = FALSE){
  
  #ukb <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_features_20200529.csv")
  #ukb_ext <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_features_20200602_FSA.csv")
  # Merge the two tables
  #ukb <- merge(ukb , ukb_ext , by=intersect(names(ukb),names(ukb_ext) )) # merge by intersect makes sure we don't duplicate columns
  
  # For readability, rename some columns
  names(ukb)[names(ukb) == 'X31.0.0']    <- 'sex'
  names(ukb)[names(ukb) == 'X21003.0.0'] <- 'Age.When.Attended.Assesment.Center'
  names(ukb)[names(ukb) == 'X21003.2.0'] <- 'Age.Imaging.Visit'
  names(ukb)[names(ukb) == 'X25000.2.0'] <- 'Volumetric.scaling.from.T1.head.image.to.standard.space'
  names(ukb)[names(ukb) == 'X53.2.0']    <- 'Scan.Date'
  names(ukb)[names(ukb) == 'X25019.2.0'] <- 'HV.Left.FSL'
  names(ukb)[names(ukb) == 'X25020.2.0'] <- 'HV.Right.FSL'
  names(ukb)[names(ukb) == 'X26562.2.0'] <- 'HV.Left'
  names(ukb)[names(ukb) == 'X26593.2.0'] <- 'HV.Right'
  names(ukb)[names(ukb) == 'X21000.0.0'] <- 'Ethnic.Background'
  # re-nameing columns with multiple instances by using regexp to grab all instances and attach the instance number to the new column name
  names(ukb)[ grepl( "^X22009.0." , names(ukb)) ] = paste( "Genetic.PC." , substring( grep( "^X22009.0." , names(ukb) , value = TRUE) , first = 10) , sep="")
  #names(ukb)[ grepl( "^X20002.2." , names(ukb)) ] = paste( "Self.Reported.Condition." , substring( grep( "^X20002.2." , names(ukb) , value = TRUE)  , first = 10) , sep="")
  
  
  # First step is to select only the samples that have imaging data
  ukb_img <- ukb[!is.na(ukb$Age.Imaging.Visit),]
  ukb_img <- ukb_img[!is.na(ukb_img$HV.Left),]
  ukb_img <- ukb_img[!is.na(ukb_img$HV.Right),]
  
  # Ages in default ukbb table are just years without months/decimals.
  # but taking exact birth date and scan date into account will give more accurate ages.
  ukb_ages <- read.csv("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/ukb_all_subjs.csv")
  # create a column with the age from the latest imaging visit (some subjects came for multiple visits)
  ukb_ages$age <- ukb_ages$AGE2
  ukb_ages$age [is.na(ukb_ages$AGE2)] <- ukb_ages$AGE[is.na(ukb_ages$AGE2)]
  
  ukb_img <- merge(ukb_img , ukb_ages[,c("eid","age")] , by="eid" , all.x = TRUE) 
  
  # ICVs that better match the adni and i46 ICVs
  ukb_icv <- read.csv("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/fsetiv.csv")
  ukb_icv$ICV <- ukb_icv[,2] * (1/ukb_icv[,4])
  ukb_img <- merge(ukb_img , ukb_icv[ , c("eid","ICV")] , by="eid" , all.x = TRUE) 
  
  # file containing more diagnosis and brain volume columns 
  more_columns <- read.csv("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/ukb_imaging_dx.csv")  
  ukb_img = merge( ukb_img , more_columns[ , c(-3,-4,-5) ] , by="eid" , all.x = TRUE)
  
  gif_table <- PREPROCESS_GIF_FILE("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/with_imaging_40k/ParcellationsDatabase_UKB_GIF.csv")  
  ukb_img <- merge(ukb_img , gif_table , by="eid" , all.x = TRUE) 
  
  # Include a score for family history of AD
  ukb_img = MERGE_FAMILY_HISTORY_AD(ukb_img)
  
  # Genetics
  ukb_img = MERGE_GENETICS(ukb_img)
  
  # Genetic Principle Components
  ukb_gpc_1 = read.table("~/OneDrive - University College London/ANDRE STAGE/Datasets/UKBB_SNPweights/EA.predpc")
  ukb_gpc_2 = read.table("~/OneDrive - University College London/ANDRE STAGE/Datasets/UKBB_SNPweights/NA.predpc")
  names(ukb_gpc_1) = c("eid" , "label" , "num_snps" , "snp_weights_gpc_EA_1" , "snp_weights_gpc_EA_2" , "prob_NW" , "prob_SE" , "prob_AJ")
  names(ukb_gpc_2) = c("eid" , "label" , "num_snps" , "snp_weights_gpc_NA_1" , "snp_weights_gpc_NA_2" , "snp_weights_gpc_NA_3" , "prob_1" , "prob_2" , "prob_3" , "prob_4")
  
  ukb_img = merge(ukb_img , ukb_gpc_1[,c("eid","snp_weights_gpc_EA_1" , "snp_weights_gpc_EA_2")] , by="eid")
  ukb_img = merge(ukb_img , ukb_gpc_2[,c("eid","snp_weights_gpc_NA_1" , "snp_weights_gpc_NA_2" , "snp_weights_gpc_NA_3")] , by="eid")
  
  
  # Exclusions
  
  # Need to exclude a few sample EID's (patient choose to be excluded)
  private_eids <- c('1090030','1091961','1230112','1325363','1503907','1546281','1779627','2072543','2129742','2328077','2381848',
                    '2562577','2682290','2850777','2991997','3279656','3303963','3629688','3640543','3819653','3850604','3874773',
                    '3983605','3999481','4217643','4278871','4362486','4588727','4652426','4684825','4838503','5122856','5203595',
                    '5494205','5708826','5912074','5954147')
  
  ukb_img_ex <- ukb_img[! ukb_img$eid %in% private_eids , ]
  
  # Need to include only White/British ethnicity.
  # value 1001 means "British". So we are only including those for now.
  ukb_img_ex <- ukb_img_ex[which(ukb_img_ex$Ethnic.Background==1001),]
  
  # Need to exclude people with Neurological or Psychiatric Disorders
  # and people with substance abuse or history of head trauma
  # and people with cardiovascular disorders
  # 1079  cardiomyopathy
  # 1240  neurological injury/trauma
  # 1243  psychological/psychiatric problem
  # 1258  chronic/degenerative neurological problem
  # 1261	multiple sclerosis
  # 1262	parkinsons disease
  # 1263	dementia/alzheimers/cognitive impairment
  # 1264	epilepsy
  # 1266	head injury
  # 1289	schizophrenia
  # 1408  alcohol dependency
  # 1409  opioid dependency
  # 1410  other substance abuse/dependency
  # 1425	cerebral aneurysm
  # 1434  other neurological problem
  # 1491	brain hemorrhage
  exclude_conditions <- c( 1079 , 1240 , 1243 , 1258 , 1261 , 1262 , 1263 , 1264 , 1266 , 1289 , 
                           1408 , 1409 , 1410 , 1425 , 1434 , 1491 )
  #excluded_rows <- ukb_img_ex
  
  ukb_img_ex$exclude = 0
  
  ukb_img_ex$exclude <- WHICH_TO_EXCLUDE(ukb_img_ex)
  
  # Filter the self reported conditions
  #ukb_img_ex[ WHICH_CONTAINS( ukb_img_ex , grep("^X20002" , names(ukb_img_ex) , value = TRUE) , exclude_conditions) , "exclude"] = 1
  
  
  #excluded_rows <- subset( excluded_rows , !( excluded_rows$eid %in% ukb_img_ex$eid) )
  
  # stratify by gender
  ukb_img_ex_outliers_male   <- ukb_img_ex[ which(ukb_img_ex$sex == 1),]
  ukb_img_ex_outliers_female <- ukb_img_ex[ which(ukb_img_ex$sex == 0),]
  
  
  # Outliers
  
  # next step is to filter out outliers.
  # This is done by filtering out all volumes that are more than 5 MAE's away from the mean. 
  # Find the mean of both HV
  ukb_img_ex_outliers_male <- MAD_FILTER(ukb_img_ex_outliers_male , c("HV.Left","HV.Right") , 5)
  ukb_img_ex_outliers_female <- MAD_FILTER(ukb_img_ex_outliers_female , c("HV.Left","HV.Right") , 5)
  
  ukb_img_ex_outliers_male <- MAD_FILTER(ukb_img_ex_outliers_male , c("HV.Left.GIF","HV.Right.GIF") , 5)
  ukb_img_ex_outliers_female <- MAD_FILTER(ukb_img_ex_outliers_female , c("HV.Left.GIF","HV.Right.GIF") , 5)
  
  
  
  # Confounders
  
  # next step is to correct for some variables by regressing them out. 
  
  # do the regressing out of variables independently for each gender
  if(GIF){
    ukb_img_ex_outliers_male$HV.Bilateral <- ( ukb_img_ex_outliers_male$HV.Left.GIF + ukb_img_ex_outliers_male$HV.Right.GIF )/2
    ukb_img_ex_outliers_female$HV.Bilateral <- ( ukb_img_ex_outliers_female$HV.Left.GIF + ukb_img_ex_outliers_female$HV.Right.GIF )/2
    ukb_img_ex_outliers_male$clean_hv_left <- ukb_img_ex_outliers_male$HV.Left.GIF
    ukb_img_ex_outliers_male$clean_hv_right <- ukb_img_ex_outliers_male$HV.Right.GIF
    ukb_img_ex_outliers_female$clean_hv_left <- ukb_img_ex_outliers_female$HV.Left.GIF
    ukb_img_ex_outliers_female$clean_hv_right <- ukb_img_ex_outliers_female$HV.Right.GIF
  }else{
    ukb_img_ex_outliers_male$HV.Bilateral <- ( ukb_img_ex_outliers_male$HV.Left + ukb_img_ex_outliers_male$HV.Right )/2
    ukb_img_ex_outliers_female$HV.Bilateral <- ( ukb_img_ex_outliers_female$HV.Left + ukb_img_ex_outliers_female$HV.Right )/2
    ukb_img_ex_outliers_male$clean_hv_left <- ukb_img_ex_outliers_male$HV.Left
    ukb_img_ex_outliers_male$clean_hv_right <- ukb_img_ex_outliers_male$HV.Right
    ukb_img_ex_outliers_female$clean_hv_left <- ukb_img_ex_outliers_female$HV.Left
    ukb_img_ex_outliers_female$clean_hv_right <- ukb_img_ex_outliers_female$HV.Right
  }

  ukb_img_ex_outliers_male$clean_hv_bilateral <- ( ukb_img_ex_outliers_male$clean_hv_left + ukb_img_ex_outliers_male$clean_hv_right )/2
  ukb_img_ex_outliers_female$clean_hv_bilateral <- (ukb_img_ex_outliers_female$clean_hv_left + ukb_img_ex_outliers_female$clean_hv_right)/2

  #ukb_img_ex_outliers_male$clean_hv_left <- REGRESS_OUT("HV.Left" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_male)
  #ukb_img_ex_outliers_male$clean_hv_right <- REGRESS_OUT("HV.Right" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_male)
  #ukb_img_ex_outliers_male$clean_hv_bilateral <- ( ukb_img_ex_outliers_male$clean_hv_left + ukb_img_ex_outliers_male$clean_hv_right )/2
  
  #ukb_img_ex_outliers_female$clean_hv_left <- REGRESS_OUT("HV.Left" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_female)
  #ukb_img_ex_outliers_female$clean_hv_right <- REGRESS_OUT("HV.Right" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_female)
  #ukb_img_ex_outliers_female$clean_hv_bilateral <- (ukb_img_ex_outliers_female$clean_hv_left + ukb_img_ex_outliers_female$clean_hv_right)/2
  
  # do the regressing out of variables independently for each gender 
  #ukb_img_ex_outliers_male$clean_hv_left_gif <- REGRESS_OUT("HV.Left.GIF" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_male)
  #ukb_img_ex_outliers_male$clean_hv_right_gif <- REGRESS_OUT("HV.Right.GIF" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_male)
  #ukb_img_ex_outliers_male$clean_hv_bilateral_gif <- ( ukb_img_ex_outliers_male$clean_hv_left_gif + ukb_img_ex_outliers_male$clean_hv_right_gif )/2
  
  #ukb_img_ex_outliers_female$clean_hv_left_gif <- REGRESS_OUT("HV.Left.GIF" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_female)
  #ukb_img_ex_outliers_female$clean_hv_right_gif <- REGRESS_OUT("HV.Right.GIF" , "Volumetric.scaling.from.T1.head.image.to.standard.space", ukb_img_ex_outliers_female)
  #ukb_img_ex_outliers_female$clean_hv_bilateral_gif <- (ukb_img_ex_outliers_female$clean_hv_left_gif + ukb_img_ex_outliers_female$clean_hv_right_gif)/2
  
  
  ukb_male = ukb_img_ex_outliers_male[ukb_img_ex_outliers_male$exclude == 0 ,]
  excluded_males = ukb_img_ex_outliers_male[ukb_img_ex_outliers_male$exclude == 1 ,]
  ukb_female = ukb_img_ex_outliers_female[ukb_img_ex_outliers_female$exclude == 0 ,]
  excluded_females = ukb_img_ex_outliers_female[ukb_img_ex_outliers_female$exclude == 1 ,]
  
  return( list(ukb_img_ex_outliers_male = ukb_male,
               ukb_img_ex_outliers_female = ukb_female,
               excluded_males = excluded_males , excluded_females = excluded_females ) )
}

PREPROCESS_UKBB_NO_HV <- function( ukb ){
  
  #ukb <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_features_20200529.csv")
  #ukb_ext <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_features_20200602_FSA.csv")
  # Merge the two tables
  #ukb <- merge(ukb , ukb_ext , by=intersect(names(ukb),names(ukb_ext) )) # merge by intersect makes sure we don't duplicate columns
  
  # For readability, rename some columns
  names(ukb)[names(ukb) == 'X31.0.0']    <- 'sex'
  names(ukb)[names(ukb) == 'X21003.0.0'] <- 'Age.When.Attended.Assesment.Center'
  names(ukb)[names(ukb) == 'X21003.2.0'] <- 'Age.Imaging.Visit'
  names(ukb)[names(ukb) == 'X25000.2.0'] <- 'Volumetric.scaling.from.T1.head.image.to.standard.space'
  names(ukb)[names(ukb) == 'X53.2.0']    <- 'Scan.Date'
  names(ukb)[names(ukb) == 'X21000.0.0'] <- 'Ethnic.Background'
  # re-nameing columns with multiple instances by using regexp to grab all instances and attach the instance number to the new column name
  names(ukb)[ grepl( "^X22009.0." , names(ukb)) ] = paste( "Genetic.PC." , substring( grep( "^X22009.0." , names(ukb) , value = TRUE) , first = 10) , sep="")
  #names(ukb)[ grepl( "^X20002.2." , names(ukb)) ] = paste( "Self.Reported.Condition." , substring( grep( "^X20002.2." , names(ukb) , value = TRUE)  , first = 10) , sep="")
  
  
  # file containing more diagnosis and brain volume columns 
  more_columns <- read.csv("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/ukb_imaging_dx.csv")  
  ukb = merge( ukb , more_columns , by="eid" , all.x = TRUE)

  # Include a score for family history of AD
  ukb = MERGE_FAMILY_HISTORY_AD(ukb)
  
  # Genetics
  ukb = MERGE_GENETICS(ukb)
  
  
  # Exclusions
  
  # Need to exclude a few sample EID's (patient choose to be excluded)
  private_eids <- c('1090030','1091961','1230112','1325363','1503907','1546281','1779627','2072543','2129742','2328077','2381848',
                    '2562577','2682290','2850777','2991997','3279656','3303963','3629688','3640543','3819653','3850604','3874773',
                    '3983605','3999481','4217643','4278871','4362486','4588727','4652426','4684825','4838503','5122856','5203595',
                    '5494205','5708826','5912074','5954147')
  
  ukb_ex <- ukb[! ukb$eid %in% private_eids , ]
  
  # Need to include only White/British ethnicity.
  # value 1001 means "British". So we are only including those for now.
  ukb_ex <- ukb_ex[which(ukb_ex$Ethnic.Background==1001),]
  
  # Need to exclude people with Neurological or Psychiatric Disorders
  # and people with substance abuse or history of head trauma
  # and people with cardiovascular disorders
  # 1079  cardiomyopathy
  # 1240  neurological injury/trauma
  # 1243  psychological/psychiatric problem
  # 1258  chronic/degenerative neurological problem
  # 1261	multiple sclerosis
  # 1262	parkinsons disease
  # 1263	dementia/alzheimers/cognitive impairment
  # 1264	epilepsy
  # 1266	head injury
  # 1289	schizophrenia
  # 1408  alcohol dependency
  # 1409  opioid dependency
  # 1410  other substance abuse/dependency
  # 1425	cerebral aneurysm
  # 1434  other neurological problem
  # 1491	brain hemorrhage
  exclude_conditions <- c( 1079 , 1240 , 1243 , 1258 , 1261 , 1262 , 1263 , 1264 , 1266 , 1289 , 
                           1408 , 1409 , 1410 , 1425 , 1434 , 1491 )
  #excluded_rows <- ukb_img_ex
  
  ukb_ex$exclude = 0
  
  ukb_ex$exclude <- WHICH_TO_EXCLUDE(ukb_ex)
  
  # Filter the self reported conditions
  #ukb_img_ex[ WHICH_CONTAINS( ukb_img_ex , grep("^X20002" , names(ukb_img_ex) , value = TRUE) , exclude_conditions) , "exclude"] = 1
  
  
  #excluded_rows <- subset( excluded_rows , !( excluded_rows$eid %in% ukb_img_ex$eid) )
  
  # stratify by gender
  ukb_ex_outliers_male   <- ukb_ex[ which(ukb_ex$sex == 1),]
  ukb_ex_outliers_female <- ukb_ex[ which(ukb_ex$sex == 0),]
  
  ukb_male = ukb_ex_outliers_male[ukb_ex_outliers_male$exclude == 0 ,]
  excluded_males = ukb_ex_outliers_male[ukb_ex_outliers_male$exclude == 1 ,]
  ukb_female = ukb_ex_outliers_female[ukb_ex_outliers_female$exclude == 0 ,]
  excluded_females = ukb_ex_outliers_female[ukb_ex_outliers_female$exclude == 1 ,]
  
  return( list(ukb_ex_outliers_male = ukb_male,
               ukb_ex_outliers_female = ukb_female,
               excluded_males = excluded_males , excluded_females = excluded_females ) )
}

PREPROCESS_ADNI <- function( ADNI_table , include_neuro_tests = FALSE){
  
  # rename some columns for readability
  # col names : 
  # ST88SV : Volume (WM Parcellation) of RightHippocampus
  # ST29SV : Volume (WM Parcellation) of LeftHippocampus
  names(ADNI_table)[names(ADNI_table) == 'ST88SV'] <- 'HV_Right'
  names(ADNI_table)[names(ADNI_table) == 'ST29SV'] <- 'HV_Left'
  ADNI_table$sex <- (ADNI_table$PTGENDER == "Male")
  
  # From duplicated columns, Keep the .x's and drop the .y's
  names(ADNI_table)[grep( "\\.x" , names(ADNI_table))] = str_replace(  grep( "\\.x" , names(ADNI_table),value = TRUE),".x","")
  ADNI_table = ADNI_table[ , !(colnames(ADNI_table) %in% grep( "\\.y" , names(ADNI_table),value = TRUE) )]
  
  ADNI_Gentic_PCs <- read.csv("~/Desktop/UKB NOMOGRAM PROJECT/adni_data/adni_genetic_pcs.csv" , header = TRUE)
  
  ADNI_table <- merge(ADNI_table , ADNI_Gentic_PCs , by = "RID" )
  
  adni_gpc_1 = read.table("~/OneDrive - University College London/ANDRE STAGE/Datasets/ADNI_SNPweights/EA.predpc")
  adni_gpc_2 = read.table("~/OneDrive - University College London/ANDRE STAGE/Datasets/ADNI_SNPweights/NA.predpc")
  
  names(adni_gpc_1) = c("PTID" , "label" , "num_snps" , "snp_weights_gpc_EA_1" , "snp_weights_gpc_EA_2" , "prob_NW" , "prob_SE" , "prob_AJ")
  names(adni_gpc_2) = c("PTID" , "label" , "num_snps" , "snp_weights_gpc_NA_1" , "snp_weights_gpc_NA_2" , "snp_weights_gpc_NA_3" , "prob_WA" , "prob_EU" , "prob_EA" , "prob_NA")
  
  ADNI_table = merge(ADNI_table , adni_gpc_1[ , c(-2,-3)] , by="PTID")
  ADNI_table = merge(ADNI_table , adni_gpc_2[ , c(-2,-3)] , by="PTID")
  
  #Rscript ../ukb_data/PRSice_mac/PRSice.R --prsice ./../ukb_data/PRSice_mac/PRSice \
  #          --base ../ukb_data/CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL \
  #          --target merged_imputed_maf001_geno010_4batches \
  #           --extract ../ukb_data/with_imaging_40k/ukb_cal_merged_maf01.rsids \
  #          --missing SET_ZERO \
  #          --geno 0.1 --maf 0.05 --no-regress --fastscore --all-score --score sum \
  #          --binary-target F --thread 1 \
  #          --beta --stat Beta --A1 Allele1 --snp RSNUMBERS --pvalue P.value \
  #          --bar-levels 0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.05,0.1,0.2,0.4,0.5,0.75,1.0 \
  #          --out PRS_HV_FREESURFER_ADNI
  
  #ADNI_PRS <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/adni_data/PRS_HV_FREESURFER_ADNI.all.score")
  #ADNI_PRS <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/adni_data/PRS_HV_FREESURFER_ADNI3.all.score")
  ADNI_HV_PRS <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/adni_data/PRS_HV_FREESURFER_ADNI_NO_EXTRACT.all.score")
  ADNI_table <- merge(ADNI_table , ADNI_HV_PRS , by.x="PTID" , by.y="IID" , all.x=TRUE)
  
  ADNI_AD_PRS <- READ_PRS_TABLE("~/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data/ADNI_AD_PRS.all.score" , "AD")
  ADNI_table <- merge(ADNI_table , ADNI_AD_PRS , by.x="PTID" , by.y="IID" , all.x=TRUE)
  
  
  ADNI_UKB_PRS <- READ_PRS_TABLE("~/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data/PRS_HV_INTERSECT_WITH_UKBB.all.score" , "UKB_INTERSECT")
  ADNI_table <- merge(ADNI_table , ADNI_UKB_PRS , by.x="PTID" , by.y="IID" , all.x=TRUE)
  

  ADNI_UKB_I46_PRS <- READ_PRS_TABLE("~/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data/PRS_HV_INTERSECT_WITH_UKB_I46.all.score" , "UKB_I46_INTERSECT")
  ADNI_table <- merge(ADNI_table , ADNI_UKB_I46_PRS , by.x="PTID" , by.y="IID" , all.x=TRUE)
  
  ADNI_UKB_I46_EPAD_PRS <- READ_PRS_TABLE("~/OneDrive - University College London/ANDRE STAGE/Datasets/adni_data/PRS_HV_INTERSECT_WITH_UKB_I46_EPAD.all.score" , "INTERSECT")
  ADNI_table <- merge(ADNI_table , ADNI_UKB_I46_EPAD_PRS , by.x="PTID" , by.y="IID" , all.x=TRUE)
  
  # We need to clean up the data before we can use it
  # exclude non-british/white samples 
  ADNI_filter <- ADNI_table[ ADNI_table$PTETHCAT == "Not Hisp/Latino" , ]
  ADNI_filter <- ADNI_filter[ ADNI_filter$PTRACCAT == "White" , ]
  
  ADNI_filter[ ADNI_filter$RHIPQC == "Fail" , "HV_Right" ] = NA
  ADNI_filter[ ADNI_filter$LHIPQC == "Fail" , "HV_Left" ] = NA
  
  ADNI_filter$clean_hv_right <- ADNI_filter$HV_Right  
  ADNI_filter$clean_hv_left <- ADNI_filter$HV_Left  
  ADNI_filter$clean_hv_bilateral <- ( ADNI_filter$clean_hv_left + ADNI_filter$clean_hv_right ) / 2
  
  if(include_neuro_tests){
    
  adni_id_ref = read.csv("Downloads/ROSTER_12Oct2023.csv")
  adni_CDR = read.csv("Downloads/Neuropsychological (1)/CDR_18Sep2023.csv")
  adni_MMSE = read.csv("Downloads/Neuropsychological (1)/MMSE_18Sep2023.csv")
  adni_ADAS = read.csv("Downloads/Neuropsychological (1)/ADASSCORES_18Sep2023.csv")
  adni_GDS = read.csv("Downloads/Neuropsychological (1)/GDSCALE_18Sep2023.csv")
  adni_MOCA = read.csv("Downloads/Neuropsychological (1)/MOCA_18Sep2023.csv")
  adni_MODH = read.csv("Downloads/Neuropsychological (1)/MODHACH_18Sep2023.csv")
  adni_ECOG = read.csv("Downloads/Neuropsychological/ECOGSP_12Jul2023.csv")
  
  # add PTID to ECOG file
  adni_ECOG = merge(adni_ECOG , adni_id_ref[,c("RID","PTID")] , by = "RID" , all.x = 1)
  
  names(ADNI_filter)[ names(ADNI_filter) == "VISCODE"] = "VISCODE2"
  names(adni_ADAS)[ names(adni_ADAS) == "VISCODE"] = "VISCODE2"
  
  adni_CDR$CDRSB = rowSums(adni_CDR[,c("CDMEMORY","CDORIENT","CDJUDGE","CDCOMMUN","CDHOME","CDCARE")])
  adni_CDR[ !is.na(adni_CDR$CDRSB) & (adni_CDR$CDRSB < 0) , "CDRSB"] = NA 
  names(adni_CDR)[ names(adni_CDR) == "CDRSB"] = "CDRSB_ref"

  names(adni_MOCA)[ names(adni_MOCA) == "MOCA"] = "MOCA_ref"
  names(adni_MMSE)[ names(adni_MMSE) == "MMSCORE"] = "MMSE_ref"
  names(adni_ECOG)[ grep("Ecog" , names(adni_ECOG))] = paste0(grep("Ecog" , names(adni_ECOG) , value = TRUE),"_ref")
  
  names(adni_CDR)[ names(adni_CDR) == "VISDATE"] = "EXAM_DATE_CDR"
  names(adni_MMSE)[ names(adni_MMSE) == "VISDATE"] = "EXAM_DATE_MMSE"
  names(adni_ADAS)[ names(adni_ADAS) == "EXAMDATE"] = "EXAM_DATE_ADAS"
  names(adni_GDS)[ names(adni_GDS) == "VISDATE"] = "EXAM_DATE_GDS"
  names(adni_MOCA)[ names(adni_MOCA) == "VISDATE"] = "EXAM_DATE_MOCA"
  names(adni_MODH)[ names(adni_MODH) == "VISDATE"] = "EXAM_DATE_MODH"
  
  # filter some rows which we cannot place in time. 
  adni_CDR = adni_CDR[!(adni_CDR$EXAM_DATE_CDR %in% "" & adni_CDR$VISCODE2 %in% "") , ]
  adni_MMSE = adni_MMSE[!(adni_MMSE$EXAM_DATE_MMSE %in% "" & adni_MMSE$VISCODE2 %in% "") , ]
  adni_GDS = adni_GDS[!(adni_GDS$EXAM_DATE_GDS %in% "" & adni_GDS$VISCODE2 %in% "") , ]
  
  # filter out any rows where we have no imaging for the subjects
  adni_CDR = adni_CDR[ adni_CDR$PTID %in% ADNI_filter$PTID , ]
  adni_MMSE = adni_MMSE[ adni_MMSE$PTID %in% ADNI_filter$PTID , ]
  adni_ADAS = adni_ADAS[ adni_ADAS$PTID %in% ADNI_filter$PTID , ]
  adni_GDS = adni_GDS[ adni_GDS$PTID %in% ADNI_filter$PTID , ]
  adni_MOCA = adni_MOCA[ adni_MOCA$PTID %in% ADNI_filter$PTID , ]
  adni_MODH = adni_MODH[ adni_MODH$PTID %in% ADNI_filter$PTID , ]
  adni_ECOG = adni_ECOG[ adni_ECOG$PTID %in% ADNI_filter$PTID , ]
  
  # merge all sub-databses with the full adni dataset
  ADNI_filter = merge(ADNI_filter , adni_CDR[ , c("PTID","VISCODE2" , grep("CD" , names(adni_CDR) , value = TRUE) ) ] , by = c("PTID","VISCODE2") , all=TRUE)
  ADNI_filter = merge(ADNI_filter , adni_ADAS[ , c("PTID","VISCODE2", "EXAM_DATE_ADAS" ,"TOTAL11","TOTALMOD" , grep("Q" , names(adni_ADAS) , value = TRUE)) ] , by = c("PTID","VISCODE2") , all=TRUE)
  ADNI_filter = merge(ADNI_filter , adni_MMSE[ , c("PTID","VISCODE2" , grep("MM" , names(adni_MMSE) , value = TRUE)) ] , by = c("PTID","VISCODE2") , all=TRUE)
  ADNI_filter = merge(ADNI_filter , adni_MOCA[ , c(4,5,8,11:54) ] , by = c("PTID","VISCODE2") , all=TRUE)
  ADNI_filter = merge(ADNI_filter , adni_GDS[ , c("PTID","VISCODE2" , grep("GD" , names(adni_GDS) , value = TRUE)) ] , by = c("PTID","VISCODE2") , all=TRUE)
  ADNI_filter = merge(ADNI_filter , adni_MODH[ , c("PTID","VISCODE2", "EXAM_DATE_MODH"  , grep("HM" , names(adni_MODH) , value = TRUE)) ] , by = c("PTID","VISCODE2") , all=TRUE)
  ADNI_filter = merge(ADNI_filter , adni_ECOG[ , c(3,6,10:ncol(adni_ECOG) )] , by = c("PTID","VISCODE2") , all=TRUE)
  
  
  ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAMDATE"] = ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAM_DATE_CDR" ] 
  ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAMDATE"] = ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAM_DATE_ADAS"]
  ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAMDATE"] = ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAM_DATE_MMSE"]
  ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAMDATE"] = ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAM_DATE_GDS"]
  ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAMDATE"] = ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAM_DATE_MOCA"]
  ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAMDATE"] = ADNI_filter[ is.na(ADNI_filter$EXAMDATE) , "EXAM_DATE_MODH"]
  
  # all baseline and characteristic columns should remain the same, so make sure they are filled in with the same values for every participant
  colunms_that_should_be_the_same = c("RID" , "sex" , "PTGENDER" , "PTEDUCAT" , "PTETHCAT" , 
                                      "PTRACCAT" , "PTMARRY" , "APOE4" , "SITE" , 
                                      "COLPROT" , "ORIGPROT" , "AGE" , 
                                      grep( "snp_weig" , names(ADNI_filter) , value = TRUE) ,
                                      grep( "INTERSECT" , names(ADNI_filter) , value = TRUE) , 
                                      grep( "ICV" , names(ADNI_filter) , value = TRUE) ,
                                      grep( ".bl$" , names(ADNI_filter) , value = TRUE) ,
                                      grep( "prob_" , names(ADNI_filter) , value = TRUE))

  ADNI_filter = as.data.frame(ADNI_filter %>% group_by(PTID) %>% mutate(across(colunms_that_should_be_the_same, ~ { if (is.Date(.)) {
      coalesce(., na.omit(.)[1])
    } else {
      ifelse(!is.na(.), ., first(na.omit(.)))
    }})))
  
  # fill in the columns from adni merge with the newly found values from the separate datasets. 
  ADNI_filter[ is.na(ADNI_filter$ADAS11) , "ADAS11"] = ADNI_filter[ is.na(ADNI_filter$ADAS11) , "TOTAL11"]
  ADNI_filter[ is.na(ADNI_filter$ADAS13) , "ADAS13"] = ADNI_filter[ is.na(ADNI_filter$ADAS13) , "TOTALMOD"]
  
  fill_cols = c("MMSE","CDRSB","MOCA","EcogSPDivatt","EcogSPLang","EcogSPMem","EcogSPOrgan","EcogSPPlan","EcogSPVisspat","EcogSPTotal")
  
  for( col in fill_cols){
    ADNI_filter[ is.na(ADNI_filter[,col]) & (ADNI_filter$VISCODE2 %in% c("bl") | ADNI_filter$VISCODE2 %in% c("sc")) , col] =
    ADNI_filter[ is.na(ADNI_filter[,col]) & (ADNI_filter$VISCODE2 %in% c("bl") | ADNI_filter$VISCODE2 %in% c("sc")) , paste0(col,".bl")]
    ADNI_filter[ is.na(ADNI_filter[,col]), col] = ADNI_filter[ is.na(ADNI_filter[,col]) , paste0(col,"_ref")]
    
  }
  
  # we can now recalculate the age and replace the Years.bl and Month.bl with their intended values.
  ADNI_filter$age <- ADNI_filter$AGE + time_length( difftime(as.Date(ADNI_filter$EXAMDATE) , as.Date(ADNI_filter$EXAMDATE.bl) ) , "years")
  ADNI_filter$Years.bl <- time_length( difftime(as.Date(ADNI_filter$EXAMDATE) , as.Date(ADNI_filter$EXAMDATE.bl) ) , "years")
  ADNI_filter$Month.bl <- time_length( difftime(as.Date(ADNI_filter$EXAMDATE) , as.Date(ADNI_filter$EXAMDATE.bl) ) , "months")
  ADNI_filter$VISCODE_NUM = as.numeric( str_replace( str_replace( str_replace( ADNI_filter$VISCODE2 , "bl" , "0")  , "sc" , "0" ) , "m" , "") )
  
  ADNI_filter = unique(ADNI_filter)
  
  
  ADNI_filter$EcogSPEF = rowMeans(ADNI_filter[ c("EcogSPOrgan" , "EcogSPPlan" , "EcogSPVisspat"  , "EcogSPDivatt")])
  
  ADNI_filter[ !is.na(ADNI_filter$PTAU) & ADNI_filter$PTAU == "<8" , "PTAU" ] = 7.9
  ADNI_filter$PTAU = as.numeric(ADNI_filter$PTAU)
  
  ADNI_filter[ !is.na(ADNI_filter$ABETA) & ADNI_filter$ABETA == ">1700" , "ABETA" ] = 1800
  ADNI_filter$ABETA = as.numeric(ADNI_filter$ABETA)
  
  ADNI_filter[ !is.na(ADNI_filter$GDTOTAL) & ADNI_filter$GDTOTAL < 0 , "GDTOTAL"] = NA

  
  
  ADNI_filter$MMSE_ORIENTATION = rowSums( ADNI_filter[ ,c("MMDATE", "MMYEAR" ,"MMMONTH" ,"MMDAY", "MMSEASON" , "MMHOSPIT" , "MMFLOOR" , "MMCITY" , "MMAREA" , "MMSTATE")])
  ADNI_filter$MMSE_REGISTRATION = rowSums( ADNI_filter[ ,c("MMBALL" , "MMFLAG" , "MMTREE")])
  ADNI_filter$MMSE_ATTCALC = rowSums( ADNI_filter[ ,c("MMTRIALS" , "MMW" , "MMO", "MMR", "MML", "MMD")])
  ADNI_filter$MMSE_RECALL = rowSums( ADNI_filter[ ,c("MMBALLDL" , "MMFLAGDL" , "MMTREEDL")])
  ADNI_filter$MMSE_LANGUAGE = rowSums( ADNI_filter[ ,c("MMWATCH" , "MMPENCIL" , "MMREPEAT" , "MMHAND" , "MMFOLD" , "MMONFLR" , "MMREAD" , "MMWRITE" , "MMDRAW")])
  
  ADNI_filter$ADAS_MEMORY = rowSums( ADNI_filter[ ,c("Q1" , "Q4" , "Q8")])
  ADNI_filter$ADAS_LANGUAGE = rowSums( ADNI_filter[ ,c("Q10","Q11","Q12")])
  ADNI_filter$ADAS_VISUOSPATIAL = ADNI_filter[ ,c("Q3")]
  ADNI_filter$ADAS_ORIENTATION = ADNI_filter[ ,c("Q7")]
  
  # Delayed reacall in MOCA is scored (0,3,2,1) for (wrong , correct with multiple choice , correct with hint , correct without hint)
  # her I flip the scores so that the same categroies correspond to the scores (0,1,2,3), so that we can sum them as normal with the other memory tasks
  # the code applies the following formula to each score x in these columns: x = (x + (x%2) * 2) % 4
  # Note: running this line of code again returns the order to what it was. 
  ADNI_filter[ , c("DELW1","DELW2","DELW3","DELW4","DELW5")] = ((ADNI_filter[ , c("DELW1","DELW2","DELW3","DELW4","DELW5")] +
                                                                (ADNI_filter[ , c("DELW1","DELW2","DELW3","DELW4","DELW5")] %% 2)*2) %% 4)
  # LETTERS task is counting errors, so make those scores minus points
  ADNI_filter$LETTERS = - ADNI_filter$LETTERS
  
  ADNI_filter$MOCA_MEMORY = rowSums( ADNI_filter[ ,c("IMMT1W1","IMMT1W2","IMMT1W3","IMMT1W4","IMMT1W5",
                                               "IMMT2W1","IMMT2W2","IMMT2W3","IMMT2W4","IMMT2W5",
                                               "DELW1","DELW2","DELW3","DELW4","DELW5")])
  ADNI_filter$MOCA_LANGUAGE = rowSums( ADNI_filter[ ,c("LION","RHINO","CAMEL","REPEAT1","REPEAT2","FFLUENCY")])
  ADNI_filter$MOCA_VISUOSPATIAL = rowSums( ADNI_filter[ ,c("TRAILS","CUBE","CLOCKCON","CLOCKNO","CLOCKHAN")])
  ADNI_filter$MOCA_ORIENTATION = rowSums( ADNI_filter[ ,c("DATE","MONTH","YEAR","DAY","PLACE","CITY")])
  ADNI_filter$MOCA_EF = rowSums( ADNI_filter[ ,c("DIGFOR","DIGBACK","LETTERS","SERIAL1","SERIAL2","SERIAL3","SERIAL4","SERIAL5","ABSTRAN","ABSMEAS")])
  }
  # Sex strata
  ADNI_filter_male <- ADNI_filter[ ADNI_filter$PTGENDER == "Male", ]
  ADNI_filter_female <- ADNI_filter[ ADNI_filter$PTGENDER == "Female", ]

  return( list( ADNI_male = ADNI_filter_male ,
                ADNI_female = ADNI_filter_female ) )
}

PREPROCESS_I46 <- function( ucl_table ){
  
  ucl_table = ucl_table[ , -(which(names(ucl_table) %in% c("IID.x","IID.y"))) ]
  
  names(ucl_table)[names(ucl_table) == 'lefthippo_i46p1'] = "clean_hv_left"
  names(ucl_table)[names(ucl_table) == 'righthippo_i46p1'] = "clean_hv_right"
  
  names(ucl_table)[names(ucl_table) == 'ageatscandate_i46p1'] = "age"
  names(ucl_table)[names(ucl_table) == 'spm_tiv_vol_i46p1'] = "ICV"
  
  ucl_table$ICV <- ucl_table$ICV * 1000  # fix the units difference in I46 ICVs to match the other datasets
  
  names( ucl_table)[names( ucl_table) == 'Pt_1e.08']  <-   'PRS_TH_1e.08'  
  names( ucl_table)[names( ucl_table) == 'Pt_1e.07']  <-   'PRS_TH_1e.07' 
  names( ucl_table)[names( ucl_table) == 'Pt_1e.06']  <-   'PRS_TH_1e.06' 
  names( ucl_table)[names( ucl_table) == 'Pt_1e.05']  <-   'PRS_TH_1e.05' 
  names( ucl_table)[names( ucl_table) == 'Pt_0.0001'] <-   'PRS_TH_1e.04'  
  names( ucl_table)[names( ucl_table) == 'Pt_0.001']  <-   'PRS_TH_1e.03'  
  names( ucl_table)[names( ucl_table) == 'Pt_0.01']   <-   'PRS_TH_0.01'  
  names( ucl_table)[names( ucl_table) == 'Pt_0.05']   <-   'PRS_TH_0.05'   
  names( ucl_table)[names( ucl_table) == 'Pt_0.1']    <-   'PRS_TH_0.1'    
  names( ucl_table)[names( ucl_table) == 'Pt_0.2']    <-   'PRS_TH_0.2'   
  names( ucl_table)[names( ucl_table) == 'Pt_0.4']    <-   'PRS_TH_0.4'  
  names( ucl_table)[names( ucl_table) == 'Pt_0.5']    <-   'PRS_TH_0.5'   
  names( ucl_table)[names( ucl_table) == 'Pt_0.75']   <-   'PRS_TH_0.75' 
  names( ucl_table)[names( ucl_table) == 'Pt_1']      <-   'PRS_TH_1'
  
  
  PRS = READ_PRS_TABLE("~/OneDrive - University College London/ANDRE STAGE/Datasets/ucl46_data/PGS_HV_MRC1946_overlap_v2.all_score" , "INTERSECT" )
  
  names( PRS)[names( PRS) == 'Pt_1e.08']  <-   'INTERSECT_PRS_TH_1e.08'  
  names( PRS)[names( PRS) == 'Pt_1e.07']  <-   'INTERSECT_PRS_TH_1e.07' 
  names( PRS)[names( PRS) == 'Pt_1e.06']  <-   'INTERSECT_PRS_TH_1e.06' 
  names( PRS)[names( PRS) == 'Pt_1e.05']  <-   'INTERSECT_PRS_TH_1e.05' 
  names( PRS)[names( PRS) == 'Pt_0.0001'] <-   'INTERSECT_PRS_TH_1e.04'  
  names( PRS)[names( PRS) == 'Pt_0.001']  <-   'INTERSECT_PRS_TH_1e.03'  
  names( PRS)[names( PRS) == 'Pt_0.01']   <-   'INTERSECT_PRS_TH_0.01'  
  names( PRS)[names( PRS) == 'Pt_0.05']   <-   'INTERSECT_PRS_TH_0.05'   
  names( PRS)[names( PRS) == 'Pt_0.1']    <-   'INTERSECT_PRS_TH_0.1'    
  names( PRS)[names( PRS) == 'Pt_0.2']    <-   'INTERSECT_PRS_TH_0.2'   
  names( PRS)[names( PRS) == 'Pt_0.4']    <-   'INTERSECT_PRS_TH_0.4'  
  names( PRS)[names( PRS) == 'Pt_0.5']    <-   'INTERSECT_PRS_TH_0.5'   
  names( PRS)[names( PRS) == 'Pt_0.75']   <-   'INTERSECT_PRS_TH_0.75' 
  names( PRS)[names( PRS) == 'Pt_1']      <-   'INTERSECT_PRS_TH_1'
  
  ucl_table = merge(ucl_table , PRS , by.x="nshdid" , by.y = "IID" , all.x = TRUE)
    
  ucl_gpc_1 = read.table("~/OneDrive - University College London/ANDRE STAGE/Datasets/I46_SNPweights/EA.predpc")
  ucl_gpc_2 = read.table("~/OneDrive - University College London/ANDRE STAGE/Datasets/I46_SNPweights/NA.predpc")
  
  names(ucl_gpc_1) = c("nshdid" , "label" , "num_snps" , "snp_weights_gpc_EA_1" , "snp_weights_gpc_EA_2" , "prob_NW" , "prob_SE" , "prob_AJ")
  names(ucl_gpc_2) = c("nshdid" , "label" , "num_snps" , "snp_weights_gpc_NA_1" , "snp_weights_gpc_NA_2" , "snp_weights_gpc_NA_3" , "prob_WA" , "prob_EU" , "prob_EA" , "prob_NA")
  
  ucl_table = merge(ucl_table , ucl_gpc_1 , by="nshdid")
  ucl_table = merge(ucl_table , ucl_gpc_2 , by="nshdid")
  
  
  ucl_male <- ucl_table[ grepl( "male" , ucl_table$sex ) , ]
  ucl_male$sex = 1
  ucl_female <- ucl_table[ !grepl( "male" , ucl_table$sex ) , ]
  ucl_female$sex = 0

  ucl_male$clean_hv_bilateral <- ( ucl_male$clean_hv_left + ucl_male$clean_hv_right ) / 2
  ucl_female$clean_hv_bilateral <- ( ucl_female$clean_hv_left + ucl_female$clean_hv_right ) / 2
  
  
  return( list( ucl_male = ucl_male , 
                ucl_female = ucl_female) )
}

PREPROCESS_EPAD <- function( epad_table ){
  
  epad_table = epad_table[ , c(2,3,8,40,6,7,4,5,12:18,21:25,27:39,41:ncol(epad_table)) ]
  epad_table = epad_table[ , -(which(names(epad_table) %in% c("X.1","ICV.x","ICV.y" , "eTIV.y"))) ]
  names(epad_table)[names(epad_table) == 'eTIV.x'] = "eTIV"
  
  names(epad_table)[names(epad_table) == 'biological_sex'] = "sex"
  epad_table$sex = (epad_table$sex == "m")
  
  names(epad_table)[names(epad_table) == 'Lhippo'] = "clean_hv_left"
  names(epad_table)[names(epad_table) == 'Rhippo'] = "clean_hv_right"
  
  

  names(epad_table)[names( epad_table) == 'Pt_1e.08']  <-   'INTERSECT_PRS_TH_1e.08'  
  names(epad_table)[names( epad_table) == 'Pt_1e.07']  <-   'INTERSECT_PRS_TH_1e.07' 
  names( epad_table)[names( epad_table) == 'Pt_1e.06']  <-   'INTERSECT_PRS_TH_1e.06' 
  names( epad_table)[names( epad_table) == 'Pt_1e.05']  <-   'INTERSECT_PRS_TH_1e.05' 
  names( epad_table)[names( epad_table) == 'Pt_0.0001'] <-   'INTERSECT_PRS_TH_1e.04'  
  names( epad_table)[names( epad_table) == 'Pt_0.001']  <-   'INTERSECT_PRS_TH_1e.03'  
  names( epad_table)[names( epad_table) == 'Pt_0.01']   <-   'INTERSECT_PRS_TH_0.01'  
  names( epad_table)[names( epad_table) == 'Pt_0.05']   <-   'INTERSECT_PRS_TH_0.05'   
  names( epad_table)[names( epad_table) == 'Pt_0.1']    <-   'INTERSECT_PRS_TH_0.1'    
  names( epad_table)[names( epad_table) == 'Pt_0.2']    <-   'INTERSECT_PRS_TH_0.2'   
  names( epad_table)[names( epad_table) == 'Pt_0.4']    <-   'INTERSECT_PRS_TH_0.4'  
  names( epad_table)[names( epad_table) == 'Pt_0.5']    <-   'INTERSECT_PRS_TH_0.5'   
  names( epad_table)[names( epad_table) == 'Pt_0.75']   <-   'INTERSECT_PRS_TH_0.75' 
  names( epad_table)[names( epad_table) == 'Pt_1']      <-   'INTERSECT_PRS_TH_1'
  
  names( epad_table)[names( epad_table) == 'EA_PC1']      <-   'snp_weights_gpc_EA_1'
  names( epad_table)[names( epad_table) == 'EA_PC2']      <-   'snp_weights_gpc_EA_2'
  names( epad_table)[names( epad_table) == 'NA_PC1']      <-   'snp_weights_gpc_NA_1'
  names( epad_table)[names( epad_table) == 'NA_PC2']      <-   'snp_weights_gpc_NA_2'
  names( epad_table)[names( epad_table) == 'NA_PC3']      <-   'snp_weights_gpc_NA_3'

  epad_gpc_1 = read.table("~/OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/EPAD_QC.EA.predpc")
  epad_gpc_2 = read.table("~/OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/EPAD_QC.NA.predpc")
  
  names(epad_gpc_1) = c("GWAS_ID" , "label" , "num_snps" , "snp_weights_gpc_EA_1" , "snp_weights_gpc_EA_2" , "prob_NW" , "prob_SE" , "prob_AJ")
  names(epad_gpc_2) = c("GWAS_ID" , "label" , "num_snps" , "snp_weights_gpc_NA_1" , "snp_weights_gpc_NA_2" , "snp_weights_gpc_NA_3" , "prob_WA" , "prob_EU" , "prob_EA" , "prob_NA")
  
  epad_table = merge(epad_table , epad_gpc_1[,c("GWAS_ID", "prob_NW" , "prob_SE" , "prob_AJ")] , by="GWAS_ID")
  epad_table = merge(epad_table , epad_gpc_2[,c("GWAS_ID" ,"prob_WA" , "prob_EU" , "prob_EA" , "prob_NA")] , by="GWAS_ID")
  

  epad_table$clean_hv_bilateral <- ( epad_table$clean_hv_left + epad_table$clean_hv_right ) / 2
  
  epad_table$visit = "V1"

  epad_evals = read.csv("OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/EPAD_EVAL_COLS.csv")
  
  dates = grep("assessment_date", names(epad_evals) , value = 1)
  for(date_col in dates){
    epad_evals[ is.na(epad_evals$assessment_date) , "assessment_date"] = epad_evals[ is.na(epad_evals$assessment_date) , date_col ] 
    epad_evals[ (epad_evals$assessment_date %in% "") , "assessment_date"] = epad_evals[ (epad_evals$assessment_date %in% "") , date_col ] 
  }

  epad_evals = epad_evals[ !(epad_evals$assessment_date %in% "") , ]
  
  
  epad_evals$initial_assessment_date = NA
  epad_evals[ epad_evals$visit %in% "V1", "initial_assessment_date" ] =  epad_evals[ epad_evals$visit %in% "V1", "assessment_date" ]
  epad_evals = as.data.frame(epad_evals %>% group_by(patient_id) %>% mutate(across("initial_assessment_date", ~ { coalesce(., na.omit(.)[1])})))
  
  epad_evals = epad_evals[ epad_evals$patient_id %in% epad_table$patient_id , ]
  epad_table = epad_table[ epad_table$patient_id %in% epad_evals$patient_id , ]
  
  
  epad_table = epad_table[!duplicated(epad_table$patient_id) , ]
  
  epad_table = epad_table[ !is.na(epad_table$age) , ]
  epad_table$initial_age = epad_table$age
  
  epad_table = merge(epad_table , epad_evals , by = c("patient_id","visit") , all = 1)
  
  epad_table = as.data.frame(epad_table %>% group_by(patient_id) %>% mutate(across("initial_age", ~ { ifelse(!is.na(.), ., first(na.omit(.))) })))
  
  epad_table$visit_year = time_length( difftime( as.Date(epad_table$assessment_date) ,as.Date(epad_table$initial_assessment_date)) , "years")
  
  epad_table$age = epad_table$initial_age + epad_table$visit_year
  
  
  # performed on epad and adni and saved as of 02_10_2023 version
  colunms_that_should_be_the_same = c("GWAS_ID" , "sex" , "ethnicity" , 
                                      grep( "snp_weig" , names(epad_table) , value = TRUE) ,
                                      grep( "INTERSECT" , names(epad_table) , value = TRUE) , 
                                      grep( "ICV" , names(epad_table) , value = TRUE) , 
                                      grep( "prob_" , names(epad_table) , value = TRUE))
  
  # matched by PTID, fill in all the columns that should be the same regardless of visit code. 
  epad_table = as.data.frame(epad_table %>% group_by(patient_id) %>% mutate(across( colunms_that_should_be_the_same , ~ifelse(!is.na(.), ., first(na.omit(.))))) )
  
  
  epad_table[  !is.na(epad_table$gds_total) & epad_table$gds_total > 400 , "gds_total" ] = NA
  
  epad_table[  !is.na(epad_table$rbans_total_scale) & epad_table$rbans_total_scale > 400 , "rbans_total_scale" ] = NA
  
  epad_table[  !is.na(epad_table$rbans_delayed_memory_index) & epad_table$rbans_delayed_memory_index > 400 , "rbans_delayed_memory_index" ] = NA
  
  epad_table[  !is.na(epad_table$rbans_immediate_memory_index) & epad_table$rbans_immediate_memory_index > 400 , "rbans_immediate_memory_index" ] = NA
  
  epad_table[  !is.na(epad_table$rbans_language_index) & epad_table$rbans_language_index > 400 , "rbans_language_index" ] = NA
  
  epad_table[  !is.na(epad_table$rbans_visuo_constructional_index) & epad_table$rbans_visuo_constructional_index > 400 , "rbans_visuo_constructional_index" ] = NA
  
  epad_table[ !is.na(epad_table$ptau_result) & epad_table$ptau_result == "<8" , "ptau_result" ] = 7.9
  epad_table$ptau_result = as.numeric(epad_table$ptau_result)
  
  epad_table[ !is.na(epad_table$ttau_result) & epad_table$ttau_result == "<80" , "ttau_result" ] = 60
  epad_table$ttau_result = as.numeric(epad_table$ttau_result)
  
  epad_table[ !is.na(epad_table$abeta_1_42_result) & epad_table$abeta_1_42_result == "<200" , "abeta_1_42_result" ] = 100
  epad_table[ !is.na(epad_table$abeta_1_42_result) & epad_table$abeta_1_42_result == ">1700" , "abeta_1_42_result" ] = 1800
  epad_table$abeta_1_42_result = as.numeric(epad_table$abeta_1_42_result)
  
  epad_table$st_trial_total = rowSums( epad_table[ , grep("st_trial_" , names(epad_table) , value = TRUE)] == "Correct")
  
  epad_table$fms_med_total = rowSums( epad_table[ , grep("fms_med_item.*mark$", names(epad_table) , value = TRUE)] == "Correct")
  
  epad_table$mmse_recall_total = rowSums( epad_table[ , grep("mmse_recall_", names(epad_table) , value = TRUE)] )
  epad_table$mmse_naming_total = rowSums( epad_table[ , grep("mmse_naming_", names(epad_table) , value = TRUE)] )
  epad_table$mmse_orientation_total = rowSums( epad_table[ , grep("mmse_orientation_", names(epad_table) , value = TRUE)])
  epad_table$mmse_calc_total = rowSums( epad_table[ , grep("mmse_calculation_", names(epad_table) , value = TRUE)])
  
  # Or redo the pre-processing done to the file above.
  #epad_base = read.csv("OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/My_EPADS.csv")
  #epad_genetics = read.csv("OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/My_EPADS_GPC_PGS.csv")
  #epad_img = read.csv("OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/FreeSurfer_crossectional_EPAD.csv")
  #epad_link = read.table("OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/linkfile.txt",header = TRUE)
  #epad_link_2 = read.csv("OneDrive - University College London/ANDRE STAGE/Datasets/epad_data/csf_to_patient_id.csv")
  
  #epad_base = epad_base[,-1]
  #epad_genetics = epad_genetics[,-1]
  #epad_link_2 = epad_link_2[,c(-3,-4,-5,-6)]
  
  #epad_link_all = merge(epad_link , epad_link_2 , by.x = "CSF_ID" , by.y = "csf_sample_id" , all.y = TRUE)
  #epad_all = merge(epad_base , epad_link_all , by.x="patient_id", by.y = "patient_id" , all.x = TRUE)
  #epad_all = merge(epad_all , epad_genetics , by.x = "GWAS_ID" , by.y = "IID" , all.x = TRUE)
  
  #epad_all[ is.na(epad_all$age_months) , "age_months"] = 0
  #epad_all$age = epad_all$age_years + (epad_all$age_months/12)
  
  #epad_all = merge(epad_all , epad_img , by.x = "patient_id" , by.y = "SubjID" , all.x = TRUE)
  
  epad_male = epad_table[ epad_table$sex, ]
  epad_female = epad_table[ !epad_table$sex, ]
  
  return( list( epad_male = epad_male , 
                epad_female = epad_female) )
  
}


MERGE_FAMILY_HISTORY_AD <- function(ukb){
  # Pulled the faily history columns and from them will estimated how many family members had AD
  # in these three columns, a value of 10	indicates Alzheimer's disease/dementia
  # X20107  Illnesses of father
  # X20110	Illnesses of mother
  # X20111	Illnesses of siblings
  
  # X1797	Father still alive
  # X2946	Father's age
  # X1807	Father's age at death
  # X1835	Mother still alive
  # X1845	Mother's age
  # X3526	Mother's age at death
  
  ukb_family_history <- read.csv("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/UKB_Family_History.csv")
  # See if mother is reported to have AD:
  # grab all the columns from the table with "20110" in there name (illness of mother)
  # check if any of them are = 10 (AD)
  # sum that number up, if the total is greater than zero, then mother had AD, and didn't otherwise
  mother_AD  <- rowSums( ukb_family_history[ , names(ukb_family_history)[grepl( "20110" , names(ukb_family_history))]] == 10 , na.rm = TRUE) > 0
  # similar for father and siblings
  father_AD  <- rowSums( ukb_family_history[ , names(ukb_family_history)[grepl( "20107" , names(ukb_family_history))]] == 10 , na.rm = TRUE) > 0
  sibling_AD <- rowSums( ukb_family_history[ , names(ukb_family_history)[grepl( "20111" , names(ukb_family_history))]] == 10 , na.rm = TRUE) > 0
  # Add these columns to the family history table before merging to full ukb table
  ukb_family_history <- cbind(ukb_family_history , mother_AD)
  ukb_family_history <- cbind(ukb_family_history , father_AD)
  ukb_family_history <- cbind(ukb_family_history , sibling_AD)
  
  # Get the age of the mother/father 
  # rather than find: if alive : "check 'age' columns" else: "check 'age at death' columns",
  # I get the maximum age across all age columns for father (ignoring NA's) and use that as the latest age
  # If the age columns checked are all NA's then will return -Inf (Also if participant doesn't know -1)
  # then I replace negative values with NA's
  father_age <- apply(ukb_family_history[ , names(ukb_family_history)[grepl( "2946|1807" , names(ukb_family_history))] ] , 1 , max , na.rm=TRUE)
  father_age[ father_age < 60 ]  <- NA
  father_age[ father_age > 101 ] <- NA
  ukb_family_history <- cbind(ukb_family_history , father_age)
  
  mother_age <- apply(ukb_family_history[ , names(ukb_family_history)[grepl( "1845|3562" , names(ukb_family_history))] ] , 1 , max , na.rm=TRUE)
  mother_age[ mother_age < 60 ]  <- NA
  mother_age[ mother_age > 101 ] <- NA
  ukb_family_history <- cbind(ukb_family_history , mother_age)
  
  # father_AD (0 , 1) and W goes from 0 to 1
  W <- (father_age - 59.5) / (100.5 - 59.5)
  # if 0 control --> following formula collapses to log(1-W) - 0.5
  father_AD_score <- father_AD
  father_AD_score[!father_AD] <- log( 1 - W[!father_AD] ) - 0.5
  # if 1 case --> following formula collapses to -log(W) + 0.5
  father_AD_score[father_AD] <- -log(W[father_AD]) + 0.5
  
  father_AD_score[ father_age < 60 ]  <- 0
  father_AD_score[ father_age > 100 ] <- 0
  
  # father_AD (0 , 1) and W goes from 0 to 1
  W <- (mother_age - 59.5) / ( 100.5 - 59.5)
  # if 0 control --> following formula collapses to log(1-W) - 0.5
  # range of score is -0.5 to -inf
  mother_AD_score <- mother_AD
  mother_AD_score[!mother_AD] <- log( 1 - W[!mother_AD] ) - 0.5
  # if 1 case --> following formula collapses to -log(W) + 0.5
  # range of score is 0.5 to inf
  mother_AD_score[mother_AD] <- -log(W[mother_AD]) + 0.5
  
  mother_AD_score[ mother_age < 60 ]  <- 0
  mother_AD_score[ mother_age > 100 ] <- 0
  
  familial_AD_score <- father_AD_score + mother_AD_score + (sibling_AD-0.5)
  
  ukb_family_history <- cbind(ukb_family_history , father_AD_score)
  ukb_family_history <- cbind(ukb_family_history , mother_AD_score)
  ukb_family_history <- cbind(ukb_family_history , familial_AD_score)
  # Merging to full ukb table
  ukb <- merge(ukb , ukb_family_history , by="eid" , all.x = TRUE) 
  
  return(ukb)
}

MERGE_GENETICS <- function( ukb_table ){
  
  # If we want PRSICE to also generate PRS bar plots, then
  # we need to export the HV in a file with 3 columns:
  # FID : the eid from our table
  # IID : again the eid from our table
  # HV : corrected bilateral hippocampal volume
  
  #full <- merge(ukb_img_ex_outliers_male , ukb_img_ex_outliers_female , by=names(ukb_img_ex_outliers_male) , all.x = TRUE , all.y = TRUE)
  #to_save <- full[ , c("eid" , "eid" , "clean_hv_bilateral")]
  #names(to_save) <- c("FID" , "IID" , "clean_hv_bilateral")
  #write.table( to_save , "/Users/mohammedjanahi/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/ukb_cal_merged_maf01.pheno", 
  #             append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE , quote = FALSE)
  
  # this table goes to PRSice where a PRS is calculated with the command:
  
  # -------->
  #         Rscript ../PRSice_mac/PRSice.R --prsice ./../PRSice_mac/PRSice \
  #                    --base ../CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL \
  #                    --target ../ukb_cal_merged_maf01 --missing SET_ZERO \
  #                    --geno 0.1 --maf 0.05 --no-regress --fastscore --all-score --score sum \
  #                    --binary-target F --thread 1 \
  #                    --beta --stat Beta --A1 Allele1 --snp RSNUMBERS --pvalue P.value \
  #                    --bar-levels 0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.05,0.1,0.2,0.4,0.5,0.75,1.0 \
  #                    --out PRS_HV_FREESURFER_40k
  #
  # <--------
  # then we read in the results of this call for analysis
  
  # PRS_PRSICE <- READ_PRS_TABLE("~/Desktop/ukb_data/PRS-FINAL.all.score")
  # PRS_PRSICE <- READ_PRS_TABLE("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/PRS_CLN_HV_40k.all.score")
  
  prs_tables <- array(              c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_FIX.all.score"                   , "FIX") , dim =c(1,2) )
  prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_HV_FREESURFER_40k.all.score"     , NA ) )
  prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/AD_PRS_HV_40k.all.score"             , "AD" ) )
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/ICV_PRS_HV_40k.all.score"            , "ICV") )
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_HV_SE_0.01.all.score"            , "SE_0.01") )
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_HV_SE_0.0087.all.score"          , "SE_0.0087"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_0_10.all.score"              , "MAF_0.1"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_0_40.all.score"              , "MAF_0.4"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_0_50.all.score"              , "MAF_0.5"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_0_5.all.score"               , "MAF_0_5"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_5_10.all.score"              , "MAF_5_10"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_10_15.all.score"             , "MAF_10_15"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_15_20.all.score"             , "MAF_15_20"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_20_25.all.score"             , "MAF_20_25"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_25_30.all.score"             , "MAF_25_30"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_30_35.all.score"             , "MAF_30_35"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_35_40.all.score"             , "MAF_35_40"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_40_45.all.score"             , "MAF_40_45"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_45_50.all.score"             , "MAF_45_50"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_BOTTOM_THIRD.all.score"      , "MAF_BOT_THIRD"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_MIDDLE_THIRD.all.score"      , "MAF_MID_THIRD"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_MAF_TOP_THIRD.all.score"         , "MAF_TOP_THIRD"))
  #prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_NC.all.score"                    , "PRS_NC"))
  prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_HV_INTERSECT_WITH_ADNI.all.score", "ADNI_INTERSECT"))
  prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_HV_INTERSECT_WITH_I46.all.score", "I46_INTERSECT"))
  prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_HV_INTERSECT_WITH_I46_AND_ADNI.all.score", "I46_ADNI_INTERSECT"))
  prs_tables <- rbind( prs_tables , c("~/OneDrive - University College London/ANDRE STAGE/Datasets/ukb_data/PGS/PRS_HV_INTERSECT_WITH_I46_ADNI_EPAD.all.score", "INTERSECT"))
  
    for( i in 1:nrow(prs_tables)){
      PRS = READ_PRS_TABLE(prs_tables[i,1] , prs_tables[i,2] )
      ukb_table <- merge(ukb_table , PRS , by.x="eid" , by.y = "IID" , all.x=TRUE )
      
    }
    return( ukb_table)
}

WHICH_TO_EXCLUDE <- function( ukb_table ){
  # directory of the excel sheet containing which columns and what codes to exclude
  ex_dir = "~/OneDrive - University College London/ANDRE STAGE/Project 2 Multi-PGS/Scripts/Exclusion Files"
  ex_spreadsheet = "Excluded_Cases.xlsx"
  
  # Excel sheet format expected: First sheet is the reference sheet that tells which column is to be 
  # filtered with which encoding list, the rest of the sheet names match the name listed in the reference sheet
  # find all sheet names inside the main excel sheet
  all_sheets = excel_sheets(paste(ex_dir,ex_spreadsheet,sep = "/"))
  
  for (sheet_name in all_sheets)
    assign( sheet_name , data.frame( read_excel(paste(ex_dir,ex_spreadsheet,sep = "/") , sheet_name) ) )
  
  ref_sheet = get( all_sheets[1] )
  
  # the prefix to all the column names in the ukbb dataset
  pre_fix = "f."
  # EXCLUDE variable will flag all rows that need to be excluded, initially all FALSE since nothing to exclude
  EXCLUDE <- rep(FALSE , nrow(ukb_table))
  for( i in 1:nrow(ref_sheet)){
     # Grab the column to exclude
     exclude_column <- ref_sheet[i,1]
     # And its corresponding codes from the sheet name listed in the third column of the reference sheet
     # I assume the codes in the code sheet are in the first column
     exclude_codes <- get( ref_sheet[i,3] )[,1]
     # add the prefix to the column name and check all columns for any codes. 
     # return the list of rows (T/F) indicating whether codes where found
     new_exclude <- WHICH_CONTAINS( ukb_table , grep( paste( pre_fix,exclude_column, sep="") , names(ukb_table) , value = TRUE) , exclude_codes)
     # combine with the full exclusion list. If any new exclusions are found, add them to full list. 
     # FALSE in the full, TRUE in the new => TRUE
     # TRUE in the full, FALSE in the new => TRUE
     EXCLUDE <- EXCLUDE | new_exclude
  }
  return(EXCLUDE)
}

WHICH_CONTAINS <- function( table , table_columns , check_for ){
  
  check_here <- table[,table_columns]
  
  # to check in each column we use apply on dimention number 2
  # in each column, we check across rows if any contain the check_for values
  # we get a table of TRUE/FALSE across all rows/columns for if that value is found in check_for
  # rowSums > 0 converts that into "this row has or doesnt have a forbidin value". 
  # and thats what we return
  ROWS = rowSums( apply( check_here, 2, function(r) r %in% check_for )) > 0 
  
  return(ROWS)
}

REGRESS_OUT_COLUMN <- function( col_name , regress_column ){
  
  mod  <- lm( col_name  ~  regress_column , na.action = na.exclude )
  mean_factor <- mean( regress_column , na.rm = TRUE)
  
  # Now shift the residuals (centered around zero) so that their values are back in the 'normal' range. 
  # we do this with the equation:  Corrected_volume = Intercept + Residuals + ( Beta of ICV * mean ICV)
  clean_col <- coefficients(mod)[1] + resid(mod) + (coefficients(mod)[2] * mean_factor)
  
  return(clean_col)
}


REGRESS_OUT_COLUMNS <- function( table , dependant_variable , covairiates_to_remove ){
  
  # build the linear regression formula
  form = paste(dependant_variable , "~" , paste0(covairiates_to_remove,collapse=" + "))
  
  # run the linear regression on the given table
  mod  <- lm( form , data = table)
  
  # find the mean value of every covariate column
  covs_mean <- colMeans( table[,covairiates_to_remove] , na.rm = TRUE)
  
  # get the beta/intercept of every covariate column
  covs_beta <- coefficients(mod)[-1]
  
  # Now shift the residuals (centered around zero) so that their values are back in the 'normal' range. 
  # corrected_dependant_variable = intercept + residuals + sum_over_covariates( covariate_intercept * mean_covariate_value )
  # note: %*% is matrix multiplication
  clean_col <- coefficients(mod)[1] + resid(mod) + c(covs_mean %*% covs_beta)
  
  return(clean_col)
}


REGRESS_OUT <- function( hv_col_name , ICV_col_name , table ){
  
  # regressing variables out. 
  mod  <- lm( table[,hv_col_name]  ~  table[,ICV_col_name]
              + as.numeric(as.Date(table[,"Scan.Date"])) ,
              na.action = na.exclude )
  
  mean_scaling_factor <- mean(table[,ICV_col_name] , na.rm = TRUE)
  mean_scan_date <- mean(as.numeric(as.Date(table[,"Scan.Date"])) , na.rm = TRUE)
  
  # Now shift the residuals (centered around zero) so their values are back in the normal hv volume range. 
  # we do this with the equation:  Corrected_volume = Intercept + Residuals + ( Beta of ICV * mean ICV)
  clean_hv <- coefficients(mod)[1] + resid(mod) + (coefficients(mod)[2] * mean_scaling_factor) + (coefficients(mod)[3] * mean_scan_date)
  
  
  return(clean_hv)
}

normalize <- function(x) {
  rr = (x - min(x,na.rm = TRUE)) / (max(x,na.rm = TRUE) - min(x,na.rm = TRUE))
  #rr = (x - mean(x,na.rm=TRUE)) / sqrt(var(x,na.rm = TRUE))
  return( rr )
}

READ_PRS_TABLE <- function( table_path , prefix=NA){
  # Read in table that has the PRS calculated at different p-value thresholds.
  #PRS <- read.table("~/Desktop/ukb_data/PRS-FINAL.all.score" , header = TRUE)
  #PRS <- read.table("~/Desktop/UKB NOMOGRAM PROJECT/ukb_data/with_imaging_40k/PRS_CLN_HV_40k.all.score" , header = TRUE)
  PRS <- read.table(table_path , header = TRUE)
  #PRS[ ,-c(1,2)] <- lapply( PRS[,-c(1,2)] , normalize )
  PRS <- PRS[,-1]
  prefix = ifelse( is.na(prefix) , "" ,  paste0(prefix,"_") )
  # rename columns before returning 
  names(PRS)[names(PRS) == 'X1e.08']  <- paste0( prefix , "PRS_TH_1e.08")
  names(PRS)[names(PRS) == 'X1e.07']  <- paste0( prefix , "PRS_TH_1e.07")
  names(PRS)[names(PRS) == 'X1e.06']  <- paste0( prefix , "PRS_TH_1e.06")
  names(PRS)[names(PRS) == 'X1e.05']  <- paste0( prefix , "PRS_TH_1e.05")
  names(PRS)[names(PRS) == 'X0.0001'] <- paste0( prefix , "PRS_TH_1e.04")
  names(PRS)[names(PRS) == 'X0.001']  <- paste0( prefix , "PRS_TH_1e.03")
  names(PRS)[names(PRS) == 'X0.01']   <- paste0( prefix , "PRS_TH_0.01")
  names(PRS)[names(PRS) == 'X0.05']   <- paste0( prefix , "PRS_TH_0.05")
  names(PRS)[names(PRS) == 'X0.1']    <- paste0( prefix , "PRS_TH_0.1")
  names(PRS)[names(PRS) == 'X0.2']    <- paste0( prefix , "PRS_TH_0.2")
  names(PRS)[names(PRS) == 'X0.4']    <- paste0( prefix , "PRS_TH_0.4")
  names(PRS)[names(PRS) == 'X0.5']    <- paste0( prefix , "PRS_TH_0.5")
  names(PRS)[names(PRS) == 'X0.75']   <- paste0( prefix , "PRS_TH_0.75")
  names(PRS)[names(PRS) == 'X1']      <- paste0( prefix , "PRS_TH_1")
  
  return(PRS)
  
}

MAD_FILTER <- function( table , columns , threshold){
  
  for (col in columns) {
    vals <- table[,col]
    mean<- mean(vals , na.rm = TRUE)
    mae <- sum(abs(vals - mean ) , na.rm = TRUE)  / sum(!is.na(vals))
    table <- table[ which( abs(vals- mean) < 5 * mae) , ]
  }
  return (table)
  
  #ukb_img_ex_hv_left_mean  <- mean(ukb_img_ex$HV.Left , na.rm = TRUE)
  # Then the MAE: 1/n * sum( x - mean(x) )
  #ukb_img_ex_hv_left_mae  <- sum(abs(ukb_img_ex$HV.Left - ukb_img_ex_hv_left_mean ) , na.rm = TRUE)  / sum(!is.na(ukb_img_ex$HV.Left))
  
  # Repeat for right hemisphere
  #ukb_img_ex_hv_right_mean <- mean(ukb_img_ex$HV.Right , na.rm = TRUE)
  #ukb_img_ex_hv_right_mae <- sum(abs(ukb_img_ex$HV.Right - ukb_img_ex_hv_right_mean ) , na.rm = TRUE) / sum(!is.na(ukb_img_ex$HV.Right))
  
  # Do the outlier filtering
  #ukb_img_ex_outliers <- ukb_img_ex[ which( (abs(ukb_img_ex$HV.Right - ukb_img_ex_hv_right_mean) < 5 * ukb_img_ex_hv_right_mae)
  #                                          | (abs(ukb_img_ex$HV.Left - ukb_img_ex_hv_left_mean) < 5 * ukb_img_ex_hv_left_mae) ) , ]
  
  
}




PREPROCESS_GIF_FILE <- function( GIF_path ){
  
  GIF <- read.csv(GIF_path)
  
  # some column names are in the first row (all except first three)
  # so use them to rename the columns of the table
  names(GIF)[-(1:3)] = GIF[1,-(1:3)]
  names(GIF)[names(GIF) == 'Left Hippocampus'] <- 'HV.Left.GIF'
  names(GIF)[names(GIF) == 'Right Hippocampus'] <- 'HV.Right.GIF'
  
  #now remove the first row
  GIF = GIF[-1,]
  
  #cloumn 52 seams to be all NA's. so remove it.
  GIF = GIF[,-c(1,52)]
  
  # now coerce all numeric entries to numbers (all except ID's)
  GIF[,-1] = sapply(GIF[,-1],as.double)
  
  
  # make the ID's formated like the ukbb eid's
  GIF$eid = as.double(substr(GIF$ID,start=4 , stop=30))
  
  return(GIF)
}