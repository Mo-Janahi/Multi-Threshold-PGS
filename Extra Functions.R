# File Containing all extra or exploratory functions that deal with Nomograms
# By: Mohammed Janahi
#
#



GPR_PREDICT_OLD <- function( ukb=NA , train_set = array( NA , dim = c(1,1)) , test_set = array( NA , dim = c(1,1)) , x_cols=c("age") , y_col="clean_hv_bilateral" , bin_training = TRUE , training_limit = NA , plt = TRUE){
  
  # if no training limit is specified, then the limit depends on the other inputs
  if( is.na(training_limit))
    training_limit = min ( nrow(train_set) , nrow(ukb)*.5 , na.rm = TRUE)
  
  # if no train/test set is given, then I expect 'ukb' table to be given, 
  # split that given table into training (50%) and testing (50%) sets
  if( is.na(test_set[1,1]) ){
    
    # remove any rows that are missing values where we need them to be complete rows. 
    ukb = ukb[ complete.cases(ukb[, c(x_cols,y_col) ]) , ]
    
    # find 50% of the sample size
    smp_size <- floor(0.5 * nrow(ukb))
    ## set the seed to make the partition reproducible
    #set.seed(123)
    train_ind <- sample(seq_len(nrow(ukb)), size = smp_size)
    # save X_train and X_test as matrices since it allows operations to be done when 1 column or more
    ukb_train <- ukb[train_ind[1:training_limit] , c(x_cols,y_col) ]
    ukb_test <- ukb[-train_ind, c(x_cols,y_col) ]
  }
  else{
    ukb_train <- train_set[ , c(x_cols,y_col)]
    ukb_test <- test_set[ ,  c(x_cols,y_col)]
  }
  # Now that train/test split handled,
  # If more than one input columns exist, scale them to the range of the first column (expected to be age)
  all_rows = rbind(ukb_train , ukb_test)
  
  max_age <-  max(all_rows[,x_cols[1]] , na.rm = TRUE)
  min_age <- min(all_rows[,x_cols[1]] , na.rm = TRUE)
  
  if(length(x_cols) > 1 ){
    # loop over non-age columns and scale each
    for( i in 2:(length(x_cols)) ){
      
      rescaled_col = rescale(all_rows[,x_cols[i]] , to = c(min_age , max_age) )
      
      # save the scales values in a new column ("old_column_name"_SCALED)
      new_col_name = paste(x_cols[i],"_SCALED",sep="")
      
      ukb_train[ , new_col_name ] <- rescaled_col[ 1:nrow(ukb_train) ]
      ukb_test[ , new_col_name ]  <- rescaled_col[-(1:nrow(ukb_train)) ]
      
      if(bin_training){
        # for the training set, bin the extra columns by 5 percentiles
        # this generally increased the training size for each point 
        # and makes the training less sensitive to the specific training set
        five_percent <- (max_age - min_age)/20
        ukb_train[ , new_col_name ] <- round(  ukb_train[ , new_col_name ] / five_percent ) * five_percent
      }
    }
    # Also replace the names of the columns to be used for training to these scaled columns
    x_cols[-1] = paste(x_cols[-1] , "_SCALED" , sep = "")
  }
  
  
  X_train <- as.matrix( ukb_train[1:training_limit, x_cols] , ncol = length(x_cols) , dimnames = list(NULL , x_cols))
  Y_train <- ukb_train[1:training_limit, y_col ]
  
  X_test <- as.matrix( ukb_test[, x_cols] , ncol = length(x_cols) , dimnames = list(NULL , x_cols))
  Y_test <- ukb_test[ , y_col ]
  
  # initialize hyper parameters of the GPR
  g2 <- garg(list(mle = TRUE , start = mean(Y_train,na.rm = TRUE)), Y_train )
  d2 <- darg(list(mle = TRUE, start = mean(dist(X_train)) ), X_train )
  # train the GPR model
  gp  <- myGPsep(X_train ,Y_train , d=d2$start , g= g2$start , dK=TRUE) 
  mle <- jmleGPsep(  gp$gpsepi , drange=c(d2$min,d2$max) , grange=c(g2$min, g2$max) , dab = d2$ab , gab=g2$ab)
  
  # keep a copy -on the R side- of the optimized hyperparameters.
  gp$d <- as.double(mle[  , grep("d\\.|d$",names(mle),value=TRUE)]) # match "d" or "d.x" for cases with >1 dimension
  gp$g <- mle$g
  
  #keep track of the scaling limits to allow future inputs to be scaled for prediction
  gp$input_range = c(min_age , max_age)
  gp$X_ncol = length(x_cols)
  
  llik = llikGPsep(gp$gpsepi)
  p_test <- myaGPsep( X_train , Y_train , XX = X_train , d = gp$d , g = gp$g , verb = 0)
  p_means <- p_test$mean
  p_sds <- sqrt(p_test$var)
  z_scores <- (Y_train - p_means) / p_sds
  
  # Find the z-scores for the training set
  #p_test <- predGPsep( gp$gpsepi , X_train )
  #p_means <- p_test$mean
  #p_sds <- sqrt(diag(p_test$Sigma))
  #z_scores <- (Y_train - p_means) / p_sds
  
  # place them in the training table
  ukb_train$zscores <- z_scores 
  ukb_train$predicted_vals <- p_means
  ukb_train$predicted_sds <- p_sds
  # keep track of performance metrics
  train_bic <- mean( (-2 * llik) + (log( dim(X_train)[1] ) * dim(X_train)[2]) )
  train_rmse <- sqrt(mean((Y_train - p_means)^2))
  train_cor <- cor(Y_train , p_means)
  
  p_test <- myaGPsep( X_train , Y_train , XX = X_test , d = gp$d , g = gp$g , verb = 0)
  p_means <- p_test$mean
  p_sds <- sqrt(p_test$var)
  z_scores <- (Y_test - p_means) / p_sds
  
  # Find the z-scores for the test set
  #p_test <- predGPsep( gp$gpsepi , X_test )
  #p_means <- p_test$mean
  #p_sds <- sqrt(diag(p_test$Sigma))
  
  ukb_test$zscores <- z_scores
  ukb_test$predicted_vals <- p_means
  ukb_test$predicted_sds <- p_sds
  test_bic <- mean( (-2 * llik) + (log( dim(X_test)[1] ) * dim(X_test)[2]) )
  test_rmse <- sqrt(mean((Y_test - p_means)^2))
  test_cor <- cor(Y_test , p_means)
  
  # plot the GPR
  if(plt){
    num_dims = length(x_cols)
    input_list = list()
    L  <- 100
    input_list = append(input_list ,  list(seq( min(X_train[,1],na.rm = TRUE), max(X_train[,1],na.rm = TRUE), length = L )) )
    if(num_dims == 1){
      # if 1D then plot nomograms as in project 1
      bins = GPR_MODEL_TO_BINS(gp , XX=matrix(unlist(input_list[1])))
      PLOT_NOMOGRAM( bins )
    }
    else{
      input_list = append(input_list ,  list(seq( min(X_train[,2],na.rm = TRUE), max(X_train[,2],na.rm = TRUE), length = L )) )
      if(num_dims > 2)
        for(i in 3:num_dims){
          if(length(unique(X_train[,i])) < 3)
            new_list = list(unique(rescale(X_train[,i] , to = c(min_age , max_age) ) ) )
          else
            new_list = list(mean(X_train[,i],na.rm = TRUE))
          
          input_list = append(input_list ,  new_list)
          
        }
      XX <- expand.grid( input_list )
      
      p <- predGPsep( gp$gpsepi , XX )
      
      fig <- plot_ly(x =~unlist(input_list[2]) ,
                     y =~unlist(input_list[1]) ,
                     z =~matrix(p$mean[1:(length(unlist(input_list[1])) * length(unlist(input_list[2]))) ] , ncol = L ),
                     showscale = FALSE )  %>%
        add_surface() %>% 
        layout(scene = list(xaxis = list(title = x_cols[2]),
                            yaxis = list(title = x_cols[1]),
                            zaxis = list(title = y_col, range = c(2700,5300) )))
      print(fig)
    }
  }
  
  
  # Once done with predictions remove GP object from the C-side back-end since leaving it causes memory leaks
  deleteGPsep(gp$gpsepi)
  
  R = list(gp = gp,
           training_set = ukb_train,
           training_bic = train_bic,
           training_rmse = train_rmse,
           training_cor = train_cor ,
           test_set = ukb_test ,
           testing_bic = test_bic,
           testing_rmse = test_rmse,
           testing_cor = test_cor )
  # return the training and testing sets
  return ( R )
}


GPR_PREDICT <- function(train_set , test_set = NA , test_set_2 = NA , x_cols=c("age") , y_col="clean_hv_bilateral_nc_gpc_icv" , out_sample_sizes = c(0,0,0) ){
  # Re-scaling / Z-Scoring the x columns
  ukb_train  <- train_set[ , x_cols , drop = FALSE]
  ukb_test   <- test_set[ ,  x_cols ,drop = FALSE]
  other_test <- test_set_2[, x_cols , drop = FALSE]
  
  means = colMeans(rbind(ukb_train,ukb_test,other_test) , na.rm = TRUE)
  sds = colSds(as.matrix(rbind(ukb_train,ukb_test,other_test)) , na.rm = TRUE)
  

  # if( length(x_cols) > 1){
  #   max_age <-  max(rbind(ukb_train,ukb_test)[,x_cols[1]] , na.rm = TRUE)
  #   min_age <- min(rbind(ukb_train,ukb_test)[,x_cols[1]] , na.rm = TRUE)
  #   for(i in 2:length(x_cols)){
  #     ukb_train[, paste0(x_cols[i],"_RESCALED") ] = rescale( ukb_train[, x_cols[i] ] , to = c(min_age , max_age) )
  #     ukb_test[, paste0(x_cols[i],"_RESCALED") ] = rescale( ukb_test[, x_cols[i] ] , to = c(min_age , max_age) )
  #     other_test[, paste0(x_cols[i],"_RESCALED") ] = rescale( other_test[, x_cols[i] ] , to = c(min_age , max_age) )
  #   }
  # }
  # 
  # x_cols = c(x_cols[1],paste0(x_cols[-1],"_RESCALED"))
  # 
  
  x_cols = paste0(x_cols,"_RESCALED")
  
  ukb_train[, x_cols]  = sweep( sweep(ukb_train  , 2 , means) , 2 , sds , "/")
  ukb_test[, x_cols]   = sweep( sweep(ukb_test   , 2 , means) , 2 , sds , "/")
  other_test[, x_cols] = sweep( sweep(other_test , 2 , means) , 2 , sds , "/")
  
  ukb_train[,y_col]   <- train_set[,y_col]
  ukb_test[,y_col]    <- test_set[,y_col]
  other_test[,y_col]  <- test_set_2[,y_col]

  X_train <- as.matrix( ukb_train[, x_cols] , ncol = length(x_cols) , dimnames = list(NULL , x_cols))
  Y_train <- ukb_train[, y_col ]
  
  X_test <- as.matrix( ukb_test[, x_cols] , ncol = length(x_cols) , dimnames = list(NULL , x_cols))
  Y_test <- ukb_test[ , y_col ]
  
  X_test_2 <- as.matrix( other_test[, x_cols] , ncol = length(x_cols) , dimnames = list(NULL , x_cols))
  Y_test_2 <- other_test[ , y_col ]
  
  #keep track of the scaling limits to allow future inputs to be scaled for prediction
  gp = list()
  gp$X_ncol = length(x_cols)
  gp$scaling_means = means
  gp$scaling_sds = sds
  
  # initialize hyper parameters of the GPR
  #g2 <- garg(list(mle = TRUE , start = mean(Y_train,na.rm = TRUE)), Y_train )
  #d2 <- darg(list(mle = TRUE, start = mean(dist(X_train)) ), X_train )
  # train the GPR model
  #gp  <- myGP(X_train ,Y_train , d=d2$start , g= g2$start , dK=TRUE) 
  #mle <- jmleGPsep(  gp$gpsepi , drange=c(d2$min,d2$max) , grange=c(g2$min, g2$max) , dab = d2$ab , gab=g2$ab)
  
  # keep a copy -on the R side- of the optimized hyperparameters.
  #gp$d <- as.double(mle[  , grep("d\\.|d$",names(mle),value=TRUE)]) # match "d" or "d.x" for cases with >1 dimension
  #gp$g <- mle$g
  p_test <- aGP( X_train , Y_train , XX = rbind(X_train,X_test,X_test_2) , method="nn" , d = list(mle=TRUE , start = 500 , min = 500 , max = 1000) , verb = 0, omp.threads = 16)
  #print( paste("time elapsed:" , p_test$time , "seconds"))

  #p_test <- aGPsep( X_train , Y_train , XX = rbind(X_train,X_test,X_test_2) , method="alc" , d=p_test$mle$d , g=p_test$mle$g , verb = 0, omp.threads = 16)
  #print( paste("time elapsed:" , p_test$time , "seconds"))
  
  # p_test <- aGP( X_train , Y_train , XX = rbind(X_train,X_test,X_test_2) , method="nn" , d = list(mle=TRUE) , g = list(mle=TRUE) , verb = 0, omp.threads = 16)
  # print( paste("time elapsed:" , p_test$time , "seconds"))
  # 
  # p_test <- aGP( X_train , Y_train , XX = rbind(X_train,X_test,X_test_2) , method="alc" , d=p_test$mle$d , g=p_test$mle$g , verb = 0, omp.threads = 16)
  # print( paste("time elapsed:" , p_test$time , "seconds"))
  # 
  #print(sqrt(mean((p_test$mean - Y_train)^2)))
  
  #p_test <- aGP( X_train , Y_train , XX = rbind(X_train,X_test,X_test_2) , method="mspe" , d=p_test$mle$d , verb = 0, omp.threads = 16)
  #print( paste("time elapsed:" , p_test$time , "seconds"))
  #print(sqrt(mean((p_test$mean - Y_train)^2)))
  
  X_train_ind = 1:nrow(X_train)
  X_test_ind =  1:nrow(X_test) + nrow(X_train)
  X_test_2_ind = 1:nrow(X_test_2) + nrow(X_train) + nrow(X_test)
  
  X_test_2_ind_1 = 1:out_sample_sizes[1] 
  X_test_2_ind_2 = 1:out_sample_sizes[2] + out_sample_sizes[1]
  X_test_2_ind_3 = 1:out_sample_sizes[3] + out_sample_sizes[1] + out_sample_sizes[2]
  
  
  ukb_train$predicted_means <- p_test$mean[X_train_ind]
  ukb_train$predicted_sds <- sqrt(p_test$var[X_train_ind])
  ukb_train$zscores <- (Y_train - ukb_train$predicted_means ) / ukb_train$predicted_sds
  
  ukb_test$predicted_means <- p_test$mean[X_test_ind]
  ukb_test$predicted_sds <- sqrt(p_test$var[X_test_ind])
  ukb_test$zscores <- (Y_test - ukb_test$predicted_means ) / ukb_test$predicted_sds
  
  other_test$predicted_means <- p_test$mean[X_test_2_ind]
  other_test$predicted_sds <- sqrt(p_test$var[X_test_2_ind])
  other_test$zscores <- (Y_test_2 - other_test$predicted_means ) / other_test$predicted_sds
  
  
  other_test_1 = other_test[ X_test_2_ind_1 , ]
  other_test_2 = other_test[ X_test_2_ind_2 , ]
  other_test_3 = other_test[ X_test_2_ind_3 , ]
  
  # keep track of performance metrics
  # train_bic <- mean( (-2 * llik) + (log( dim(X_train)[1] ) * dim(X_train)[2]) )
  train_rmse <- sqrt(mean((Y_train - ukb_train$predicted_means)^2))
  train_cor <- cor(Y_train , ukb_train$predicted_means)
  
  test_rmse <- sqrt(mean((Y_test - ukb_test$predicted_means)^2))
  test_cor <- cor(Y_test , ukb_test$predicted_means)
  
  test_2_rmse <- sqrt(mean((Y_test_2 - other_test$predicted_means)^2))
  test_2_cor <- cor(Y_test_2 , other_test$predicted_means)

  
  test_2_rmse_1 <- sqrt(mean((Y_test_2[X_test_2_ind_1] - other_test_1$predicted_means)^2))
  test_2_cor_1 <- cor(Y_test_2[X_test_2_ind_1] , other_test_1$predicted_means)
  
  test_2_rmse_2 <- sqrt(mean((Y_test_2[X_test_2_ind_2] - other_test_2$predicted_means)^2))
  test_2_cor_2 <- cor(Y_test_2[X_test_2_ind_2] , other_test_2$predicted_means)
  
  test_2_rmse_3 <- sqrt(mean((Y_test_2[X_test_2_ind_3] - other_test_3$predicted_means)^2))
  test_2_cor_3 <- cor(Y_test_2[X_test_2_ind_3] , other_test_3$predicted_means)
  
  #test_2_rmse <- mean( sqrt(mean((other_test_1$clean_hv_bilateral_nc_gpc_icv - other_test_1$predicted_means)^2)) ,
  #                     sqrt(mean((other_test_3$clean_hv_bilateral_nc_gpc_icv - other_test_3$predicted_means)^2)) )
  
  #test_2_cor <- mean( cor( other_test_1$predicted_means , other_test_1$clean_hv_bilateral_nc_gpc_icv ) ,
  #                    cor( other_test_3$predicted_means , other_test_3$clean_hv_bilateral_nc_gpc_icv) )
  
   # plot(Y_train , ukb_train$predicted_means , xlim = c(min(c(Y_train,Y_test,Y_test_2)) , max(c(Y_train,Y_test,Y_test_2))) , ylim = c( min(c(ukb_train$predicted_means , ukb_test$predicted_means , other_test$predicted_means)) , max(c(ukb_train$predicted_means , ukb_test$predicted_means , other_test$predicted_means)) ) , pch = 19)
   # points(Y_test , ukb_test$predicted_means , col = rgb(red = 0 , green = 0 , blue = 1 , alpha = 0.5) , pch = 19)
   # points(Y_test_2 , other_test$predicted_means , col = rgb(red = 1 , green = 0 , blue = 0 , alpha = 0.5) , pch = 19)
   # points(other_test[ other_test$age > max(ukb_train$age) ,"clean_hv_bilateral_nc_gpc_icv"] , other_test[ other_test$age > max(ukb_train$age) , "predicted_means"] , col = "yellow")
   # 
   # abline(lm(other_test[ other_test$age <= max(ukb_train$age) , "predicted_means"] ~ other_test[ other_test$age <= max(ukb_train$age) ,"clean_hv_bilateral_nc_gpc_icv"]) , col = "red")
   # abline(lm(other_test$predicted_means ~ Y_test_2) , col = "yellow")
   # abline(lm(ukb_test$predicted_means ~ Y_test) , col = "blue")
   # abline(lm(ukb_train$predicted_means ~ Y_train) )
   # 
   # TEST = (merge(other_test , test_set_2[,c("age","DX")] , by = "age"))
   # points(TEST[ TEST$DX %in% c("Dementia","MCI","e4/e4") ,"clean_hv_bilateral_nc_gpc_icv"] , TEST[ TEST$DX %in% c("Dementia","MCI","e4/e4") , "predicted_means"] , col = "yellow")
   # 
   # 
   # plot( ukb_train$age , Y_train , ylim = c(min(c(Y_train,Y_test,Y_test_2)) , max(c(Y_train,Y_test,Y_test_2))) , xlim = c( min(c(ukb_train$age , ukb_test$age , other_test$age)) , max(c(ukb_train$age , ukb_test$age , other_test$age)) ) , pch = 19)
   # points( ukb_test$age , Y_test , col = rgb(red = 0 , green = 0 , blue = 1 , alpha = 0.5) , pch = 19)
   # points( other_test_1$age , Y_test_2[X_test_2_ind_1] , col = rgb(red = 1 , green = 0 , blue = 0 , alpha = 0.5) , pch = 19)
   # points( other_test_2$age , Y_test_2[X_test_2_ind_2] , col = rgb(red = 0 , green = 1 , blue = 0 , alpha = 0.5) , pch = 19)
   # points( other_test_3$age , Y_test_2[X_test_2_ind_3] , col = rgb(red = 1 , green = 0 , blue = 1 , alpha = 0.5) , pch = 19)
   # 
   # plot(other_test$age , other_test[,y_col] , col = "red")
   # points(other_test$age , other_test$predicted_means)
   # plot(ukb_test$age , ukb_test[,y_col] , col = "red")
   # points(ukb_test$age , ukb_test$predicted_means)
   # 
  print(paste0("Training RMSE:  ",train_rmse, "Testing RMSE: ", test_rmse , "Out of Sample RMSE: ",test_2_rmse))
  # gp$ m n X Z d g dK gpsepi input range x_ncol
  R = list(gp = gp,
           training_set = ukb_train,
           training_rmse = train_rmse,
           training_cor = train_cor ,
           test_set = ukb_test ,
           testing_rmse = test_rmse,
           testing_cor = test_cor,
           test_2_set = other_test ,
           testing_2_rmse = test_2_rmse,
           testing_2_cor = test_2_cor ,
           
           testing_2_rmse_1 = test_2_rmse_1,
           testing_2_cor_1 = test_2_cor_1 , 
           testing_2_rmse_2 = test_2_rmse_2,
           testing_2_cor_2 = test_2_cor_2 , 
           testing_2_rmse_3 = test_2_rmse_3,
           testing_2_cor_3 = test_2_cor_3  
  )
  # return the training and testing sets
  return ( R )
}


GPR_ZSCORE <- function( gp , test_set , pred_col , covs){
  # make a new gp object
  gpsepi = newGPsep( matrix(gp$X , ncol=gp$X_ncol) , gp$Z , gp$d , gp$g , gp$dK)
  
  #scale the input covariates to the range that was used to train the gp
  if(length(covs)>1)
  for( cov in covs[-1]){
    test_set[,cov] = rescale( test_set[,cov], to = gp$input_range )
  }
  # predict the zscores of the test set
  p_test <- predGPsep(gpsepi , as.matrix(test_set[,covs] , ncol=gp$X_ncol) )
  p_means <- p_test$mean
  p_sds <- sqrt(diag(p_test$Sigma))
  z_scores <- (test_set[,pred_col] - p_means) / p_sds
  
  deleteGPsep(gpsepi)
  
  return(z_scores)
}


PLOT_GPR <- function( gp , covs ){
  
  # make a new gp object
  gp = myGPsep( matrix(gp$X , ncol=gp$X_ncol) , gp$Z , gp$d , gp$g , gp$dK)
  
  num_dims = length(covs)
  input_list = list()
  L  <- 100
  input_list = append(input_list ,  list(seq( gp$input_range[2] , gp$input_range[2], length = L )) )
  if(num_dims == 1){
    # if 1D then plot nomograms as in project 1
    bins = GPR_MODEL_TO_BINS(gp , XX=matrix(unlist(input_list[1])))
    PLOT_NOMOGRAM( bins )
  }
  else{
    input_list = append(input_list ,  list(seq( min(X_train[,2],na.rm = TRUE), max(X_train[,2],na.rm = TRUE), length = L )) )
    if(num_dims > 2)
      for(i in 3:num_dims){
        if(length(unique(X_train[,i])) < 3)
          new_list = list(unique(rescale(X_train[,i] , to = c(min_age , max_age) ) ) )
        else
          new_list = list(mean(X_train[,i],na.rm = TRUE))
        
        input_list = append(input_list ,  new_list)
        
      }
    XX <- expand.grid( input_list )
    
    p <- predGPsep( gp$gpsepi , XX )
    
    fig <- plot_ly(x =~unlist(input_list[2]) ,
                   y =~unlist(input_list[1]) ,
                   z =~matrix(p$mean[1:(length(unlist(input_list[1])) * length(unlist(input_list[2]))) ] , ncol = L ),
                   showscale = FALSE )  %>%
      add_surface() %>% 
      layout(scene = list(xaxis = list(title = x_cols[2]),
                          yaxis = list(title = x_cols[1]),
                          zaxis = list(title = y_col, range = c(2700,5300) )))
    print(fig)
  
  return()
  }
}

GPR_MODEL_NEW <- function(ukb , x_cols=c("AGE_Latest") , y_col="clean_hv_bilateral" ){
  
  ## 70% of the sample size
  smp_size <- floor(0.7 * nrow(ukb))
  
  ## set the seed to make the partition reproducible
  #set.seed(123)
  train_ind <- sample(seq_len(nrow(ukb)), size = smp_size)
  
  ukb_train <- ukb[train_ind, ]
  ukb_test <- ukb[-train_ind, ]
  
  X = ukb_train[ !is.na(ukb[,y_col]) , x_cols]
  y = ukb_train[ !is.na(ukb[,y_col]) , y_col]
  
  X = X[1:500]
  y = y[1:500]
  
  g2 <- garg(list(mle = TRUE , start=mean(y) , min=0), y )
  d2 <- darg(list(mle = TRUE, start=mean(distance(X))), matrix(X , ncol=1) )
  
  gp  <- myGPsep(matrix(X , ncol=1) ,y , d=d2$start , g= g2$start , dK=TRUE) 
  #mle <- mleGPsep(  gp$gpsepi , param="both" , tmin=c(d2$min,g2$min) , tmax=c(d2$max, g2$max) )
  mle <- jmleGPsep(  gp$gpsepi , drange=c(d2$min,d2$max) , grange=c(g2$min, g2$max) )
  
  gp$d <- mle$theta[1]
  gp$g <- mle$theta[2]
  
  p_test <- predGPsep( gp$gpsepi , matrix(ukb_test[, x_cols] , ncol = 1) )
  p_means <- p_test$mean
  p_sds <- sqrt(diag(p_test$Sigma))
  y_test <- ukb_test[ , c(y_col)]
  z_scores <- (y_test - p_means) / p_sds
  ukb_test$scores <- z_scores
  
  return(ukb_test)
}

GPR_MODEL_2D <- function( ukb , age_column = "AGE_Latest" , THRESH , y_col="clean_hv_bilateral"  ){
  
  ukb = ukb[!is.na(ukb[,age_column]) , ]
  ukb = ukb[!is.na(ukb[,THRESH]) , ]
  
  ## 70% of the sample size
  smp_size <- floor(0.7 * nrow(ukb))
  
  ## set the seed to make the partition reproducible
  set.seed(123)
  train_ind <- sample(seq_len(nrow(ukb)), size = smp_size)
  
  ukb_train <- ukb[train_ind, ]
  ukb_test <- ukb[-train_ind, ]
  
  X = ukb_train[ , c(age_column,THRESH)]
  
  # round to the nearest 5%
  X[,THRESH] = round( X[,THRESH] / 0.5 ) * 0.5
  
  max_a <- max(X[,THRESH] , na.rm = TRUE)
  min_a <- min(X[,THRESH] , na.rm = TRUE)
  max_b <- max(X[,age_column] , na.rm = TRUE)
  min_b <- min(X[,age_column] , na.rm = TRUE)
  #  min_b <- -50
  #  max_b <- 50
  
  X[,THRESH] = ( (max_b-min_b) * ((X[,THRESH] - min_a) / (max_a-min_a)) ) + min_b
  #X[,THRESH] = X[,THRESH]*100
  
  y = ukb_train$clean_hv_bilateral
  
  X = X[1:500, ]
  y = y[1:500]
  
  #gpi <- newGPsep( as.matrix(X[1:1000 , ] ) , y[1:1000] , d=0.1, g=0.1*var(y) , dK=TRUE) 
  #eps <- sqrt(.Machine$double.eps)
  #mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  
  g2 <- garg(list(mle = TRUE , start = mean(y,na.rm = TRUE)), y )
  
  d2 <- darg(list(mle = TRUE, start = mean(dist(X)) ), X )
  # use dist function (mean of uper diag of dist matrix)
  # potentially standardize both dimentions and then rescale
  
  gp  <- myGPsep(as.matrix(X) ,y , d=d2$start , g= g2$start , dK=TRUE) 
  #mle <- mleGPsep(  gp$gpsepi , param="both" , tmin=c(d2$min,g2$min) , tmax=c(d2$max, g2$max) )
  mle <- jmleGPsep(  gp$gpsepi , drange=c(d2$min,d2$max) , grange=c(g2$min, g2$max) )
  
  gp$d <- mle$theta[1]
  gp$g <- mle$theta[2]
  
  L  <- 20
  X1 <- seq( min(X[,1],na.rm = TRUE)-10, max(X[,1],na.rm = TRUE)+10, length = L )
  X2 <- seq( min(X[,2],na.rm = TRUE)-10, max(X[,2],na.rm = TRUE)+10, length = L )
  #X2 <- min(X[,2],na.rm = TRUE)
  XX <- expand.grid(  X1 , X2 )
  
  X2_unscaled = ( (max_a-min_a) * ((X2 - min_b) / (max_b-min_b)) ) + min_a
  
  p <- predGPsep( gp$gpsepi , XX )
  
  p_test <- predGPsep( gp$gpsepi , ukb_test[, c(age_column,THRESH)] )
  p_means <- p_test$mean
  p_sds <- sqrt(diag(p_test$Sigma))
  y_test <- ukb_test$clean_hv_bilateral
  z_scores <- (y_test - p_means) / p_sds
  ukb_test$scores <- z_scores
  fig <- plot_ly(x =~rev(X1) ,
                 y =~rev(X2_unscaled) ,
                 z = ~matrix(p$mean , ncol = L),
                 showscale = FALSE )  %>%
    add_surface() %>%
    layout(scene = list(xaxis = list(title = 'AGE'),
                        yaxis = list(title = 'PRS_SCORE'),
                        zaxis = list(title = 'Bilateral HV', range = c(2700,5300)))) %>%
    add_trace( x = ~ukb_test[, c(age_column)] , 
               y = ~ukb_test[, c(THRESH)] ,
               z = ~ukb_test[, c("clean_hv_bilateral")] )
  print(fig)
  return( ukb_test )
  
  #  YY <- rmvnorm( 500, p$mean , p$Sigma)
  #  meanYY <- colMeans(YY)
  #  sdYY <- sqrt(colVars(YY))
  
  #  persp( X1 , X2 , matrix(p$mean, ncol=L) , theta=-30, phi=30, xlab="AGE", ylab="PRS_SCORE", zlab="HV_left")
  
  #matplot(XX , t(YY) , type = "l" , col="grey" , lty = 1)
  #lines(XX, p$mean, col = 2, lwd = 2)
  #points(X,y , col=rgb(0,1,0,0.1) )
  #q1 <- qnorm(0.05, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  #q2 <- qnorm(0.95, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  #lines(XX, q1, col = 4, lty = 2, lwd = 2 )
  #lines(XX, q2, col = 4, lty = 2, lwd = 2 )
  
  # Contributions of my son Jaber stepping on the keyboard to get to his toys
  #ÛS~Û~~~CC~~~~~~C ZZ§13    w3§§§§§44444444444444444
  #444444444444444444444444444444444444444444444444444
  #4444444444444444444444444444444444444
  #44444444423q131§3e§41  ae````FRRRRRRRRRRRRRRHU
  
  #  fig <- plot_ly(x =~X1 ,
  #                 y =~X2 ,
  #                 showscale = FALSE )  %>%
  #                 z = ~matrix(meanYY , ncol = L),
  #    add_surface() %>%
  #    layout(scene = list(xaxis = list(title = 'AGE'),
  #                        yaxis = list(title = 'PRS_SCORE'),
  #                        zaxis = list(title = 'Bilateral HV' , range = c(2700,5300) ))) %>%
  #    add_surface(z = ~matrix(meanYY + sdYY , ncol = L), opacity = 0.5) %>%
  #    add_surface(z = ~matrix(meanYY - sdYY , ncol = L), opacity = 0.5)
  #  fig
  
  
  
}

GPR_MODEL_2D_2 <- function( ukb , x_cols=c("AGE_Latest","PRS_TH_1e.07") , y_col="clean_hv_bilateral"  ){
  
  # revome any rows that are missing values where we need them to be complete rows. 
  ukb = ukb[ complete.cases(ukb[, c(x_cols,y_col) ]) , ]
  
  ## 70% of the sample size
  smp_size <- floor(0.7 * nrow(ukb))
  
  ## set the seed to make the partition reproducible
  #set.seed(123)
  train_ind <- sample(seq_len(nrow(ukb)), size = smp_size)
  
  ukb_train <- ukb[train_ind, ]
  ukb_test <- ukb[-train_ind, ]
  
  X = ukb_train[ , x_cols]
  
  # round all but age to the nearest 5%
  X[,-1] = round( X[,-1] / 0.5 ) * 0.5
  
  #maxs <- apply( ukb_male[ , c("PRS_TH_0.1","PRS_TH_1")] , 2 , max , na.rm=TRUE)
  #mins <- 
  max_a <- max(X[,2] , na.rm = TRUE)
  min_a <- min(X[,2] , na.rm = TRUE)
  max_b <- max(X[,1] , na.rm = TRUE)
  min_b <- min(X[,1] , na.rm = TRUE)
  #  min_b <- -50
  #  max_b <- 50
  
  X[,2] = ( (max_b-min_b) * ((X[,2] - min_a) / (max_a-min_a)) ) + min_b
  #X[,THRESH] = X[,THRESH]*100
  
  y = ukb_train$clean_hv_bilateral
  
  X = X[1:500, ]
  y = y[1:500]
  
  #gpi <- newGPsep( as.matrix(X[1:1000 , ] ) , y[1:1000] , d=0.1, g=0.1*var(y) , dK=TRUE) 
  #eps <- sqrt(.Machine$double.eps)
  #mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  
  g2 <- garg(list(mle = TRUE , start = mean(y,na.rm = TRUE)), y )
  
  d2 <- darg(list(mle = TRUE, start = mean(as.vector(dist(X))) ), X )
  # use dist function (mean of uper diag of dist matrix)
  # potentially standardize both dimentions and then rescale
  
  gp  <- myGPsep(as.matrix(X) ,y , d=d2$start , g= g2$start , dK=TRUE) 
  #mle <- mleGPsep(  gp$gpsepi , param="both" , tmin=c(d2$min,g2$min) , tmax=c(d2$max, g2$max) )
  mle <- jmleGPsep(  gp$gpsepi , drange=c(d2$min,d2$max) , grange=c(g2$min, g2$max) )
  
  gp$d <- mle$theta[1]
  gp$g <- mle$theta[2]
  
  L  <- 20
  X1 <- seq( min(X[,1],na.rm = TRUE)-10, max(X[,1],na.rm = TRUE)+10, length = L )
  X2 <- seq( min(X[,2],na.rm = TRUE)-10, max(X[,2],na.rm = TRUE)+10, length = L )
  #X2 <- min(X[,2],na.rm = TRUE)
  XX <- expand.grid(  X1 , X2 )
  
  X2_unscaled = ( (max_a-min_a) * ((X2 - min_b) / (max_b-min_b)) ) + min_a
  
  p <- predGPsep( gp$gpsepi , XX )
  
  min_a <- min(ukb_test[,2] , na.rm = TRUE)
  max_a <- max(ukb_test[,2] , na.rm = TRUE)
  
  ukb_test$scaled_x = ( (max_b-min_b) * ((ukb_test[,2] - min_a) / (max_a-min_a)) ) + min_b
  
  p_test <- predGPsep( gp$gpsepi , ukb_test[, c(x_cols[1],"scaled_x")] )
  p_means <- p_test$mean
  p_sds <- sqrt(diag(p_test$Sigma))
  y_test <- ukb_test$clean_hv_bilateral
  z_scores <- (y_test - p_means) / p_sds
  ukb_test$scores <- z_scores
  
  #fig <- plot_ly(x =~rev(X1) ,
  #               y =~rev(X2_unscaled) ,
  #               z = ~matrix(p$mean , ncol = L),
  #               showscale = FALSE )  %>%
  #  add_surface() %>%
  #  layout(scene = list(xaxis = list(title = 'AGE'),
  #                      yaxis = list(title = 'PRS_SCORE'),
  #                      zaxis = list(title = 'Bilateral HV', range = c(2700,5300)))) %>%
  #  add_trace( x = ~ukb_test[, x_cols[1] ] , 
  #             y = ~ukb_test[, x_cols[2] ] ,
  #             z = ~ukb_test[, y_col ] )
  #print(fig)
  
  return( ukb_test )
  
  #  YY <- rmvnorm( 500, p$mean , p$Sigma)
  #  meanYY <- colMeans(YY)
  #  sdYY <- sqrt(colVars(YY))
  
  #  persp( X1 , X2 , matrix(p$mean, ncol=L) , theta=-30, phi=30, xlab="AGE", ylab="PRS_SCORE", zlab="HV_left")
  
  #matplot(XX , t(YY) , type = "l" , col="grey" , lty = 1)
  #lines(XX, p$mean, col = 2, lwd = 2)
  #points(X,y , col=rgb(0,1,0,0.1) )
  #q1 <- qnorm(0.05, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  #q2 <- qnorm(0.95, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  #lines(XX, q1, col = 4, lty = 2, lwd = 2 )
  #lines(XX, q2, col = 4, lty = 2, lwd = 2 )
  
  # Contributions of my son Jaber stepping on the keyboard to get to his toys
  #ÛS~Û~~~CC~~~~~~C ZZ§13    w3§§§§§44444444444444444
  #444444444444444444444444444444444444444444444444444
  #4444444444444444444444444444444444444
  #44444444423q131§3e§41  ae````FRRRRRRRRRRRRRRHU
  
  #  fig <- plot_ly(x =~X1 ,
  #                 y =~X2 ,
  #                 showscale = FALSE )  %>%
  #                 z = ~matrix(meanYY , ncol = L),
  #    add_surface() %>%
  #    layout(scene = list(xaxis = list(title = 'AGE'),
  #                        yaxis = list(title = 'PRS_SCORE'),
  #                        zaxis = list(title = 'Bilateral HV' , range = c(2700,5300) ))) %>%
  #    add_surface(z = ~matrix(meanYY + sdYY , ncol = L), opacity = 0.5) %>%
  #    add_surface(z = ~matrix(meanYY - sdYY , ncol = L), opacity = 0.5)
  #  fig
  
  
  
}


GPR_MODEL_XD <- function(ukb , x_cols = c("AGE_Latest","PRS_TH_1","PRS_TH_1e.07","PRS_TH_0.05") , y_col="clean_hv_bilateral" ){
  
  library(laGP)
  library(Rfast)
  library(plotly)
  
  ukb = ukb[!is.na(ukb[,age_column]) , ]
  ukb = ukb[!is.na(ukb[,THRESH]) , ]
  
  X = ukb[ , x_cols]
  y = ukb$clean_hv_left
  
  X = X[!is.na(y) , ]
  y = y[!is.na(y)]
  
  y = y[!is.na(X$PRS_TH_1)]
  X = X[!is.na(X$PRS_TH_1) , ]
  
  X = X[1:500 , ]
  y = y[1:500]
  
  #gpi <- newGPsep( as.matrix(X[1:1000 , ] ) , y[1:1000] , d=0.1, g=0.1*var(y) , dK=TRUE) 
  #eps <- sqrt(.Machine$double.eps)
  #mle <- mleGPsep(gpi , param="both" , tmin=c(eps,eps) , tmax=c(10,var(y)) )
  
  g2 <- garg(list(mle = TRUE , start=mean(y,na.rm = TRUE)), y )
  
  #d2 <- darg(list(mle = TRUE, start=mean(distance(X[,1]),na.rm=TRUE)), X )
  d2 <- darg(list(mle = TRUE, start = mean(as.vector(dist(X))) ), X )
  # potentially standardize both dimentions and then rescale
  
  gp  <- myGPsep(as.matrix(X) ,y , d=d2$start , g= g2$start , dK=TRUE) 
  #mle <- mleGPsep(  gp$gpsepi , param="both" , tmin=c(d2$min,g2$min) , tmax=c(d2$max, g2$max) )
  mle <- jmleGPsep(  gp$gpsepi , drange=c(d2$min,d2$max) , grange=c(g2$min, g2$max) )
  
  gp$d <- mle$theta[1]
  gp$g <- mle$theta[2]
  
  L  <- 10
  X1 <- seq( min(X[,1],na.rm = TRUE), max(X[,1],na.rm = TRUE), length = L )
  X2 <- seq( min(X[,2],na.rm = TRUE), max(X[,2],na.rm = TRUE), length = L )
  X3 <- seq( min(X[,3],na.rm = TRUE), max(X[,3],na.rm = TRUE), length = L )
  X4 <- seq( min(X[,4],na.rm = TRUE), max(X[,4],na.rm = TRUE), length = L )
  XX <- expand.grid(  X1 , X2 , X3 , X4)
  
  p <- predGPsep( gp$gpsepi , XX )
  
  # YY <- rmvnorm( 500, p$mean , p$Sigma)
  # meanYY <- colMeans(YY)
  # sdYY <- sqrt(colVars(YY))
  
  # persp( X1 , X2 , matrix(p$mean, ncol=L) , theta=-30, phi=30, xlab="AGE", ylab="PRS_SCORE", zlab="HV_left")
  
  #matplot(XX , t(YY) , type = "l" , col="grey" , lty = 1)
  #lines(XX, p$mean, col = 2, lwd = 2)
  #points(X,y , col=rgb(0,1,0,0.1) )
  #q1 <- qnorm(0.05, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  #q2 <- qnorm(0.95, mean = p$mean, sd = sqrt(diag(p$Sigma)))
  #lines(XX, q1, col = 4, lty = 2, lwd = 2 )
  #lines(XX, q2, col = 4, lty = 2, lwd = 2 )
  
  # Contributions of my son Jaber stepping on the keyboard to get to his toys
  #ÛS~Û~~~CC~~~~~~C ZZ§13    w3§§§§§44444444444444444
  #444444444444444444444444444444444444444444444444444
  #4444444444444444444444444444444444444
  #44444444423q131§3e§41  ae````FRRRRRRRRRRRRRRHU
  
  #  fig <- plot_ly(x =~X1 ,
  #                 y =~X2 ,
  #                 z = ~matrix(meanYY , ncol = L),
  #                 showscale = FALSE )  %>%
  #    add_surface() %>%
  #    layout(scene = list(xaxis = list(title = 'AGE'),
  #                        yaxis = list(title = 'PRS_SCORE'),
  #                        zaxis = list(title = 'Bilateral HV' , range = c(2700,5300) ))) %>%
  #    add_surface(z = ~matrix(meanYY + sdYY , ncol = L), opacity = 0.5) %>%
  #    add_surface(z = ~matrix(meanYY - sdYY , ncol = L), opacity = 0.5)
  #  fig
  
  fig <- plot_ly(x =~X1 ,
                 y =~X2 ,
                 z = ~matrix(p$mean[1:1000] , ncol = L ),
                 showscale = FALSE )  %>%
    add_surface() %>%
    layout(scene = list(xaxis = list(title = 'AGE'),
                        yaxis = list(title = 'PRS_SCORE'),
                        zaxis = list(title = 'Bilateral HV', range = c(2700,5300) )))
  print(fig)
  
}


WINDOW_ANALYSIS_SEPARATED_LASSO <- function(ukb ,  sep_columns , sep_percents , percent_samples_per_bin=10 , kernel_width = 20,
                                            age_column="AGE_Latest" , hv_mean_column="clean_hv_bilateral" ,
                                            hv_left_column="clean_hv_left", hv_right_column="clean_hv_right" ){
  ukb <- ukb[ !is.na(ukb[,sep_columns[1]]) , ]
  x = as.matrix( ukb[  , sep_columns ] )
  y = as.matrix( ukb[  , c(hv_mean_column)] )
  
  n = nrow(x)
  train_rows <- sample(1:n, .66*n)
  x.train <- x[train_rows, ]
  x.test <- x[-train_rows, ]
  
  y.train <- y[train_rows]
  y.test <- y[-train_rows]
  
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1, family="gaussian")
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  
  best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
  
  #use fitted best model to make predictions
  #y_predicted <- predict(cv_model, s = best_lambda, newx = x.test)
  #mean((y.test - y_predicted)^2)
  #print(coef(best_model))
  
  W <- data.frame(matrix(ncol=1,nrow= length(sep_columns)))
  colnames(W) <- c("weights")
  #W$weights <-  (unlist(as.list(coef(best_model)))[-1] != 0 )
  W$weights <- unlist(as.list(coef(best_model)))[-1]
  
  NEW_SCORE <- rowSums(t(t(x) * W$weights))
  print(paste( "Combined score correlation with HV:" , cor(NEW_SCORE , y )))
  
  ukb[ !is.na(ukb[,sep_columns[1]]) , "NEW_SCORE"] <- NEW_SCORE
  xx <- WINDOW_ANALYSIS_SEPARATED(ukb,c("NEW_SCORE"), sep_percents , percent_samples_per_bin , kernel_width, age_column ,
                                  hv_mean_column, hv_left_column, hv_right_column)
  
  PLOT_NOMOGRAM_COMPARE(xx$upper,xx$lower)
  print(paste( "mean nomogram difference:", NOMOGRAM_DIFF_INTERPOLATE(xx$upper,xx$lower)))
  
  #return(list(upper = xx$upper, lower = xx$lower))
  
}


GPR_ANALYSIS_SEPARATED_LASSO <- function(ukb ,  sep_columns , sep_percents , XX=NA , age_column="AGE_Latest" , id_column="eid" ,
                                         hv_column="clean_hv_bilateral" , hv_left_column=NA, hv_right_column=NA , kernel_width=20 , seed = NA , plot=FALSE , print = FALSE){
  ukb <- ukb[ !is.na(ukb[,sep_columns[1]]) , ]
  x = as.matrix( ukb[  , sep_columns ] )
  y = as.matrix( ukb[  , c(hv_column)] )
  
  n = nrow(x)
  if(is.na(seed)) set.seed(4567) else set.seed(seed)
  train_rows <- sample(1:n, .5*n)
  x.train <- x[train_rows, ]
  x.test <- x[-train_rows, ]
  
  y.train <- y[train_rows]
  y.test <- y[-train_rows]
  
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1, family="gaussian")
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min
  
  best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
  
  #use fitted best model to make predictions
  #y_predicted <- predict(cv_model, s = best_lambda, newx = x.test)
  #mean((y.test - y_predicted)^2)
  
  W <- data.frame(matrix(ncol=1,nrow= length(sep_columns)))
  colnames(W) <- c("weights")
  W$weights <- unlist(as.list(coef(best_model)))[-1]
  
  NEW_SCORE <- rowSums(t(t(x) * W$weights))
  
  ukb[ !is.na(ukb[,sep_columns[1]]) , "NEW_SCORE"] <- NEW_SCORE
  
  xx <- GPR_ANALYSIS_SEPARATED(ukb[-train_rows , ] ,c("NEW_SCORE"), sep_percents,
                               XX, age_column, hv_column, hv_left_column, hv_right_column, id_column , kernel_width , P = print)
  if(plot)
    PLOT_NOMOGRAM_COMPARE(xx$upper,xx$lower , xlim = range(ukb[,age_column]))
  
  if (print){
    print(coef(best_model))
    print(paste( "Combined score correlation with HV:" , cor(NEW_SCORE , y )))
    print(paste( "mean nomogram difference:", NOMOGRAM_DIFFERENCE(xx$upper,xx$lower , age_range = range(ukb[,age_column]) )))
  }
  
  return(list(upper = xx$upper, lower = xx$lower))
  
}


LASSO_SCORE <- function( table , y_column , x_columns , print_weights=FALSE , seed = NA){
  
  x = as.matrix( table[, x_columns] )
  y = as.matrix( y_column )
  
  n = nrow(x)
  if( !is.na(seed) ) set.seed(seed)
  train_rows <- sample(1:n, .5*n)
  x.train <- x[train_rows,]
  x.test <- x[-train_rows,]
  
  y.train <- y[train_rows]
  y.test <- y[-train_rows]
  
  #perform k-fold cross-validation to find optimal lambda value
  cv_model <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1, family="gaussian")
  #find optimal lambda value that minimizes test MSE
  best_lambda <- cv_model$lambda.min

  best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)

  
  #use fitted best model to make predictions
  #y_predicted <- predict(cv_model, s = best_lambda, newx = x.test)
  #y_predicted_train <- predict(cv_model, s = best_lambda, newx = x.train)
  
  #mean(abs(y.test - y_predicted))
  #mean(abs(y.train - y_predicted_train))
  
  
  W <- data.frame(matrix(ncol=2,nrow=length(x_columns)))
  colnames(W) <- c("names","weights")
  W$names <- x_columns
  W$weights <- unlist(as.list(coef(best_model)))[-1]
  
  NEW_SCORE <- rowSums(t(t(as.matrix(table[, x_columns ])) * W$weights))
  
  if(print_weights){
    #print(paste(W$names,W$weights))
    return(list(score = NEW_SCORE , weights = W$weights))
  }
  return(NEW_SCORE)
}



GIF_IMAGE <- function(){
  
  ## GIF IMAGE INCORPERATION
  
  plot(ukb_male$clean_hv_bilateral , ukb_male$clean_hv_bilateral_gif , xlim = c(0,5500) , ylim = c(0,5000))
  abline(lm(ukb_male$clean_hv_bilateral_gif ~ ukb_male$clean_hv_bilateral ) )
  
  mean(ukb_male$clean_hv_bilateral)
  mean(ukb_male$clean_hv_bilateral_gif)
  
  
  bins_male_gif <- WINDOW_ANALYSIS(ukb_male , hv_left_column = "clean_hv_left_gif" , hv_right_column = "clean_hv_right_gif" , hv_mean_column = "clean_hv_bilateral_gif")
  bins_female_gif <- WINDOW_ANALYSIS(ukb_female , hv_left_column = "clean_hv_left_gif", hv_right_column = "clean_hv_right_gif" , hv_mean_column = "clean_hv_bilateral_gif")
  
  
  gp_model_male_mean_gif <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_bilateral_gif")
  gp_model_male_left_gif <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_left_gif")
  gp_model_male_right_gif <- GPR_MODEL(ukb_male[1:1000,] , y_col = "clean_hv_right_gif")
  
  ages_male <- seq( min(ukb_male$AGE_Latest) , max(ukb_male$AGE_Latest) , 0.25)
  bins_male_gpr_gif <- GPR_MODEL_TO_BINS(gp_model_male_mean_gif , gp_model_male_left_gif , gp_model_male_right_gif , ages_male)
  
  # or re-train models 
  gp_model_female_mean_gif <- GPR_MODEL(ukb_female[1:1000,] , y_col = "clean_hv_bilateral_gif")
  gp_model_female_left_gif <- GPR_MODEL(ukb_female[1:1000,] , y_col = "clean_hv_left_gif")
  gp_model_female_right_gif <- GPR_MODEL(ukb_female[1:1000,] , y_col = "clean_hv_right_gif")
  
  ages_female <- seq( min(ukb_female$AGE_Latest) , max(ukb_female$AGE_Latest) , 0.25)
  bins_female_gpr_gif <- GPR_MODEL_TO_BINS(gp_model_female_mean_gif , gp_model_female_left_gif , gp_model_female_right_gif , ages_female)
  
  
  bins_male_low_gif  <- WINDOW_ANALYSIS( ukb_male_low , hv_left_column = "clean_hv_left_gif", hv_right_column = "clean_hv_right_gif", hv_mean_column = "clean_hv_bilateral_gif")
  bins_male_high_gif <- WINDOW_ANALYSIS( ukb_male_high ,hv_left_column = "clean_hv_left_gif" ,hv_right_column = "clean_hv_right_gif" , hv_mean_column = "clean_hv_bilateral_gif") 
  bins_female_low_gif <- WINDOW_ANALYSIS( ukb_female_low ,hv_left_column = "clean_hv_left_gif" ,  hv_right_column = "clean_hv_right_gif" , hv_mean_column = "clean_hv_bilateral_gif") 
  bins_female_high_gif <- WINDOW_ANALYSIS( ukb_female_high , hv_left_column = "clean_hv_left_gif" , hv_right_column = "clean_hv_right_gif" , hv_mean_column = "clean_hv_bilateral_gif") 
  
  
  gp_model_male_mean_low_gif <- GPR_MODEL(ukb_male_low[1:1000,] , y_col = "clean_hv_bilateral_gif")
  gp_model_male_left_low_gif <- GPR_MODEL(ukb_male_low[1:1000,] , y_col = "clean_hv_left_gif")
  gp_model_male_right_low_gif <- GPR_MODEL(ukb_male_low[1:1000,] , y_col = "clean_hv_right_gif")
  bins_male_gpr_low_gif <- GPR_MODEL_TO_BINS(gp_model_male_mean_low_gif , gp_model_male_left_low_gif , gp_model_male_right_low_gif , ages_male)
  
  # or re-train models 
  gp_model_male_mean_high_gif <- GPR_MODEL(ukb_male_high[1:1000,] , y_col = "clean_hv_bilateral_gif")
  gp_model_male_left_high_gif <- GPR_MODEL(ukb_male_high[1:1000,] , y_col = "clean_hv_left_gif")
  gp_model_male_right_high_gif <- GPR_MODEL(ukb_male_high[1:1000,] , y_col = "clean_hv_right_gif")
  bins_male_gpr_high_gif <- GPR_MODEL_TO_BINS(gp_model_male_mean_high_gif , gp_model_male_left_high_gif , gp_model_male_right_high_gif , ages_male)
  
  
  gp_model_female_mean_low_gif <- GPR_MODEL(ukb_female_low[1:1000,] , y_col = "clean_hv_bilateral_gif")
  gp_model_female_left_low_gif <- GPR_MODEL(ukb_female_low[1:1000,] , y_col = "clean_hv_left_gif")
  gp_model_female_right_low_gif <- GPR_MODEL(ukb_female_low[1:1000,] , y_col = "clean_hv_right_gif")
  bins_female_gpr_low_gif <- GPR_MODEL_TO_BINS(gp_model_female_mean_low_gif , gp_model_female_left_low_gif , gp_model_female_right_low_gif , ages_female)
  
  # or re-train models 
  gp_model_female_mean_high_gif <- GPR_MODEL(ukb_female_high[1:1000,] , y_col = "clean_hv_bilateral_gif")
  gp_model_female_left_high_gif <- GPR_MODEL(ukb_female_high[1:1000,] , y_col = "clean_hv_left_gif")
  gp_model_female_right_high_gif <- GPR_MODEL(ukb_female_high[1:1000,] , y_col = "clean_hv_right_gif")
  bins_female_gpr_high_gif <- GPR_MODEL_TO_BINS(gp_model_female_mean_high_gif , gp_model_female_left_high_gif , gp_model_female_right_high_gif , ages_female)
  
  
  par(mfrow=c(1,2))
  
  NOMOGRAM_DIFF(bins_male , bins_male_gif)
  NOMOGRAM_DIFF(bins_female , bins_female_gif)
  PLOT_NOMOGRAM(bins_male_gpr)
  PLOT_NOMOGRAM(bins_male_gpr_gif)
  
  NOMOGRAM_DIFF(bins_male_gpr_high , bins_male_gpr_low)
  NOMOGRAM_DIFF(bins_male_gpr_high_gif , bins_male_gpr_low_gif)
  PLOT_NOMOGRAM_COMPARE(bins_male_gpr_high , bins_male_gpr_low)
  PLOT_NOMOGRAM_COMPARE(bins_male_gpr_high_gif , bins_male_gpr_low_gif)
  
}


MAF_EXPLORATIONS <- function(){
  
  ## MAF experiments
  
  #
  #
  #
  #
  
  #for( i in seq(0 , 0.45 , 0.05) )
  #write(GWAS[ (i < GWAS$MAF) & (GWAS$MAF <= i+0.05)  , "RSNUMBERS"] , 
  #                     paste("/Users/mohammedjanahi/Desktop/UKB NOMOGRAM PROJECT/ukb_data/MY-HV-GWAS-",i,"to",(i+0.05),"-MAF.snp",sep=""))
  
  
  ukb_male_copy <- ukb_male
  ukb_male <- ukb_male[ !is.na(ukb_male$MAF_0.1_PRS_TH_1) , ]
  
  par(mfrow=c(1,3))
  hist(ukb_male$MAF_0.1_PRS_TH_1)
  hist(ukb_male$MAF_0.4_PRS_TH_1)
  hist(ukb_male$MAF_0.5_PRS_TH_1)
  
  mean(ukb_male$MAF_0.1_PRS_TH_1 , na.rm = TRUE)
  mean(ukb_male$MAF_0.4_PRS_TH_1 , na.rm = TRUE)
  mean(ukb_male$MAF_0.5_PRS_TH_1 , na.rm = TRUE)
  
  ggpairs(ukb_male[ , c("MAF_0.1_PRS_TH_1","MAF_0.4_PRS_TH_1","MAF_0.5_PRS_TH_1") ])
  
  percent = nrow(ukb_male)*0.3
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_0.1_PRS_TH_1) , ]
  upper_prs_maf_10 <- tail( ukb_male , percent)
  lower_prs_maf_10 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_0.4_PRS_TH_1) , ]
  upper_prs_maf_40 <- tail( ukb_male , percent)
  lower_prs_maf_40 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_0.5_PRS_TH_1) , ]
  upper_prs_maf_50 <- tail( ukb_male , percent)
  lower_prs_maf_50 <- head( ukb_male , percent)
  
  length( intersect(upper_prs_maf_10$eid , upper_prs_maf_40$eid))
  sum( upper_prs_maf_10$eid == upper_prs_maf_50$eid)
  sum( upper_prs_maf_40$eid == upper_prs_maf_50$eid)
  
  sum( lower_prs_maf_10$eid == lower_prs_maf_40$eid)
  sum( lower_prs_maf_10$eid == lower_prs_maf_50$eid)
  sum( lower_prs_maf_40$eid == lower_prs_maf_50$eid)
  
  sum( upper_prs_maf_10$eid == lower_prs_maf_40$eid)
  sum( upper_prs_maf_10$eid == lower_prs_maf_50$eid)
  
  sum( upper_prs_maf_40$eid == lower_prs_maf_10$eid)
  sum( upper_prs_maf_40$eid == lower_prs_maf_50$eid)
  
  sum( upper_prs_maf_50$eid == lower_prs_maf_40$eid)
  sum( upper_prs_maf_50$eid == lower_prs_maf_10$eid)
  
  PLOT_NOMOGRAM_COMPARE(WINDOW_ANALYSIS(upper_prs_maf_10) , WINDOW_ANALYSIS(lower_prs_maf_10), title = paste("MAF SNPs 00-10% \n mean separation:",floor(NOMOGRAM_DIFF(WINDOW_ANALYSIS(upper_prs_maf_10) , WINDOW_ANALYSIS(lower_prs_maf_10))),"mm3"))
  PLOT_NOMOGRAM_COMPARE(WINDOW_ANALYSIS(upper_prs_maf_40) , WINDOW_ANALYSIS(lower_prs_maf_40), title = paste("MAF SNPs 10-40% \n mean separation:",floor(NOMOGRAM_DIFF(WINDOW_ANALYSIS(upper_prs_maf_40) , WINDOW_ANALYSIS(lower_prs_maf_40))),"mm3"))
  PLOT_NOMOGRAM_COMPARE(WINDOW_ANALYSIS(upper_prs_maf_50) , WINDOW_ANALYSIS(lower_prs_maf_50), title = paste("MAF SNPs 40-50% \n mean separation:",floor(NOMOGRAM_DIFF(WINDOW_ANALYSIS(upper_prs_maf_50) , WINDOW_ANALYSIS(lower_prs_maf_50))),"mm3"))
  
  upper_ids_10_40 <- intersect(upper_prs_maf_10$eid , upper_prs_maf_40$eid)
  upper_ids_10_50 <- intersect(upper_prs_maf_10$eid , upper_prs_maf_50$eid)
  upper_ids_40_50 <- intersect(upper_prs_maf_40$eid , upper_prs_maf_50$eid)
  upper_ids_all <- intersect(upper_prs_maf_10$eid , intersect( upper_prs_maf_40$eid , upper_prs_maf_50$eid))
  
  lower_ids_10_40 <- intersect(lower_prs_maf_10$eid , lower_prs_maf_40$eid)
  lower_ids_10_50 <- intersect(lower_prs_maf_10$eid , lower_prs_maf_50$eid)
  lower_ids_40_50 <- intersect(lower_prs_maf_40$eid , lower_prs_maf_50$eid)
  lower_ids_all <- intersect(lower_prs_maf_10$eid , intersect( lower_prs_maf_40$eid , lower_prs_maf_50$eid))
  
  bins_upper_10_40 <- GPR_ANALYSIS( upper_prs_maf_10[  match(upper_ids_10_40 , upper_prs_maf_10$eid) ,])
  bins_upper_10_50 <- GPR_ANALYSIS( upper_prs_maf_10[  match(upper_ids_10_50 , upper_prs_maf_10$eid) ,])
  bins_upper_40_50 <- GPR_ANALYSIS( upper_prs_maf_40[  match(upper_ids_40_50 , upper_prs_maf_40$eid) ,])
  bins_upper_all <- GPR_ANALYSIS(upper_prs_maf_10[  match(upper_ids_all , upper_prs_maf_10$eid) ,])
  
  bins_lower_10_40 <- GPR_ANALYSIS( lower_prs_maf_10[  match(lower_ids_10_40 , lower_prs_maf_10$eid) ,])
  bins_lower_10_50 <- GPR_ANALYSIS( lower_prs_maf_10[  match(lower_ids_10_50 , lower_prs_maf_10$eid) ,])
  bins_lower_40_50 <- GPR_ANALYSIS( lower_prs_maf_40[  match(lower_ids_40_50 , lower_prs_maf_40$eid) ,])
  bins_lower_all <- GPR_ANALYSIS(lower_prs_maf_10[  match(lower_ids_all , lower_prs_maf_10$eid) ,])
  
  NOMOGRAM_DIFF(bins_upper_10_40 , bins_lower_10_40)
  NOMOGRAM_DIFF(bins_upper_10_50 , bins_lower_10_50)
  NOMOGRAM_DIFF(bins_upper_40_50 , bins_lower_40_50)
  NOMOGRAM_DIFF(bins_upper_all , bins_lower_all)
  
  par(mfrow=c(1,3))
  PLOT_NOMOGRAM_COMPARE(bins_upper_10_40 , bins_lower_10_40, title = paste("MAF SNPs 00-10% & 10-40% \n mean separation:",floor(NOMOGRAM_DIFF_INTERPOLATE(bins_upper_10_40 , bins_lower_10_40, age_range = c(47,83) )),"mm3"))
  PLOT_NOMOGRAM_COMPARE(bins_upper_10_50 , bins_lower_10_50, title = paste("MAF SNPs 00-10% & 40-50% \n mean separation:",floor(NOMOGRAM_DIFF_INTERPOLATE(bins_upper_10_50 , bins_lower_10_50, age_range = c(47,83))),"mm3"))
  PLOT_NOMOGRAM_COMPARE(bins_upper_40_50 , bins_lower_40_50, title = paste("MAF SNPs 10-40% & 40-50% \n mean separation:",floor(NOMOGRAM_DIFF_INTERPOLATE(bins_upper_40_50 , bins_lower_40_50, age_range = c(47,83))),"mm3"))
  par(mfrow=c(1,1))
  PLOT_NOMOGRAM_COMPARE(bins_upper_all , bins_lower_all, title = paste("MAF SNPs 10-40% & 10-40% & 40-50%\n mean separation:",floor(NOMOGRAM_DIFF_INTERPOLATE(bins_upper_all , bins_lower_all, age_range = c(47,83))),"mm3"))
  
  
  par(mfrow=c(2,5))
  hist(ukb_male$MAF_0_5_PRS_TH_1e.03)
  hist(ukb_male$MAF_5_10_PRS_TH_1e.03)
  hist(ukb_male$MAF_10_15_PRS_TH_1e.03)
  hist(ukb_male$MAF_15_20_PRS_TH_1e.03)
  hist(ukb_male$MAF_20_25_PRS_TH_1e.03)
  hist(ukb_male$MAF_25_30_PRS_TH_1e.03)
  hist(ukb_male$MAF_30_35_PRS_TH_1e.03)
  hist(ukb_male$MAF_35_40_PRS_TH_1e.03)
  hist(ukb_male$MAF_40_45_PRS_TH_1e.03)
  hist(ukb_male$MAF_45_50_PRS_TH_1e.03)
  
  
  mean(ukb_male$MAF_0_5_PRS_TH_1e.03,na.rm = TRUE)
  mean(ukb_male$MAF_5_10_PRS_TH_1e.03,na.rm = TRUE)
  mean(ukb_male$MAF_10_15_PRS_TH_1e.03,na.rm = TRUE)
  mean(ukb_male$MAF_15_20_PRS_TH_1e.03,na.rm = TRUE)
  mean(ukb_male$MAF_20_25_PRS_TH_1e.03,na.rm = TRUE)
  mean(ukb_male$MAF_25_30_PRS_TH_1e.03,na.rm = TRUE)
  mean(ukb_male$MAF_30_35_PRS_TH_1e.03,na.rm = TRUE)
  mean(ukb_male$MAF_35_40_PRS_TH_1e.03,na.rm = TRUE)
  mean(ukb_male$MAF_40_45_PRS_TH_1e.03,na.rm = TRUE)
  mean(ukb_male$MAF_45_50_PRS_TH_1e.03,na.rm = TRUE)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_0_5_PRS_TH_1) , ]
  upper_prs_maf_0_5 <- tail( ukb_male , percent)
  lower_prs_maf_0_5 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_5_10_PRS_TH_1) , ]
  upper_prs_maf_5_10 <- tail( ukb_male , percent)
  lower_prs_maf_5_10 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_10_15_PRS_TH_1) , ]
  upper_prs_maf_10_15 <- tail( ukb_male , percent)
  lower_prs_maf_10_15 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_15_20_PRS_TH_1) , ]
  upper_prs_maf_15_20 <- tail( ukb_male , percent)
  lower_prs_maf_15_20 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_20_25_PRS_TH_1) , ]
  upper_prs_maf_20_25 <- tail( ukb_male , percent)
  lower_prs_maf_20_25 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_25_30_PRS_TH_1) , ]
  upper_prs_maf_25_30 <- tail( ukb_male , percent)
  lower_prs_maf_25_30 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_30_35_PRS_TH_1) , ]
  upper_prs_maf_30_35 <- tail( ukb_male , percent)
  lower_prs_maf_30_35 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_35_40_PRS_TH_1) , ]
  upper_prs_maf_35_40 <- tail( ukb_male , percent)
  lower_prs_maf_35_40 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_40_45_PRS_TH_1) , ]
  upper_prs_maf_40_45 <- tail( ukb_male , percent)
  lower_prs_maf_40_45 <- head( ukb_male , percent)
  
  ukb_male <- ukb_male[ order(ukb_male$MAF_45_50_PRS_TH_1) , ]
  upper_prs_maf_45_50 <- tail( ukb_male , percent)
  lower_prs_maf_45_50 <- head( ukb_male , percent)
  
  
  for (p in c(1e-08,1e-07,1e-06,1e-05,1e-04,1e-03,1e-02,0.05,0.1,0.2,0.4,0.5,0.75,1) ){
    print(paste(sum(GWAS$P.value < p & GWAS$MAF > 0.00 & GWAS$MAF <= 0.05) ,
                sum(GWAS$P.value < p & GWAS$MAF > 0.05 & GWAS$MAF <= 0.10) ,
                sum(GWAS$P.value < p & GWAS$MAF > 0.10 & GWAS$MAF <= 0.15) ,
                sum(GWAS$P.value < p & GWAS$MAF > 0.15 & GWAS$MAF <= 0.20) ,
                sum(GWAS$P.value < p & GWAS$MAF > 0.20 & GWAS$MAF <= 0.25) ,
                sum(GWAS$P.value < p & GWAS$MAF > 0.25 & GWAS$MAF <= 0.30) ,
                sum(GWAS$P.value < p & GWAS$MAF > 0.30 & GWAS$MAF <= 0.35) ,
                sum(GWAS$P.value < p & GWAS$MAF > 0.35 & GWAS$MAF <= 0.40) ,
                sum(GWAS$P.value < p & GWAS$MAF > 0.40 & GWAS$MAF <= 0.45) ,
                sum(GWAS$P.value < p & GWAS$MAF > 0.45 & GWAS$MAF <= 0.5)))
  }
  
  
  for (p in c("1e.04","1e.03","0.01","0.05","0.1","0.2","0.4","0.5","0.75","1") ){
    print(paste( cor(ukb_male[,paste("MAF_0_5_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_5_10_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_10_15_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_15_20_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_20_25_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_25_30_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_30_35_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_35_40_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_40_45_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_45_50_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs") ))
  }
  
  
  
  th = c("1e.08","1e.07","1e.06","1e.05","1e.04","1e.03","0.01","0.05","0.1","0.2","0.4","0.5","0.75","1")
  maf = c("MAF_0_5_PRS_TH_","MAF_5_10_PRS_TH_", "MAF_10_15_PRS_TH_","MAF_15_20_PRS_TH_","MAF_20_25_PRS_TH_", "MAF_25_30_PRS_TH_", "MAF_30_35_PRS_TH_", "MAF_35_40_PRS_TH_", "MAF_40_45_PRS_TH_","MAF_45_50_PRS_TH_")
  for (t1 in 1:length(th) )
    for(t2 in c(1) )
      for(m1 in 1:length(maf))
        for(m2 in c(6) )
          try(print(paste( th[t1]," at ",maf[m1]," and ",th[t2]," at ",maf[m2]," : ",cor(ukb_male[,paste(maf[m1],th[t1],sep='')] , ukb_male[,paste(maf[m2],th[t2],sep='')] , use="complete.obs") )) , silent = TRUE)
  
  
  for (p in c("1e.08","1e.04","1e.03","0.01","0.05","0.1","0.2","0.4","0.5","0.75","1") ){
    print(paste( cor(ukb_male[,paste("MAF_0.1_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_0.4_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs"),
                 cor(ukb_male[,paste("MAF_0.5_PRS_TH_",p,sep='')] , ukb_male$clean_hv_bilateral , use="complete.obs") ))
  }
  
  
  
  
  th = c("1e.08","1e.07","1e.06","1e.05","1e.04","1e.03","0.01","0.05","0.1","0.2","0.4","0.5","0.75","1")
  maf = c("MAF_0.1_PRS_TH_","MAF_0.4_PRS_TH_", "MAF_0.5_PRS_TH_")
  for ( t1 in 1:length(th) )
    for( m1 in 1:length(maf) ){
      cat(paste(maf[m1],th[t1],"\t",sep="") )
      for (t2 in 1:length(th) )
        for(m2 in 1:length(maf) ){
          tryCatch( {
            cat( cor(ukb_male[,paste(maf[m1],th[t1],sep='')] , ukb_male[,paste(maf[m2],th[t2],sep='')] , use="complete.obs") )
            cat("\t")
          },
          error=function(e){
            cat("X")
            cat("\t")
          } )
        }
      cat("\n")
    }
  
  th = c("1e.08","1e.07","1e.06","1e.05","1e.04","1e.03","0.01","0.05","0.1","0.2","0.4","0.5","0.75","1")
  maf = c("MAF_BOT_THIRD_PRS_TH_","MAF_MID_THIRD_PRS_TH_", "MAF_TOP_THIRD_PRS_TH_")
  for ( t1 in 1:length(th) )
    for( m1 in 1:length(maf) ){
      cat(paste(maf[m1],th[t1],"\t",sep="") )
      for (t2 in 1:length(th) )
        for(m2 in 1:length(maf) ){
          tryCatch( {
            cat( cor(ukb_male[,paste(maf[m1],th[t1],sep='')] , ukb_male[,paste(maf[m2],th[t2],sep='')] , use="complete.obs") )
            cat("\t")
          },
          error=function(e){
            cat("X")
            cat("\t")
          } )
        }
      cat("\n")
    }
  
  
  x =  ukb_male[ !is.na(ukb_male$PRS_TH_1) , columns ] 
  y =  ukb_male[ !is.na(ukb_male$PRS_TH_1) , c("clean_hv_bilateral")] 
  
  
  NEW_SCORE <- LASSO_SCORE( ukb_male[  !is.na(ukb_male$PRS_TH_1) , ]  , y, columns)
  ukb_male[!is.na(ukb_male$PRS_TH_1), "NEW_SCORE"] <- NEW_SCORE
  
  table = ukb_male[  !is.na(ukb_male$PRS_TH_1),]
  y = resid(lm( clean_hv_bilateral ~ AGE_Latest, data = table))
  for(col in c(columns,"NEW_SCORE")){
    cc = cor( table[,col]  , y )
    print(paste(col,"corelation with hv bilateral:",cc))
  }
  
  
  #xx <- WINDOW_ANALYSIS_SEPARATED(ukb_male,c("NEW_SCORE"),c(0.3))
  #PLOT_NOMOGRAM_COMPARE(xx$upper,xx$lower)
  #NOMOGRAM_DIFF_INTERPOLATE(xx$upper,xx$lower)
  
  # match all column names that start with "PRS_NC_TH"
  PRS_NC_columns <- names(ukb_male)[ grepl("^PRS_NC_PRS_TH" , names(ukb_male)) ]
  # match all column names that start with "MAF_" and then contain one digit and then not a 'dot'
  # this matches all MAF_XX_XX_PRS_TH_XX formats and none of the MAF_BOT or MAF_0.X formats
  MAF_detailed_columns <- names(ukb_male)[ grepl("^MAF_\\d[^.]" , names(ukb_male)) ]
  # match all column names that start with "MAF_0."
  MAF_simple_columns <- names(ukb_male)[ grepl("MAF_0\\." , names(ukb_male)) ]
  # match all column names that start with "MAF_TOP" or "MAF_BOT" or "MAF_MID"
  MAF_thirds_columns <- names(ukb_male)[ grepl("^MAF_[BOT|TOP|MID]" , names(ukb_male)) ]
  
}


GWAS_EXPLORATIONS <- function(){
  
  GWAS <- read.table( "~/OneDrive - University College London/ADNRE STAGE/Datasets/CHARGE-ENIGMA-HV-METAANALYSIS-201311141.TBL.FINAL", header = TRUE)
  
  #plot( GWAS$Beta[1:100000] , GWAS$SE[1:100000] , pch = '.' )
  #plot( GWAS$Beta[1:100000] , -log(GWAS$P.value[1:100000]) , pch = '.' )
  #plot( -log(GWAS$P.value[1:100000]) , GWAS$SE[1:100000] , pch = '.' )
  
  
  #ggplot(GWAS[1:100000,], aes(x = -log(P.value), y = SE)) + geom_point(aes(colour = abs(Beta))) + scale_colour_gradient(low="red",high="blue")
  
  HV_SNPS <- read.table("/Users/mohammedjanahi/Desktop/UKB NOMOGRAM PROJECT/ukb_data/PRSice_mac/PRS_HV_FREESURFER_40k.snp" , header = TRUE)
  
  names(HV_SNPS) <- c("CHR","RSNUMBERS","BP","P","IN_PRSICE","Base2")
  
  GWAS <- merge(GWAS , HV_SNPS[ , c("RSNUMBERS","IN_PRSICE")] , by="RSNUMBERS" , all.x = TRUE)
  
  GWAS[ is.na(GWAS$IN_PRSICE), "IN_PRSICE"] <- 0
  
  GWAS$MAF <- (-abs(GWAS$Freq1-0.5)+0.5)
  
  ggplot(GWAS[GWAS$IN_PRSICE==1,], aes(x = -log(P.value), y = SE)) + geom_point(aes(colour = IN_PRSICE)) + scale_colour_gradient(low="red",high="blue")
  
  ggplot(GWAS[1:10000,], aes(x = -log(P.value), y = SE)) + geom_point(aes(colour = abs(Beta))) + scale_colour_gradient(low="red",high="blue") + geom_point(data=GWAS[GWAS$IN_PRSICE==1,])
  
  summary( lm(Beta ~ P.value + SE , data = GWAS) )
  
  #Call:
  #  lm(formula = Beta ~ P.value + SE, data = GWAS)
  
  #Residuals:
  #  Min       1Q   Median       3Q      Max 
  #-0.54822 -0.01094  0.00004  0.01094  0.45553 
  
  #Coefficients:
  #  Estimate Std. Error t value Pr(>|t|)    
  #(Intercept) -8.241e-06  2.356e-05  -0.350   0.7265    
  #P.value     -7.246e-05  3.615e-05  -2.004   0.0451 *  
  #  SE           2.796e-03  4.794e-04   5.831 5.51e-09 ***
  #  ---
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  #
  #Residual standard error: 0.03277 on 9707618 degrees of freedom
  #Multiple R-squared:  3.887e-06,	Adjusted R-squared:  3.681e-06 
  #F-statistic: 18.86 on 2 and 9707618 DF,  p-value: 6.416e-09
  
  
  summary( lm(abs(Beta) ~ P.value + SE , data = GWAS))
  #Residuals:
  #  Min       1Q   Median       3Q      Max 
  #-0.08443 -0.00641 -0.00201  0.00517  0.42136 
  #
  #Coefficients:
  #  Estimate Std. Error t value Pr(>|t|)    
  #(Intercept)  2.365e-02  1.012e-05    2338   <2e-16 ***
  #  P.value     -4.780e-02  1.553e-05   -3079   <2e-16 ***
  #  SE           8.054e-01  2.059e-04    3912   <2e-16 ***
  #  ---
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  #
  #Residual standard error: 0.01407 on 9707618 degrees of freedom
  #Multiple R-squared:  0.7161,	Adjusted R-squared:  0.7161 
  #F-statistic: 1.224e+07 on 2 and 9707618 DF,  p-value: < 2.2e-16
  
  
  summary( lm(Beta ~ P.value + SE , data = GWAS[GWAS$IN_PRSICE==1,]) )
  #Residuals:
  #  Min        1Q    Median        3Q       Max 
  #-0.087588 -0.011568  0.000134  0.011548  0.161559 
  #
  #Coefficients:
  #  Estimate Std. Error t value Pr(>|t|)
  #(Intercept)  5.185e-05  2.262e-04   0.229    0.819
  #P.value     -6.445e-05  2.246e-04  -0.287    0.774
  #SE           2.686e-03  1.542e-02   0.174    0.862
  #
  #Residual standard error: 0.01662 on 70248 degrees of freedom
  #Multiple R-squared:  1.459e-06,	Adjusted R-squared:  -2.701e-05 
  #F-statistic: 0.05123 on 2 and 70248 DF,  p-value: 0.9501
  
  
  summary( lm(abs(Beta) ~ P.value + SE , data = GWAS[GWAS$IN_PRSICE==1,]) )
  #Residuals:
  #  Min        1Q    Median        3Q       Max 
  #-0.018460 -0.002420 -0.001314  0.001288  0.134203 
  #
  #Coefficients:
  #  Estimate Std. Error t value Pr(>|t|)    
  #(Intercept)  1.218e-02  5.488e-05   222.0   <2e-16 ***
  #  P.value     -2.990e-02  5.450e-05  -548.7   <2e-16 ***
  #  SE           9.779e-01  3.742e-03   261.3   <2e-16 ***
  #  ---
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  #
  #Residual standard error: 0.004032 on 70248 degrees of freedom
  #Multiple R-squared:  0.829,	Adjusted R-squared:  0.829 
  #F-statistic: 1.703e+05 on 2 and 70248 DF,  p-value: < 2.2e-16
  
  
  
  ggplot(GWAS[1:10000,], aes(x = -log(P.value), y = SE)) + geom_point(aes(colour = abs(Beta))) + scale_colour_gradient(low="red",high="blue") 
  + geom_point(data=head(GWAS , n = 10000)[ head(GWAS , n = 10000)$MinFreq < 0.05,])
  
  
  
  write(GWAS[GWAS$SE < 0.0087 , "RSNUMBERS"] , "/Users/mohammedjanahi/Desktop/UKB NOMOGRAM PROJECT/ukb_data/HV_GWAS_SNPS_SE_0.0087")
  write(GWAS[GWAS$SE < 0.01 , "RSNUMBERS"] , "/Users/mohammedjanahi/Desktop/UKB NOMOGRAM PROJECT/ukb_data/HV_GWAS_SNPS_SE_0.01")
  
  
}


HARMONIZATION_EXPLORATIONS <- function(){
  
  #### ATTEMPTS AT HARMONIZATION ####
  
  # genetic principle components 
  ukb_gpcs = read.csv("~/OneDrive - University College London/ADNRE STAGE/Datasets/UKBB_predpc.csv")
  adni_gpcs = read.csv("~/OneDrive - University College London/ADNRE STAGE/Datasets/ADNI_predpc.csv")
  
  ukb_all = merge(ukb_all , ukb_gpcs , by="eid")
  adni_all = merge(adni_all , adni_gpcs , by = "PTID")
  
  for( pre in c("","INTERSECT_")){
    # NueroCombat
    ukb_all_nc = ukb_all[!is.na(ukb_all$PRS_TH_1) & !is.na(ukb_all$SEE) , ]
    adni_all_nc = adni_all[!is.na(adni_all$PRS_TH_1) & !is.na(adni_all$DX) & !is.na(adni_all$SEE) & (adni_all$VISCODE.x == "bl") , ]
    
    d = t(as.matrix(rbind(ukb_all_nc[ , paste0(pre,PRS_columns)] , adni_all_nc[ , paste0(pre,PRS_columns)])))
    
    #b = c( rep(0 , nrow(ukb_all_nc)) , adni_all_nc$SITE )
    b = c( rep(0 , nrow(ukb_all_nc)) , rep(1 , nrow(adni_all_nc)) )
    
    ukb_ICVs <- (ukb_all_nc$ICV.x - mean(ukb_all_nc$ICV.x))/sd(ukb_all_nc$ICV.x) 
    adni_ICVs <- (adni_all_nc$ICV - mean(adni_all_nc$ICV))/sd(adni_all_nc$ICV) 
    ICVs <- c(ukb_ICVs , adni_ICVs)
    
    ukb_sex <- ukb_all_nc$Sex == 1
    adni_sex <- adni_all_nc$Sex
    sexs <- as.factor(c(ukb_sex , adni_sex))
    
    ages <- c(ukb_all_nc$AGE_Latest , adni_all_nc$TRUE.AGE)
    
    ukb_DX <- rep("1" , nrow(ukb_all_nc))
    adni_DX <- adni_all_nc$DX
    diagnosis <- as.factor(c(ukb_DX , adni_DX))
    
    gpc_1 <- c(ukb_all_nc$NWE , adni_all_nc$NWE)
    
    gpc_2 <- c(ukb_all_nc$SEE , adni_all_nc$SEE)
    
    covs <- model.matrix(~ICVs + sexs + ages + diagnosis + gpc_1 + gpc_2)
    nc = neuroCombat(dat=d , batch = b , ref.batch = 0 , mod = covs)
    
    i= 1
    for( col in PRS_columns){
      ukb_all_nc[ , paste0(pre,col,"_NC") ] = (nc$dat.combat[i, (1:length(ukb_ICVs)) ])
      adni_all_nc[, paste0(pre,col,"_NC") ] = (nc$dat.combat[i, -(1:length(ukb_ICVs)) ])
      
      ukb_all[ ukb_all$eid %in% ukb_all_nc$eid , paste0(pre,col,"_NC") ] = (nc$dat.combat[i, (1:length(ukb_ICVs)) ])
      adni_all[ paste0(adni_all$PTID , adni_all$VISCODE.x) %in% paste0(adni_all_nc$PTID , adni_all_nc$VISCODE.x)  , paste0(pre,col,"_NC") ] = (nc$dat.combat[i, -(1:length(ukb_ICVs)) ])
      i = i + 1
    }
    
    for(neuro in c("","_NC")){
      for( col in PRS_columns){
        ukb_all[  ukb_all$eid %in% ukb_all_nc$eid , paste0(pre,col,neuro,"_NO_GPC") ] = REGRESS_OUT_COLUMNS(ukb_all_nc , paste0(pre,col,neuro) , c("NWE" , "SEE"))
        adni_all[ paste0(adni_all$PTID , adni_all$VISCODE.x) %in% paste0(adni_all_nc$PTID , adni_all_nc$VISCODE.x), paste0(pre,col,neuro,"_NO_GPC") ] = REGRESS_OUT_COLUMNS(adni_all_nc , paste0(pre,col,neuro) , c("NWE" , "SEE"))
      }
    }
  }
  
  # Check if harmonization worked. 
  
  p.vals = data.frame(array(dim = c(0,9)))
  par(mfrow=c(2,5))
  for( col in PRS_columns){
    #x = hist(ukb_all[ , col]  , plot = FALSE)
    #y = hist(adni_all[ , col] , plot = FALSE)
    #z = hist(ucl_all[ , col]  , plot = FALSE)
    
    #max_y = max(x$density , y$density , z$density)
    #plot( x$mids , x$density , type="l", ylim = c(0,max_y) , col="blue" , main = col)
    #lines(y$mids , y$density , col="red")
    #lines(z$mids , z$density , col="orange")
    ukb_tmp = ukb_all # [ !is.na(ukb_all$SEE) , ]
    adni_tmp = adni_all # [ !is.na(adni_all$SEE)& !is.na(adni_all$DX) , ]
    p.vals = rbind(p.vals , c( col , 
                               #t.test(ukb_tmp[ , col], adni_tmp[ , col])$p.value , 
                               #t.test(ukb_tmp[ , paste0(col,"_NO_GPC")], adni_tmp[, paste0(col,"_NO_GPC")])$p.value ,
                               #t.test(ukb_tmp[ , paste0(col,"_NC")], adni_tmp[, paste0(col,"_NC")])$p.value,
                               #t.test(ukb_tmp[ , paste0(col,"_NC_NO_GPC")], adni_tmp[, paste0(col,"_NC_NO_GPC")])$p.value , 
                               
                               #t.test(ukb_tmp[ , paste0("INTERSECT_",col)], adni_tmp[ , paste0("INTERSECT_",col)])$p.value,
                               #t.test(ukb_tmp[ , paste0("INTERSECT_",col,"_NO_GPC")], adni_tmp[, paste0("INTERSECT_",col,"_NO_GPC")])$p.value ,
                               #t.test(ukb_tmp[ , paste0("INTERSECT_",col,"_NC")], adni_tmp[, paste0("INTERSECT_",col,"_NC")])$p.value ,
                               #t.test(ukb_tmp[ , paste0("INTERSECT_",col,"_NC_NO_GPC")], adni_tmp[, paste0("INTERSECT_",col,"_NC_NO_GPC")])$p.value 
                               
                               cor(ukb_tmp[ , "clean_hv_bilateral"], ukb_tmp[ , col], use = "complete.obs") , 
                               cor(ukb_tmp[ , "clean_hv_bilateral"], ukb_tmp[, paste0(col,"_NO_GPC")] , use = "complete.obs") ,
                               cor(ukb_tmp[ , "clean_hv_bilateral"], ukb_tmp[, paste0(col,"_NC")], use = "complete.obs"),
                               cor(ukb_tmp[ , "clean_hv_bilateral"], ukb_tmp[, paste0(col,"_NC_NO_GPC")], use = "complete.obs") , 
                               
                               cor(ukb_tmp[ , "clean_hv_bilateral"], ukb_tmp[ , paste0("INTERSECT_",col)], use = "complete.obs"),
                               cor(ukb_tmp[ , "clean_hv_bilateral"], ukb_tmp[, paste0("INTERSECT_",col,"_NO_GPC")], use = "complete.obs") ,
                               cor(ukb_tmp[ , "clean_hv_bilateral"], ukb_tmp[, paste0("INTERSECT_",col,"_NC")], use = "complete.obs") ,
                               cor(ukb_tmp[ , "clean_hv_bilateral"], ukb_tmp[, paste0("INTERSECT_",col,"_NC_NO_GPC")], use = "complete.obs") 
    )) #,
    #                            t.test(ukb_all[ , col], ucl_all[  , col])$p.value,
    #                            t.test(adni_all[ , col], ucl_all[ , col])$p.value ) )
    
    
  }
  names(p.vals) = c( "thresh" , "raw_prs" , "raw_gpc" ,"raw_nc" , "raw_nc_gpc" , "int_prs" ,"int_gpc" , "int_nc" , "int_nc_gpc")# , "adni/ucl")
  print(p.vals)
  
  par(mfrow=c(4,5))
  for(pre in c("","INTERSECT_") )
    for(nc in c("","_NC"))
      for(gpc in c("","_NO_GPC") )
        for( col in PRS_columns[c(1,6,8,11,14)]){
          x = hist( ukb_all[ , paste0(pre,col,nc,gpc) ] , plot = FALSE)
          y = hist(adni_all[ , paste0(pre,col,nc,gpc) ] , plot = FALSE)
          
          max_y = max(x$counts , y$counts )
          max_x = max(abs(x$breaks) , abs(y$breaks))
          plot( x$mids , x$counts , type="l", ylim = c(0,max_y) , xlim=c(-max_x , max_x), col="blue" , main = paste0(pre,col,nc,gpc))
          lines(y$mids , y$counts , col="red")
        }
  
  # here I do some experiments to check PGS's generated from intersecting sets of snps and how that changes the raw scores
  p.vals = data.frame(array(dim = c(0,7)))
  for( col in PRS_columns){
    p.vals = rbind(p.vals , c(col, 
                              t.test( ukb_all[,col] , ukb_all[,paste0("ADNI_INTERSECT_",col)])$p.value ,
                              abs(mean(ukb_all[,col] , na.rm=TRUE) - mean(ukb_all[,paste0("ADNI_INTERSECT_",col)] , na.rm=TRUE)) / 
                                t.test( ukb_all[,col] , ukb_all[,paste0("I46_INTERSECT_",col)])$p.value ,
                              t.test( ukb_all[,col] , ukb_all[,paste0("I46_ADNI_INTERSECT_",col)])$p.value ,
                              
                              t.test( ukb_all[,paste0("ADNI_INTERSECT_",col)] , ukb_all[,paste0("I46_INTERSECT_",col)])$p.value ,
                              t.test( ukb_all[,paste0("ADNI_INTERSECT_",col)] , ukb_all[,paste0("I46_ADNI_INTERSECT_",col)])$p.value ,
                              
                              t.test( ukb_all[,paste0("I46_INTERSECT_",col)] , ukb_all[,paste0("I46_ADNI_INTERSECT_",col)])$p.value
                              
    ))
  }
  
  names(p.vals) = c( "thresh" , "raw vs adni" , "raw vs i46" , "raw vs all" , "adni vs i46" ,"adni vs all" , "i46 vs all" )
  print(p.vals)
  
  
  ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , "LASSO"] = LASSO_SCORE(ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , ] , ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , "clean_hv_bilateral"] , PRS_columns )
  ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , "ADNI_INTERSECT_LASSO"] = LASSO_SCORE(ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , ] , ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , "clean_hv_bilateral"] , paste0("ADNI_INTERSECT_",PRS_columns) )
  ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , "I46_INTERSECT_LASSO"] = LASSO_SCORE(ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1),] , ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , "clean_hv_bilateral"] , paste0("I46_INTERSECT_",PRS_columns) )
  ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , "I46_ADNI_INTERSECT_LASSO"] = LASSO_SCORE(ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1),] , ukb_all[!is.na(ukb_all$clean_hv_bilateral) & !is.na(ukb_all$PRS_TH_1) , "clean_hv_bilateral"] , paste0("I46_ADNI_INTERSECT_",PRS_columns) )
  
  par(mfrow=c(4,4))
  for( col in PRS_columns){
    
    
    a = hist(ukb_all[ , col ] , plot = FALSE)
    b = hist(ukb_all[ , paste0("ADNI_INTERSECT_",col) ] , plot = FALSE)
    c = hist(ukb_all[ , paste0("I46_INTERSECT_",col) ] , plot = FALSE)
    d = hist(ukb_all[ , paste0("I46_ADNI_INTERSECT_",col) ] , plot = FALSE)
    
    max_y = max(a$counts , b$counts , c$counts , d$counts)
    max_x = max(abs(a$breaks) , abs(b$breaks) , abs(c$breaks) , abs(d$breaks))
    plot( a$mids , a$counts , type="l", ylim = c(0,max_y) , xlim=c(-max_x , max_x), col="black" , main = col )
    lines(b$mids , b$counts , col="red")
    lines(c$mids , c$counts , col="blue")
    lines(d$mids , d$counts , col="orange")
  }
  
  #### END
}


PLOT_OUTLIERS <- function( table , ICD10_columns){
  # training set 
  group_0 = table$training_set
  # A80-A89 Viral infections of the central nervous system
  group_1 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^A8[0-9]"))) 
  # C69-C72 Malignant neoplasms of eye, brain and other parts of central nervous system
  group_2 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^C(69|70|71|72)")))
  # F10-F19 Mental and behavioural disorders due to psychoactive substance use
  group_3 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F1[0-9]"))) 
  # F20-F29 Schizophrenia, schizotypal and delusional disorders
  group_4 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F2[0-9]"))) 
  # F30-F39 Mood [affective] disorders
  group_5 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F3[0-9]"))) 
  # F40-F48 Neurotic, stress-related and somatoform disorders
  group_6 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F4[0-9]"))) 
  # F50-F59 Behavioural syndromes associated with physiological disturbances and physical factors
  group_7 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F5[0-9]"))) 
  # F60-F69 Disorders of adult personality and behaviour
  group_8 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F6[0-9]"))) 
  # F70-F79 Mental retardation
  group_9 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F7[0-9]"))) 
  # F80-F89 Disorders of psychological development
  group_10 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F8[0-9]"))) 
  # F90-F98 Behavioural and emotional disorders with onset usually occurring in childhood and adolescence
  group_11 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F9[0-8]"))) 
  # F99-F99 Unspecified mental disorder
  group_12 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F99"))) 
  # G00-G09 Inflammatory diseases of the central nervous system
  group_13 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^G0[0-9]"))) 
  # G10-G14 Systemic atrophies primarily affecting the central nervous system
  group_14 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^G1[0-4]"))) 
  # G35-G37 Demyelinating diseases of the central nervous system
  group_15 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^G3[5-7]"))) 
  # I60-I69 Cerebrovascular diseases
  group_16 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^I6[0-9]"))) 
  # F00 Dementia in Alzheimer's disease
  group_17 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F00"))) 
  # F01 Vascular dementia
  group_18 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F01"))) 
  # F02 Dementia in other diseases classified elsewhere
  group_19 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F02"))) 
  # F03 Unspecified dementia
  group_20 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F03"))) 
  # F06 Other mental disorders due to brain damage and dysfunction and to physical disease
  group_21 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F06"))) 
  # S00-S09 Injuries to the head
  group_22 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^S0[0-9]"))) 

  group_23 = tibble(table$test_set ) %>% filter_at(vars(ICD10_columns), any_vars(stringi::stri_detect_regex(., "^F0[0-3]"))) 
  
  
  group_names = c("TRAINING SET",
                  "A80-A89 Viral infections of the central nervous system",
                  "C69-C72 Malignant neoplasms of eye, brain and other parts of central nervous system",
                  "F10-F19 Mental and behavioural disorders due to psychoactive substance use",
                  "F20-F29 Schizophrenia, schizotypal and delusional disorders",
                  "F30-F39 Mood [affective] disorders",
                  "F40-F48 Neurotic, stress-related and somatoform disorders",
                  "F50-F59 Behavioural syndromes associated with physiological disturbances and physical factors",
                  "F60-F69 Disorders of adult personality and behaviour",
                  "F70-F79 Mental retardation",
                  "F80-F89 Disorders of psychological development",
                  "F90-F98 Behavioural and emotional disorders with onset usually occurring in childhood and adolescence",
                  "F99-F99 Unspecified mental disorder",
                  "G00-G09 Inflammatory diseases of the central nervous system",
                  "G10-G14 Systemic atrophies primarily affecting the central nervous system",
                  "G35-G37 Demyelinating diseases of the central nervous system",
                  "I60-I69 Cerebrovascular diseases",
                  "F00 Dementia in Alzheimer's disease",
                  "F01 Vascular dementia",
                  "F02 Dementia in other diseases classified elsewhere",
                  "F03 Unspecified dementia",
                  "F06 Other mental disorders due to brain damage and dysfunction and to physical disease",
                  "S00-S09 Injuries to the head",
                  "F Dementia Group")
  
  TH = 9
  all_groups_pass_threshold = data_frame(matrix(ncol = 2, nrow = 0))
  colnames(all_groups_pass_threshold) <- c("zscore" , "groupIndex")
  for( i in 0:23){
    group <- get(paste( "group_" , i, sep=""))
    group_df <- tibble(group$zscores , i)
    names(group_df) <- c("zscore" , "groupIndex")
    #print(paste( nrow( group ) , group_names[i+1]))
    if (nrow( group ) > TH ){
      all_groups_pass_threshold <- rbind(all_groups_pass_threshold , group_df )
    }
  }
  
  boxplot( all_groups_pass_threshold$zscore ~ all_groups_pass_threshold$groupIndex , 
           ylab = "Z-Scores" , xlab = "" , names = 1:length(unique(all_groups_pass_threshold$groupIndex)) )
  print(group_names[ unique(all_groups_pass_threshold$groupIndex) + 1])
  
  abline(h = quantile(table$training_set$zscores , 0.75))
  abline(h = quantile(table$training_set$zscores , 0.50))
  abline(h = quantile(table$training_set$zscores , 0.25))
}


SUBSAMPLE_DIFF <- function( table , order_cols , eval_col , seed = NA){
  
  n = nrow(table)
  
  if(! is.na(seed)) set.seed(seed)
  
  train_rows <- sample(1:n , .5*n )
  test_set <- table[ -train_rows, ]
  
  if(length(order_cols) > 1) 
    ord = order(rowSums(test_set[,order_cols]) )
  else
    ord = order(test_set[,order_cols])
  
  test_set = test_set[ ord  , ]
  
  table_high <- tail( test_set , .3*n )
  table_low <- head( test_set , .3*n )
  
  diff <- (mean(table_high[ , eval_col],na.rm=TRUE)) - (mean(table_low[ , eval_col],na.rm=TRUE))
  
  return(diff)
  
}



SUBSAMPLE_DIFF_NOMOGRAM <- function( table , eval_col , age_col , id_col , seed , mean_res = TRUE , plot = FALSE){
  
  set.seed(seed)
  train_rows <- sample(1:n, .5*n)
  XX <- GPR_ANALYSIS_SEPARATED( table[-train_rows ,] , eval_col , c(0.3) , P = FALSE , age_column = age_col  , id_column = id_col ,
                                XX = seq(range(table[, age_col])[1], range(table[, age_col])[2] , 0.1) )
  
  diff <- NOMOGRAM_DIFFERENCE(XX$upper , XX$lower , mean_res = mean_res)
  #print(paste(eval_col,diff))
  if(plot){
    PLOT_NOMOGRAM_COMPARE(XX$upper , XX$lower , shade=FALSE , xlim = c(40,90) , ylim=c( 1000 , 6000) , title = paste(eval_col,diff) )
    table_2 = table[-train_rows,]
    table <- table_2[ order(table_2[,eval_col]) , ]
    upper_eids <- tail( table_2[,id_col] , 0.3 * nrow(table_2) )
    lower_eids <- head( table_2[,id_col] , 0.3 * nrow(table_2) )
    
    abline(v = min(max(table_2[ table_2[,id_col] %in% upper_eids, age_col]) , max( table_2[ table_2[,id_col] %in% lower_eids, age_col] ) )  ) 
    abline(v = max(min(table_2[ table_2[,id_col] %in% upper_eids, age_col]) , min( table_2[ table_2[,id_col] %in% lower_eids, age_col] ) )  ) 
    
    points( table_2[table_2[,id_col] %in% upper_eids , age_col] , table_2[table_2[,id_col] %in% upper_eids , "clean_hv_bilateral"] , col = "black")
    points( table_2[table_2[,id_col] %in% lower_eids , age_col] , table_2[table_2[,id_col] %in% lower_eids , "clean_hv_bilateral"] , col = "red")
    
    #print(paste("outliers at age:" , boxplot.stats(table[-train_rows , age_col])$out))
}
  
  return(diff)
}


CORRECT_PRS <- function( table , column , gpcs){
  
  table[!is.na(table[,column]), paste0(column,"_NO_GPC")] <- REGRESS_OUT_COLUMNS(table , column , gpcs)
  
  col_mean <- mean(table[,paste0(column,"_NO_GPC")] , na.rm = TRUE)
  col_sds <- sd(table[,paste0(column,"_NO_GPC")] , na.rm = TRUE)
  z_scores <- ( table[,paste0(column,"_NO_GPC")] - col_mean) / col_sds
  
  #table[ , paste0(column,"_NO_GPC","_zscore")] =  sweep( sweep(table[ , paste0(column,"_NO_GPC")]  , 2 , means) , 2 , sds , "/")
  
  table[, paste(column,"_NO_GPC_ZSCORE",sep = "")] <- z_scores
  
  return( table)
}



TEST <- function( PRS_columns , table , AGE_COL , ID_COL , SEEDS ){
  diffs_1 = c()
  diffs_2 = c()
  
  pairs_1 = combn( PRS_columns , 2 , simplify = FALSE)
  pairs_2 = combn(paste0(PRS_columns , "_ZSCORE") , 2 , simplify = FALSE)
  for (i in 1:length(pairs_1) ) {
    #for( prs in PRS_columns ){
    PRS_1 = pairs_1[[i]][1] 
    PRS_2 = pairs_1[[i]][2] 
    
    score_col_name = paste0(PRS_1 ,"&",PRS_2)
    score_col_name_2 = paste0(PRS_1 ,"&",PRS_2 ,"_ZSCORE")
    
    table[ , score_col_name ] = table[, PRS_1 ] + table[, PRS_2 ]
    table[ , score_col_name_2 ] = table[, pairs_2[[i]][1] ] + table[, pairs_2[[i]][2] ] 
    
    for (seed in SEEDS) {
      #diff = SUBSAMPLE_DIFF_NOMOGRAM(table , c(score_col_name) , AGE_COL , ID_COL , seed , plot = FALSE)
      diff_2 = SUBSAMPLE_DIFF_NOMOGRAM(table , c(score_col_name_2 ) , AGE_COL , ID_COL , seed , plot = TRUE)
      
      #diffs_1 = c(diff , diffs_1)
      #diffs_2 = c(diff_2 , diffs_2)
      
      #print(paste(PRS_1 ,"&", PRS_2 ,"seed:",seed,"normal diff:" , diff , "Zscore diff:",diff_2))
      #print(paste(prs ,"seed:",seed,"normal diff:" , diff , "Zscore diff:",diff_2))
    }
  }
  
  print(mean(diffs_1))
  print(mean(diffs_2))
}




####### OLD TWO PGS METHOD ######
# all_diffs = data.frame(matrix(ncol = 5, nrow = 0))
# iter_ind = 1
# for( sex in c("male","female") ){
#   print(paste("######",sex,"######"))
#   table = get(paste(sex,"table",sep="_"))
#   
#   # We will find which scores have upper and lower strata 
#   inter_table = get(paste(sex,"intersections_table",sep="_"))
#   inter_table = (inter_table > .5) & (inter_table < .8)
#   seq_AD_cols =  1:(length(PRS_AD_columns))
#   tt = inter_table[ seq_AD_cols, seq_AD_cols] & inter_table[-seq_AD_cols , -seq_AD_cols]
#   tt = lower.tri(tt)*tt
#   wh = which( tt == TRUE, arr.in= TRUE)
#   x = cbind.data.frame( colIDs = colnames(tt)[ wh[, 1] ], rowIDs = rownames(tt)[ wh[, 2] ])
#   
#   for(row in 1:nrow(x)){
#     upper_1 = x[row,1]
#     upper_2 = x[row,2]
#     
#     lower_1 = gsub( "upper" , "lower" , upper_1 )
#     lower_2 = gsub( "upper" , "lower" , upper_2 )
#     
#     upper_intersection = intersect(get(upper_1)[,ID_COL] , get(upper_2)[,ID_COL])
#     lower_intersection = intersect(get(lower_1)[,ID_COL] , get(lower_2)[,ID_COL])
#     
#     lower_limit = floor( nrow(table) * .5 * .3 )
#     
#     for(seed in SEEDS){
#       set.seed(seed)
#       upper_rows <- sample(upper_intersection, lower_limit)
#       lower_rows <- sample(lower_intersection, lower_limit)
#       
#       diff <- (mean(table[ table[ ,ID_COL] %in% upper_rows, "familial_AD_score"],na.rm=TRUE)) - (mean(table[ table[ ,ID_COL] %in% lower_rows , "familial_AD_score"],na.rm=TRUE))
#       
#      all_diffs <- rbind(all_diffs , list(iter_ind , sex, gsub("_upper_" , "", gsub(sex , "" , upper_1)) ,gsub("_upper_" , "", gsub(sex , "" , upper_2)) , diff) )
#    }
#    iter_ind = iter_ind + 1
#  }
#}
#
#all_diffs = setNames( all_diffs , c("index" , "sex", "prs1", "prs2" , "diff") )
#
#for( x in unique(all_diffs$index))
#  print(paste( all_diffs[x*length(SEEDS) , "sex"] , 
#               all_diffs[x*length(SEEDS) , "prs1"] , 
#               all_diffs[x*length(SEEDS) , "prs2"] , 
#               mean(all_diffs[ all_diffs$index == x , "diff" ]), 
#               sd(all_diffs[ all_diffs$index == x , "diff" ])
#  ))
##




##### OLD ATTEMPTS AT EVALUATING ZSCORES ####
#
# Find the samples that cross the 90% z-score threshold from one score to another
# z_th = qnorm(.9)
# 
# more_extreme = adni_all[ abs(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean) < z_th  & abs(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean) > z_th , ]
# more_normal = adni_all[ abs(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean) > z_th  & abs(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean) < z_th , ]
# extreme_extreme = adni_all[ abs(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean) > z_th  & abs(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean) > z_th , ]
# normal_normal = adni_all[ abs(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean) < z_th  & abs(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean)  < z_th , ]
# 
# print(paste("more_extreme:",sum(!is.na(more_extreme$EcogPtMem)) , "more_normal:",sum(!is.na(more_normal$EcogPtMem)) , 
#             "extreme_extreme:",sum(!is.na(extreme_extreme$EcogPtMem)) , "normal_normal:",sum(!is.na(normal_normal$EcogPtMem))))
# 
# par(mfrow = c(4,4))
# for( table in list(more_extreme , more_normal , extreme_extreme , normal_normal ) ){
#   for( col in eval_cols_zscore[1:4] ){
#     my_boxplot(table , grep(paste0(col,"_model_") , names(table) , value = TRUE) , ylab = col)
#     abline(h=1)
#   }
# }
# 
# # Find the samples that cross the shift in terms of z-score by more than a set percentile amount (20%)
# adni_all$non_genetic_model_mean_zscore = rowMeans(adni_all[, c("clean_hv_bilateral_nc_gpc_zscore_1_mean","clean_hv_bilateral_nc_gpc_zscore_2_mean","clean_hv_bilateral_nc_gpc_zscore_3_mean")])
# adni_all$genetic_model_mean_zscore = rowMeans(adni_all[, c("clean_hv_bilateral_nc_gpc_zscore_4_mean","clean_hv_bilateral_nc_gpc_zscore_5_mean","clean_hv_bilateral_nc_gpc_zscore_6_mean","clean_hv_bilateral_nc_gpc_zscore_7_mean")])
# adni_all$single_pgs_model_mean_zscore = rowMeans(adni_all[, c("clean_hv_bilateral_nc_gpc_zscore_4_mean","clean_hv_bilateral_nc_gpc_zscore_5_mean","clean_hv_bilateral_nc_gpc_zscore_6_mean")])
# 
# 
# adni_all$clean_hv_bilateral_nc_gpc_zscore_3_to_7_shift_percentile =  
#   pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean) - pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean)
# 
# shift = 0.20
# 
# shifted = adni_all[ abs( pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean) - pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean))  > shift , ]
# 
# shifted_in = adni_all[ (pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean)) - (pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean))  > shift , ]
# shifted_out = adni_all[ (pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean)) - (pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean))  > shift , ]
# 
# 
# non_shifted = adni_all[ abs( pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean) - pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean))  < shift , ]
# 
# print(paste("shifters:",sum(!is.na(shifted$EcogPtMem)) , "non-shifters:",sum(!is.na(non_shifted$EcogPtMem))))
# 
# 
# par(mfrow = c(3,4))
# for( table in list(non_shifted , shifted_in , shifted_out) ){
#   for( col in eval_cols_zscore[4:7] ){
#     my_boxplot(table , grep(paste0(col,"_model_") , names(table) , value = TRUE)[c(3,7)],
#                ylab = paste(strsplit(col , "_")[[1]][1]," Z-score") , xlab = "HV Z-score Model" , xaxt = c(1:7))
#     abline(h=1)
#   }
# }
# 
# par(mfrow=c(2,2))
# for( table in list(shifted_in , shifted_out) ){
#   # Melt the data frame into a longer format suitable for plotting
#   mean_cols = grep( "_mean" , names(table) , value = TRUE)
#   table_tmp = table[ complete.cases(table[, mean_cols]) , ]
#   melted_data <- reshape2::melt(table_tmp[ , c("id" , mean_cols ) ] )
#   # Plot using ggplot2
#   p = ggplot(melted_data, aes(x = variable, y = value, group = factor(id) )) +
#     geom_line() +
#     labs(x = "Model", y = "Z-Score", title = "Z-Score Changes Across Models")
#   plot(p)
#   
# }
# 
# 
# par(mfrow=c(1,1))
# plot( (pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean) - pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean)) ,
#       adni_all$EcogPtMem_zscore ,
#       xlab = "amount of shift" , ylab = "ecog mem")
# 
# plot( abs(pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_3_mean) - pnorm(adni_all$clean_hv_bilateral_nc_gpc_zscore_7_mean)) ,
#       adni_all$EcogPtMem_zscore_model_7 ,
#       xlab = "amount of shift" , ylab = "ecog mem")
# 
# 
# 
# 
# par(mfrow = c(2,3))
# for( table in list(adni_all) ){
#   for( col in eval_cols_zscore[4:8] ){
#     plot( adni_all$clean_hv_bilateral_nc_gpc_zscore_3_to_7_shift_percentile , adni_all[ , col] ,
#           ylab = paste(strsplit(col , "_")[[1]][1]," Z-score") , xlab = "Percentile shift between model 3 & 7" )
#   }
# }
#
####