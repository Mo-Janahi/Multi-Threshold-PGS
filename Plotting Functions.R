# File Containing all functions that deal with Plotting Nomograms
# By: Mohammed Janahi
#
#


# Basic Nomogram Plotting function
PLOT_NOMOGRAM <- function( bins , hem = "bilateral" , title = " " , ylim = c(2700,5300) , xlim = c(40,85) , scatter = NA ){
  
  centers <- bins$mid_age
  
  suf <- ifelse( hem == "bilateral" , "" , paste("_",hem,sep="") )
  
  mn <- bins[ , paste("mean",suf,sep="") ]
  q2.5 <- bins[ , paste("q2.5",suf,sep="")]
  q5 <- bins[ , paste("q5",suf,sep="")]
  q10 <- bins[ , paste("q10",suf,sep="")]
  q25 <- bins[ , paste("q25",suf,sep="")]
  q50 <- bins[ , paste("q50",suf,sep="")]
  q75 <- bins[ , paste("q75",suf,sep="")]
  q90 <- bins[ , paste("q90",suf,sep="")]
  q95 <- bins[ , paste("q95",suf,sep="")]
  q97.5 <- bins[ , paste("q97.5",suf,sep="")]
  
  
  if( is.na(xlim[1])){
    min_x <- floor(head(centers , 1 )) - 2
    max_x <- floor(tail(centers , 1 )) + 2
  }
  else{ 
    min_x <- xlim[1]
    max_x <- xlim[2]
  }
  if( is.na(ylim[1])){
    min_y <- floor(min(q2.5[(centers < max_x) & (centers > min_x)]  )) - 100
    max_y <- floor(max(q97.5[(centers < max_x) & (centers > min_x)] )) + 100
  }
  else{
    min_y <- ylim[1]
    max_y <- ylim[2]
  }
  # artificial max for testing. comment out later
  #max_x <- 90
  
  plot( main = title , xlab = "Age" , ylab = "Hippocampal Volume" ,
        centers , mn,
        ylim=c(min_y , max_y), xlim=c(  min_x , max_x ) , 
        xaxt="n" , yaxt="none", type='l')
  
  axis( 1, seq(min_x , max_x , 1) )
  axis( 2, seq(min_y , max_y,100) , las=2)
  
  abline(v=seq(min_x, max_x, 5) , col='grey' , lty=2)
  abline(h=seq(min_y , max_y, 200) , col='grey' , lty=2)
  
  if( !is.na(scatter[1]))
    points( scatter$AGE_Latest, scatter[,paste("clean_hv_",hem,sep="")] , col="grey")
  
  polygon( c( centers , rev(centers) ) , c(q2.5 , rev(q97.5)) , col=rgb(1,0.61,0.16,0.2) , border = NA)
  polygon( c( centers , rev(centers) ) , c(q5 , rev(q95)) , col=rgb(1,0.61,0.16,0.2) , border = NA)
  polygon( c( centers , rev(centers) ) , c(q10 , rev(q90)) , col=rgb(1,0.61,0.16,0.2) , border = NA)
  polygon( c( centers , rev(centers) ) , c(q25 , rev(q75)) , col=rgb(1,0.61,0.16,0.2) , border = NA)
  
  lines( centers , q2.5 )   
  lines( centers , q5)
  lines( centers , q10)
  lines( centers , q25)
  lines( centers , mn)
  lines( centers , q75)
  lines( centers , q90)
  lines( centers , q95)
  lines( centers , q97.5)
}

# Plot two nomograms over each other with an option to shade the distance between them
PLOT_NOMOGRAM_COMPARE <- function( bins_1 , bins_2 , hem="bilateral" , col_1="black" , col_2="red" , title = " " , ylim = c(2700,5300) , xlim = c(40,85) , shade=TRUE){
  
  suf <- ifelse( hem == "bilateral" , "" , paste("_",hem,sep="") )
  percentiles <- paste("q",LEVELS , suf , sep = "")
  
  
  if( is.na(xlim[1])){
    min_x <- floor( min( head(bins_1$mid_age,1) , head(bins_2$mid_age,1) )) - 2
    max_x <- floor( max( tail(bins_1$mid_age,1) , tail(bins_2$mid_age,1) )) + 2
  }
  else{ 
    min_x <- xlim[1]
    max_x <- xlim[2]
  }
  if( is.na(ylim[1])){
    # we want the y scale to be outside the most extreme values between the bins being plotted
    # we want the y scale to be calculated within the x limit given above
    bins_1_in_range <- bins_1[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , ]
    bins_2_in_range <- bins_2[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , ]
    min_y <- floor(min( c(bins_1_in_range[ , percentiles[1] ] , bins_2_in_range[ , percentiles[1] ]) )) - 100
    max_y <- floor(max( c(bins_1_in_range[ , percentiles[length(LEVELS)] ] , bins_2_in_range[ , percentiles[length(LEVELS)] ] ) )) + 100
  }
  else{
    min_y <- ylim[1]
    max_y <- ylim[2]
  }
  
  plot( main = title , bins_1$mid_age , bins_1[, percentiles[floor(length(LEVELS)/2)] ], xlab = "Age" , ylab = "Hippocampal Volume" ,
        type='l' , col=col_1 , lwd=2 , ylim=c(min_y , max_y), xlim=c(  min_x , max_x ) , xaxt="n" , yaxt="none")
  
  axis( 1, seq(min_x , max_x , 1) )
  axis( 2, seq(min_y , max_y,100) , las=2)
  
  abline(v=seq(min_x, max_x, 5) , col='grey' , lty=2)
  abline(h=seq(min_y , max_y, 200) , col='grey' , lty=2)
  
  for ( i in percentiles ){
    
    lines(bins_1$mid_age , bins_1[,i] , col=col_1 , lwd=2 ) 
    lines(bins_2$mid_age , bins_2[,i] , col=col_2 , lwd=2 ) 
    if(shade){
      polygon( c( bins_1$mid_age  ,   rev(bins_2$mid_age)  ), 
               c(bins_1[,i], rev(bins_2[,i]) ), 
               col=rgb(0,0,0,0.5) , border = NA , density = 40 )
    }
  }
  
}

# Plot three nomograms against each other with the assumption that one will be in the middle and one higher and one lower, 
# Shade accordingly
PLOT_NOMOGRAM_COMPARE_BIDIRECTIONAL <- function( bins_1 , bins_2 , bins_3 , hem="bilateral" , col_1 , col_2 , col_3, title, ylim = c(2700,5300) , xlim = c(40,85) , shade=TRUE){
  
  suf <- ifelse( hem == "bilateral" , "" , paste("_",hem,sep="") )
  percentiles <- paste("q", LEVELS , suf , sep = "")
  
  if( is.na(xlim[1])){
    min_x <- floor( min( head(bins_1$mid_age,1) , head(bins_2$mid_age,1) , head(bins_3$mid_age,1) )) - 2
    max_x <- floor( max( tail(bins_1$mid_age,1) , tail(bins_2$mid_age,1) , tail(bins_3$mid_age,1) )) + 2
  }
  else{ 
    min_x <- xlim[1]
    max_x <- xlim[2]
  }
  if( is.na(ylim[1])){
    bins_1_in_range <- bins_1[ (bins_1$mid_age < max_x) & (bins_1$mid_age > min_x) , ]
    bins_2_in_range <- bins_2[ (bins_2$mid_age < max_x) & (bins_2$mid_age > min_x) , ]
    bins_3_in_range <- bins_3[ (bins_3$mid_age < max_x) & (bins_3$mid_age > min_x) , ]
    min_y <- floor(min( c(bins_1_in_range[ , percentiles[1] ] ,
                          bins_2_in_range[ , percentiles[1] ] ,
                          bins_3_in_range[ , percentiles[1] ]) )) - 100
    max_y <- floor(max( c(bins_1_in_range[ , percentiles[length(percentiles)] ] ,
                          bins_2_in_range[ , percentiles[length(percentiles)] ] ,
                          bins_3_in_range[ , percentiles[length(percentiles)] ]) )) + 100
  }
  else{
    min_y <- ylim[1]
    max_y <- ylim[2]
  }
  
  plot( main = title , bins_1$mid_age , bins_1[, paste("q50",suf,sep="") ], xlab = "Age" , ylab = "HV" ,
        type='l' , col=col_1 , lwd=2 , ylim=c(min_y , max_y), xlim=c(  min_x , max_x ) , xaxt="n" , yaxt="none")
  
  axis( 1, seq(min_x , max_x , 1) )
  axis( 2, seq(min_y , max_y,100) , las=2)
  
  abline(v=seq(min_x, max_x, 1) , col='grey' , lty=2)
  abline(h=seq(min_y , max_y, 100) , col='grey' , lty=2)
  
  for ( i in percentiles ){
    
    lines(bins_1$mid_age , bins_1[,i] , col=col_1 , lwd=2 ) 
    lines(bins_2$mid_age , bins_2[,i] , col=col_2 , lwd=2 ) 
    lines(bins_3$mid_age , bins_3[,i] , col=col_3 , lwd=2 ) 
    if(shade){
      polygon( c( bins_1$mid_age  ,   rev(bins_1$mid_age)  ), 
               c(bins_1[,i], rev(bins_2[,i]) ), 
               col=rgb(0,0,0,0.5) , border = NA , density = 40 )
      
      polygon( c( bins_1$mid_age  ,   rev(bins_1$mid_age)  ), 
               c(bins_1[,i], rev(bins_3[,i]) ), 
               col=rgb(0,0,0,0.5) , border = NA , density = 40 )
    }
  }
  
}

# Plot the nomogram just like the PLOT_NOMOGRAM function but add scatter plots of all the different diagnostic groups from ADNI
PLOT_NOMOGRAM_ADNI <- function( bins , adni , hem =  "bilateral" , title = "" , ylim = c(2700,5300) , xlim = c(40,85) ){
  
  adni_AD <- filter(adni , VISCODE.y=="scmri" & DX=="Dementia" )
  adni_MCI <- filter(adni , VISCODE.y=="scmri" & DX=="MCI" )
  adni_CN <- filter(adni , VISCODE.y=="scmri" & DX=="CN" ) 
  adni_NA <- filter(adni , VISCODE.y=="scmri" & is.na(DX) )
  hv_col <- paste( "clean_hv_",hem , sep="")
  
  PLOT_NOMOGRAM( bins , hem = hem , title = title , xlim = xlim , ylim = ylim)
  
  points((adni_AD$TRUE.AGE) , (adni_AD[,hv_col]) , pch=17 , col="red")
  points((adni_MCI$TRUE.AGE) , (adni_MCI[,hv_col]) , pch=16 , col="green")
  points((adni_CN$TRUE.AGE) , (adni_CN[,hv_col]) , pch=15 , col="blue")
  points((adni_NA$TRUE.AGE) , (adni_NA[,hv_col]) , pch=15 , col="grey")
  
}

# Plot the nomogram like we did for PLOT_NOMOGRAM_ADNI but connect all the subsequent visits of the same individual
PLOT_NOMOGRAM_ADNI_LONGITUDINAL <- function(bins , adni , hem="bilateral" , title="" , ylim = c(2700,5300) , xlim = c(40,85) , IDS = NA){
  
  PLOT_NOMOGRAM( bins , hem=hem , title=title, xlim = xlim, ylim=ylim)
  table <- adni[ c("RID.x" , "TRUE.AGE" , "clean_hv_left", "clean_hv_right" , "clean_hv_bilateral" , "DX")]
  hv_col <- paste("clean_hv_", hem, sep="")
  if( is.na(IDS[1]))
    IDS <- unique(table$RID.x)
  DISPLAY_ID = 0
  for (ID in IDS ){
    p_table <- table[ table$RID.x==ID , ]
    p_table <- p_table[ order(p_table$TRUE.AGE) , ]
    #if( length(p_table[ , hv_col]) > 1)
    #if(p_table[1,"TRUE.AGE"] < 73 )
    #if( sum(!is.na(unique(p_table$DX))) > 1 ){      
    if(  !is.na(match("MCI", p_table$DX)) ){
      lines( p_table$TRUE.AGE , p_table[,hv_col] , col=rgb(0,0,0, alpha = 0.5))
      points((p_table[ p_table$DX=="Dementia","TRUE.AGE"]) , (p_table[ p_table$DX=="Dementia",hv_col]) , pch=17 , col="red")
      points((p_table[ p_table$DX=="MCI","TRUE.AGE"]) , (p_table[ p_table$DX=="MCI",hv_col]) , pch=16 , col="green")
      points((p_table[ p_table$DX=="CN","TRUE.AGE"]) , (p_table[ p_table$DX=="CN",hv_col]), pch=15 , col="blue")
      points((p_table[ is.na(p_table$DX),"TRUE.AGE"]) , (p_table[ is.na(p_table$DX),hv_col]) , pch=15 , col="grey")
      set.seed(p_table[1,"RID.x"]) # attempt to annonymize the data labels by using the actual ID as a seed for random number generation
      text( x = p_table[1,"TRUE.AGE"] , y = p_table[1,hv_col] , labels = floor(rnorm(n=1 , mean = 5000 , sd = 1000)) )
    }
  }
  
}



MAKE_PRS_PLOTS <- function( table_male , table_female , confounder_columns , 
                            hv_l_column = "clean_hv_left", hv_r_column = "clean_hv_right",
                            hv_bl_column = "clean_hv_bilateral"){
  
  
  thresholds <- c('PRS_TH_1e.08','PRS_TH_1e.07','PRS_TH_1e.06','PRS_TH_1e.05','PRS_TH_1e.04','PRS_TH_1e.03',
                  'PRS_TH_0.01','PRS_TH_0.05','PRS_TH_0.1','PRS_TH_0.2','PRS_TH_0.4','PRS_TH_0.5','PRS_TH_0.75','PRS_TH_1')
  
  #thresholds <- c('PRS_TH_0.75','PRS_TH_0.75')
  
  #thresholds <- paste("NO_APOE",thresholds,sep="_")
  
  full <- merge(table_male , table_female , by=c("eid",hv_l_column,hv_r_column,hv_bl_column,confounder_columns,thresholds), all.x = TRUE , all.y = TRUE)
  # show what the prs percentiles look like 
  par(mfrow=c(3,6))
  #par(mfrow=c(1,2))
  for( GENDER in c("Male","Female","Both Sexes")){
    #for( GENDER in c("Both Sexes")){
    if (GENDER=="Male")         table <- table_male
    if (GENDER=="Female")       table <- table_female
    if (GENDER=="Both Sexes") table <- full
    for(HEMISPHERE in c("Left","Right","Bilateral")){
      #for(HEMISPHERE in c("Both Hemispheres")){
      if (HEMISPHERE=="Left")             hv <- hv_l_column
      if (HEMISPHERE=="Right")            hv <- hv_r_column
      if (HEMISPHERE=="Bilateral") hv <- hv_bl_column
      
      #table[,hv] <- (table[,hv] - min(table[,hv],na.rm = TRUE)) / (max(table[,hv],na.rm = TRUE) - min(table[,hv],na.rm = TRUE))
      
      HV_SNP <- matrix(c(0,0,0,0,0,0),ncol=6,nrow=length(thresholds))
      
      colnames(HV_SNP) <- c("Slope" , "Range" , "P-Value" , "R-Squared" , "95% CI Lower Bound" , "95% CI Upper Bound")
      rownames(HV_SNP) <- thresholds
      
      Base_r_squared = summary( lm( paste( hv , " ~ ",paste( confounder_columns , collapse = " + ") ) , data = table , na.action = na.exclude))$r.squared
      print(Base_r_squared)
      i=1
      for( col_name in thresholds ){
        #table[,col_name] <- (table[,col_name] - min(table[,col_name],na.rm = TRUE)) / (max(table[,col_name],na.rm = TRUE) - min(table[,col_name],na.rm = TRUE))
        #table[,col_name] <- (table[,hv] - min(table[,hv],na.rm = TRUE)) / (max(table[,hv],na.rm = TRUE) - min(table[,hv],na.rm = TRUE))
        
        x <- lm( paste( hv , " ~ ",col_name," + " , paste( confounder_columns , collapse = " + ") ) , data = table , na.action = na.exclude)
        
        HV_SNP[i,1] <- summary(x)$coefficients[2,1]
        HV_SNP[i,2] <- max(table[,col_name] , na.rm=TRUE) - min(table[,col_name] , na.rm=TRUE)
        HV_SNP[i,3] <- summary(x)$coefficients[2,4]
        HV_SNP[i,4] <- summary(x)$r.squared
        HV_SNP[i,5] <- as.double(CI.Rsq(rsq = summary(x)$r.squared , n = nrow(table) , k = 13)[3])
        HV_SNP[i,6] <- as.double(CI.Rsq(rsq = summary(x)$r.squared , n = nrow(table) , k = 13)[4])
        i<-i+1
      }
      
      #print(paste(GENDER , HEMISPHERE))
      #print(HV_SNP)
      
      #      if (GENDER=="Both Genders")
      #        if (HEMISPHERE=="Both Hemispheres"){
      
      bplot <- barplot(main = paste(GENDER , HEMISPHERE) , xpd = FALSE ,ylim = c(min(HV_SNP[,4])-0.001,max(HV_SNP[,4])+0.0025),
                       HV_SNP[,4] , xlab = "P-Value Threshold" , ylab = "R^2" , 
                       names= substring(thresholds, 8) , las=2)
      
      text( x=c(1:14)*1.2,  y=HV_SNP[,4]+0.0001, formatC(HV_SNP[,3], format = "e", digits = 2) , srt=90 , adj=c(-0.1,-0.5))
      #        }
      max_col_name <-  names(HV_SNP[,4])[match(max(HV_SNP[,4]),HV_SNP[,4])]
      table <- table[  order(table[,max_col_name]) , ]
      hv <- table[ , hv]
      
      #      if (GENDER=="Both Genders")
      #        if (HEMISPHERE=="Both Hemispheres"){
      means <- 1:100
      percs <- 1:100
      for ( i in 1:100 ){
        # THIS METHOD CAlCULATESS MEANS AS PERCENTILE OF SAMPLES
        start_ind = nrow(table) * ((i-1)/100)
        end_ind = nrow(table) * (i/100)
        window <- hv[start_ind:end_ind]
        means[i] <- mean(window , na.rm = TRUE)
        
      }
      
      
      #          var_vals <- 1:45
      #          var_diff_old <- 0
      #          for( i in 1:45 ){
      #            means_without <- tail( head(means , length(means)-i) , length(means)-(2*i))
      #            var_diff <- ( (var(means) - var(means_without)) / var(means) ) * 100
      #            var_vals[i] <- var_diff-var_diff_old
      #            print(paste("at",i,"incremental var diff is",var_vals[i],"%"))
      #            var_diff_old <- var_diff
      #          }
      
      plot(means[1:96] , xlab = "percentile" , ylab = "HV" , main = paste("Percentiles for Threshold: " , "0.75"))# max_col_name) )
      lines( predict(lm( head(means,-3) ~ poly(head(percs,-3),3) )) , col="grey" )
      
      
      #          means <- means[!is.na(means)]
      #          percent_of_range <- seq(0,50)*0
      #          for( i in seq(1,45) ){
      #          top_range <- range(head(means,i))[2] - range(head(means,i))[1]
      #          bot_range <- range(tail(means,i))[2] - range(tail(means,i))[1]
      #          full_range <- range(means)[2] - range(means)[1]
      #          percent_of_range[i] <- (top_range + bot_range) / full_range
      #          }
      #          par(mfrow=c(1,1))
      #          plot(percent_of_range , main = "PERCENT OF VARINACE EXPLAINED")
      #        }
    }
  }
  return(HV_SNP)
}



GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


shift_model = function(table , eval_col , shift_zscore_col , base_zscore_col , extreme_exclude_th = 11 , plt = TRUE , title = "" , xlab = "Percentile Shift" , ylab = ""  ,  yaxt = "s" , ylimit = NULL , legend = FALSE){
  
  table$shift = ( pnorm(unlist(table[ , shift_zscore_col ])) -
                    pnorm(unlist(table[ , base_zscore_col ])) )* 100
  
  table = table[ complete.cases(table[,c(eval_col, base_zscore_col , "age" , "sex" , "ICV_nc_gpc.bl" , "shift")]) , ]
  
  model_shift = lm( as.formula(paste0( eval_col," ~ ", base_zscore_col , " + age + sex + ICV_nc_gpc.bl + shift" )) , data = table)

  intercept = summary(model_shift)$coefficients["(Intercept)" , "Estimate"]
  slope = summary(model_shift)$coefficients["shift" , "Estimate"]
  pval = summary(model_shift)$coefficients["shift" , "Pr(>|t|)"]
  
  model_shift_no_extremes = lm( as.formula(paste0( eval_col," ~ ", base_zscore_col , " + age + sex + ICV_nc_gpc.bl + shift" )) , data = table[ abs(table$shift) < extreme_exclude_th , ])
  
  intercept_no_extreme = summary(model_shift_no_extremes)$coefficients["(Intercept)" , "Estimate"]
  slope_no_extreme = summary(model_shift_no_extremes)$coefficients["shift" , "Estimate"]
  pval_no_extreme = summary(model_shift_no_extremes)$coefficients["shift" , "Pr(>|t|)"]
  
  if(plt){
    sig_th = 0.05
    
    higher_better = (mean(resid(model_shift)) - median(resid(model_shift))) < 0 
    
    table$residuals = intercept + resid(model_shift) + (slope * table$shift)
    max_shift =  max(table$shift, na.rm = TRUE)
    min_shift =  min(table$shift, na.rm = TRUE)
    
    table$shift_group = factor( ifelse(table$shift < -extreme_exclude_th , paste0("Shift < -",extreme_exclude_th) ,
                                       ifelse(table$shift < 0 , paste0("-",extreme_exclude_th," < Shift < 0") ,
                                              ifelse(table$shift < extreme_exclude_th , paste0("0 < Shift < ",extreme_exclude_th) ,
                                                     ifelse( is.na(table$shift),NA, paste0(extreme_exclude_th," < Shift") )) ) ),
                                levels = c( paste0("Shift < -",extreme_exclude_th) ,  paste0("-",extreme_exclude_th," < Shift < 0") ,
                                            paste0("0 < Shift < ",extreme_exclude_th) , paste0(extreme_exclude_th," < Shift") ) )
    
    get_y = function( x ) { intercept + (slope * x) }
    get_y_no_extreme = function( x ) { intercept_no_extreme + (slope_no_extreme * x) }
    
    lm_lines = data.frame(
      x = c(0 , min_shift , 0 , min_shift) , 
      xend = c( max_shift , 0 , max_shift , 0) , 
      y = c( get_y(0) , get_y(min_shift) , get_y_no_extreme(0) , get_y_no_extreme(min_shift)) , 
      yend = c( get_y(max_shift) , get_y(0) , get_y_no_extreme(max_shift) , get_y_no_extreme(0)) ,
      color = factor( c("postive shift / all subjects" , "negative shift / all subjects" , "postive shift / grey subjects" , "negative shift / grey subjects") ,
                     levels = c("postive shift / all subjects" , "negative shift / all subjects" , "postive shift / grey subjects" , "negative shift / grey subjects") ) 
    )
    
    p = ggplot() +
        geom_point( data = table , aes(x = shift, y = residuals , color = shift_group , shape = shift_group , size = shift_group) ) + 
        geom_segment(data = lm_lines , aes(x = x, xend = xend, y = y , yend = yend, color = color ,  linetype = color) , linewidth = 0.6) +
      
        scale_color_manual( values = c("grey50" , "grey50" , "lightslateblue" , "tomato"  , "tomato" , "lightslateblue" , "lightslateblue" , "tomato") ) + 
        scale_shape_manual( values = c(17 , 2 , 1 , 16 , NA , NA , NA , NA) ) + 
        scale_size_manual( values = c(1.75 , 1.5 , 1.5 , 1.75 ,NA , NA , NA , NA) ) + 
        scale_linetype_manual(values = c("solid"  , "solid" , "dashed" , "dashed" , NA , NA , NA , NA ) ) +

        xlab(xlab) + ylab(ylab) +
        ggtitle( title ) +
        theme(plot.title = element_text(size = 22) )
    if(legend)
      p = p + guides(
        shape = guide_legend(title = "Percentile Shift Group" , override.aes = list(color = c("tomato" , "grey50"  , "grey50" , "lightslateblue" ))),
        linetype = guide_legend(title = "Regression Lines", override.aes = list(color = c( "lightslateblue" , "tomato" , "lightslateblue" , "tomato" ) )),
        color = FALSE , 
        size = FALSE
      ) + theme( legend.position = c(.15, (.8 - higher_better*.6) ),
                 legend.key.height = unit(.25, 'cm'),
                 legend.key.width = unit(1, 'cm'),
                 legend.key = element_blank() , 
                 legend.title = element_text(size=5),
                 legend.text = element_text(size=4),
                 legend.background=element_blank(), 
                 legend.margin = margin(-0.25,0,0,0, unit="cm")
                 )
    if( !is.null(ylimit) )
        p = p + ylim( ylimit )

    
    return(p)
  }
  
  return( list(slope , pval , slope_no_extreme , pval_no_extreme) )
  
  
}




plot_longitudinal = function( table, eval_col , shift_zscore_col , base_zscore_col , th = 11 , title="" , xlab = "visit_year" , ylab = paste("residual",eval_col) , ylimit = NULL ){
  
  table$shift = ( pnorm(unlist(table[ , shift_zscore_col ])) -
                    pnorm(unlist(table[ , base_zscore_col ])) )* 100
  
  table = table[ complete.cases(table[,c(eval_col, base_zscore_col , shift_zscore_col , "age" , "sex" , "ICV_nc_gpc" , "ICV_nc_gpc.bl" , "shift")]) , ]
  # add visit_year here?
  model_shift = lm( as.formula(paste0( eval_col," ~ ", base_zscore_col , " + age + sex + ICV_nc_gpc + shift " )) , data = table)
  
  intercept = summary(model_shift)$coefficients["(Intercept)" , "Estimate"]
  slope = summary(model_shift)$coefficients["shift" , "Estimate"]
  pval = summary(model_shift)$coefficients["shift" , "Pr(>|t|)"]
  
  table$eval_resid = intercept + resid(model_shift) + (slope * table$shift) 
  
  
  table$ICV_new = (table$ICV_nc_gpc.bl - mean(table$ICV_nc_gpc.bl , na.rm = TRUE)) / sd(table$ICV_nc_gpc.bl , na.rm = TRUE)
  
  #table$shift = -(table$shift)
  
  model = lmer( as.formula(paste0(eval_col," ~ age + sex + ICV_new + ( ",base_zscore_col," + shift)*visit_year + (1 + visit_year | id )")),
                data = table ,
                control = lmerControl(optimizer ="nmkbw"))
  
  model_no_extreme = lmer( as.formula(paste0( eval_col," ~ age + sex + ICV_new + ( ",base_zscore_col," + shift)*visit_year + (1 + visit_year | id )" ))  ,
                           data = table[ abs(table$shift) < th , ] ,
                           control = lmerControl(optimizer ="nmkbw"))
  
  # Generate a new data frame for prediction
  time_range <- seq(min(table$visit_year), max(table$visit_year), length.out = 100)
  genetic_component_range <-c( 20 , 10 , -10 , -20)
  
  new_data <- expand.grid(visit_year = time_range, shift = genetic_component_range)
  new_data$sex <- TRUE
  new_data$id <- table[ table[ table$sex  , eval_col] == median(table[ table$sex , eval_col] , na.rm = TRUE) , "id"][1] 
  new_data$age <- mean(table[table$sex  , "age"] , na.rm = TRUE)
  new_data$ICV_new <- mean(table[table$sex , "ICV_new"] , na.rm = TRUE)
  new_data[ , base_zscore_col] <- mean(table[ table$sex  , base_zscore_col] , na.rm = TRUE)
  
  # Predict
  new_data$predicted_values <- predict(model, newdata = new_data , re.form = NA)
  new_data$predicted_values_no_extreme <- predict(model_no_extreme, newdata = new_data , re.form = NA)
  
  
  m_10 = table[ table$shift <= -th , ]
  p_10 = table[ table$shift >= th , ]
  oth  = table[ abs(table$shift) < th , ]
  oth_p = oth[ oth$shift > 0 , ]
  oth_n = oth[ oth$shift < 0 , ]
  all_p = table[ table$shift > 0 , ]
  all_n = table[ table$shift < 0 , ]
  
  p = ggplot( data = rbind(m_10,p_10) , aes(x= visit_year , y = eval_resid)) +
    geom_point(data = oth_p , aes(x= visit_year , y = eval_resid ) , shape=2 , col =  "grey50" ) +
    geom_point(data = oth_n , aes(x= visit_year , y = eval_resid ) , shape=1 , col =  "grey50" ) +
    geom_point( data = m_10 , aes(x= visit_year , y = eval_resid) , shape=17 , size = 1.75 , col = "tomato" ) +
    geom_point(data = p_10 , aes(x= visit_year , y = eval_resid ) , shape=16 , size = 1.75 , col = "lightslateblue") +
    
    geom_line(data = new_data[ new_data$shift == 10, ], aes(x = visit_year, y = predicted_values) , color = "blue") +
    geom_line(data = new_data[ new_data$shift == 10, ], aes(x = visit_year, y = predicted_values_no_extreme) , color = "blue", linetype = "dashed") +
    geom_line(data = new_data[ new_data$shift == -10, ], aes(x = visit_year, y = predicted_values_no_extreme) , color = "red", linetype = "dashed") +
    geom_line(data = new_data[ new_data$shift == -10, ], aes(x = visit_year, y = predicted_values) , color = "red") +
    
    # geom_smooth( data = all_p , aes(x= visit_year , y = eval_resid) , method = "lm" , formula = y~x ,  se = FALSE , col = "blue" , linetype = "dashed") +
    # geom_smooth( data = all_n , aes(x= visit_year , y = eval_resid), method = "lm" ,  formula = y~x , se = FALSE , col = "red" , linetype = "dashed") +
    # geom_smooth( data = p_10 , aes(x= visit_year , y = eval_resid), method = "lm" ,  formula = y~x , se = FALSE , col = "blue") +
    # geom_smooth( data = m_10 , aes(x= visit_year , y = eval_resid) , method = "lm" , formula = y~x ,  se = FALSE , col = "red") +
    
    ylab(ylab) + xlab(xlab) + ggtitle(title) +
    theme(plot.title = element_text(size = 22))
  
  if( !is.null(ylimit) )
    p = p + ylim( ylimit )
  
  
  return(p)
}


my_boxplot = function( table , cols , xlab = "models" , ylab , outline = FALSE , xaxt = NULL , ylim = NULL ){
  if( is.null(xaxt))
    xaxt = cols
  array = unlist(table[ , cols])
  counts = xaxt[ unlist(c(t(matrix(rep(1:length(cols) , nrow(table)) , ncol = nrow(table))))) ]
  boxplot( array ~ counts , outline = outline, ylab = ylab , xlab = xlab , ylim = ylim)
  
}

make_boxplot = function( t , col , name , labs = c("A" , "A+S" ,"A+S+I" , "A+S+I+P" , "A+S+I+L") , title = "" , yl = NULL){
  
  p = ggplot( t , aes(x = factor(as.integer(factor(exp_grp))) , y = get(col) , group = factor(as.integer(factor(exp_grp))))) + 
    geom_bar( colour="grey" , stat = "summary" , fun = base::mean) +
    geom_boxplot() + 
    ylab(name) +
    scale_x_discrete( name = "Experiment Group",  labels= labs ) + 
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  
  if( is.null(yl) )
    p = p + coord_cartesian(ylim = c(min(t[,col]) , max(t[,col])) )
  else
    p = p + ylim(yl)
  
  
  return(p)
}


make_boxplot_2 = function( t , cols , grp_col = "exp_grp" , ylab , labs = c("A" , "A+S" ,"A+S+I" , "A+S+I+P" , "A+S+I+L") , xlab = c("training","testing","out of sample"), title = "" , yl = NULL , legend_name = "Zscore Group" , rm_legend = FALSE){
  
  t[ , grp_col] = factor(as.integer(factor( t[ , grp_col] )))

  t_big = data.frame(matrix(ncol = 3, nrow = 0))
  names(t_big) = c("exp_grp" , "test_res" , "test_grp")
  for( i in 1:length(cols) ){
    t_temp =  cbind( t[ , c( grp_col , cols[i]) ] , i) 
    names(t_temp) = c("exp_grp" , "test_res" , "test_grp")
    t_big  = rbind(t_big , t_temp )
  }
  
  t_big$test_grp = as.factor(t_big$test_grp)
    
  p = ggplot( t_big , aes(x = test_grp  , y = test_res , fill = exp_grp ) ) + 
    geom_boxplot() + 
    ylab(ylab) +
    scale_x_discrete( name = "Experiment Group",  labels= xlab ) + 
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_discrete(name = legend_name, labels =  labs )
  
  if( is.null(yl) )
    p = p + coord_cartesian(ylim = c(min(t_big$test_res ) , max(t_big$test_res )) )
  else
    p = p + ylim(yl)
  
  if(rm_legend)
    p = p +  scale_fill_discrete(guide="none")
  
  return(p)
}
