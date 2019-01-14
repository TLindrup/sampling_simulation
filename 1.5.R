# lav eventuelt flere gentagelser af samplings

        ##########       ####        
      ##########       ####### 
      ###            #####  ### 
        #####       #####  ####
           ###     ############
   ###########    #####   ####
  ##########     #####   #####


rm(list=ls())

.libPaths("/Users/thomaslindrup/Desktop/_/FNM/R/pkgs/")

pack<-c("car","sandwich","lmtest","RColorBrewer","mgcv","foreign","xtable","magrittr",
        "AER","stargazer", "MASS", "foreign", "httr", "XML","xml2","dplyr","gmodels",
        "rvest","purrr","rgdal", "classInt", "GISTools", "maps", "dismo", "raster", 
        "shape","ggmap","knitr", "reshape", "rgeos","spdep","sp","maptools","corrgram",
        "tidyr","stringi","spatstat")

lapply(pack, require, character.only=T)

# setting directories
data <-"/Users/thomaslindrup/Desktop/_/FNM/Applied Statistics/Data/Projekt datas\303\246t"

setwd(data)

# looking in directory (data)
dir()

## Loading data

#10 -12 yrs old forest

data <- read.csv(file = "Afd110pTrolleholm1991.csv")

######### PLOTTING FOR OVERVIEW #########

### the forest ...
forest_plot <- ggplot(data, aes(X_dm, Y_dm))
  forest_plot + geom_point(aes(color = Species), size = 1.5, alpha = 0.9) + 
    labs(title = NULL) + labs(x = NULL, y = NULL) + 
      scale_color_brewer(palette="YlGn")

## distribution of damage - using facets
sp.comp <- ggplot(data, aes(X_dm, Y_dm))
  sp.comp + geom_point(aes(color = Species), size = .5, alpha = 1) + 
    facet_grid(Dmg~.) + labs(title = NULL) + labs(x = NULL, y = NULL) +
    scale_color_brewer(palette="Set2")
  
table(data$Species,data$Dmg)
  
# It appears that it would be beneficial to solely focus on spruce - as it includes a lot of observations 
# which are distributed nicely irt. damage sustained
# we also need to transform H, X and Y into meters:

spruce <- transform(subset(data,Species=="Spruce"),X=X_dm/10,Y=Y_dm/10,H=H_dm/10)
spruce <- spruce[,c("Dmg","X","Y","H")]
summary(spruce)
  
## distribution of damage - using facets = ONLY spruce
sp.comp <- ggplot(spruce, aes(X, Y))
sp.comp + geom_point(aes(col = "Orange"),size = 2, alpha = 0.1) + 
  facet_grid(Dmg~.) + labs(title = NULL) + labs(x = NULL, y = NULL)

######### Model creation and validation #########

#### Modelling using generalised additative model - for all of the area!! ####

# Using Dmg as a continual response
library("gam")
m0 <- gam(Dmg~H+lo(X,Y),data=spruce)

# Using Dmg as a binary response
m1 <- gam(I(Dmg>0)~H+lo(X,Y),data=spruce,family=binomial())

# Residual plots
plot(predict(m0),residuals(m0))

qqnorm(residuals(m0))

# This proves that m0 is not valid. As Dmg is not continuous. 
# GAM is not able to handle 'ordinal' data, but we decide that this is enough for us atm.
# Otherwise we would need to get into some very complicated stuff in order to be able to 
# Mathematically describe the spatial relationship
# m1 can not be validated using the 'gof' package and we must 
# assume that it is valid

#### Effect of H on Dmg

# Slope coefficients
coef(m0)

coef(m0)["H"]

summary(m0)

summary(m0)$parametric.anova["H","Pr(>F)"]

# Spatial prediction
# what would the spatial dmg distribution look like if all trees were 2 meters 
# and distributed in a grid
x <- seq(min(spruce$X),max(spruce$X),length.out=50) # del X op i 50 lige store dele
y <- seq(min(spruce$Y),max(spruce$Y),length.out=50)
test.grid <- data.frame(H=2,X=c(outer(x,rep(1,length(y)))),Y=c(outer(rep(1,length(x)),y)))
summary(test.grid)

hat.Dmg <- predict(m0,newdata=test.grid)

## plotting to illustrate differences in predictions and model

library(ggplot2)
ggplot(cbind(test.grid,hat.Dmg),aes(x=X,y=Y)) +
  geom_tile(aes(fill=hat.Dmg)) +
  ggtitle("Spruce: Damage")

# Subsampling example
# taking out 10% of the data ...
smalldata <- spruce[sample.int(nrow(spruce),round(0.1*nrow(spruce))),]

# and running the model on this data
m0.small <- gam(Dmg~H+lo(X,Y),data=smalldata)

# effect of the height
coef(m0.small)

# plot
hat.Dmg.small <- predict(m0.small,newdata=test.grid)
ggplot(cbind(test.grid,hat.Dmg.small),aes(x=X,y=Y)) + 
  geom_tile(aes(fill=hat.Dmg.small)) + 
  ggtitle("Spruce: Damage udfra 10% af data")

# how much has the prediction changed?
ggplot(cbind(test.grid,fejl=hat.Dmg-hat.Dmg.small),aes(x=X,y=Y)) + 
  geom_tile(aes(fill=fejl)) + 
  ggtitle("Spruce: Fejl i pr??diktion af damage")

# calculating root mean squared error => (0 => perfect fit)
sqrt(mean((hat.Dmg-hat.Dmg.small)^2))

############### ############### ############### ####
############### SAMPLING SIMULATIONS ###############
#### ############### ############### ############### 

## setting up 'test.grid'
x <- seq(min(spruce$X),max(spruce$X),length.out=50) # del X op i 50 lige store dele
y <- seq(min(spruce$Y),max(spruce$Y),length.out=50)
test.grid <- data.frame(H=2,X=c(outer(x,rep(1,length(y)))),Y=c(outer(rep(1,length(x)),y)))

# lav pr??diktion med alle data
hat.Dmg <- predict(m0,newdata=test.grid)

             #### R A N D O M ####

   #####################################
  #### RANDOM PLOT SAMPLING LOOP 1.5 ####
   #####################################

estimator <- matrix(0,400,4)
runs <- 400
psize <- 1.5 # i meter
nsize <- 10 # n of samples

for (i in 1:runs) {
  x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
  y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
  
  ii <- sapply(1:nrow(spruce),
               function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
  smalldata_rnd_1.5m <- spruce[ii,]
  
  # Model from samples
  m0_rnd_1.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_1.5m)

  # Saving results
  estimator[i,] <- coef(m0_rnd_1.5m) ## coef(m0_rnd_1.5m)["H"]

}

  # P value for "H"
  summary(m0_rnd_1.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_1.5m <- predict(m0_rnd_1.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_1.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_1.5 <- estimator
  summary(rnd_1.5[,2])

   #####################################
  #### RANDOM PLOT SAMPLING LOOP 2.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 2.5 # i meter
  nsize <- 10 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_2.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_2.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_2.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_2.5m) ## coef(m0_rnd_2.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_2.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_2.5m <- predict(m0_rnd_2.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_2.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_2.5 <- estimator
  summary(rnd_2.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 3.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 3.5 # i meter
  nsize <- 10 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_3.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_3.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_3.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_3.5m) ## coef(m0_rnd_3.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_3.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_3.5m <- predict(m0_rnd_3.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_3.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_3.5 <- estimator
  summary(rnd_3.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 4.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 4.5 # i meter
  nsize <- 10 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    # Filter command (filter %>% might be shorter) - 'spruce %>% filter(min(abs(X_dm-x)<1),min(abs(Y_dm-y)<1)))'
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_4.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_4.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_4.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_4.5m) ## coef(m0_rnd_4.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_4.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_4.5m <- predict(m0_rnd_4.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_4.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_4.5 <- estimator
  summary(rnd_4.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 5.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 5.5 # i meter
  nsize <- 10 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_5.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_5.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_5.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_5.5m) ## coef(m0_rnd_5.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_5.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_5.5m <- predict(m0_rnd_5.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_5.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_5.5 <- estimator
  summary(rnd_5.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 7.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 7.5 # i meter
  nsize <- 10 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_7.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_7.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_7.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_7.5m) ## coef(m0_rnd_7.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_7.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_7.5m <- predict(m0_rnd_7.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_7.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_7.5 <- estimator
  summary(rnd_7.5[,2])
  
   #################################### 
  #### RANDOM PLOT SAMPLING LOOP 10 ####
   #################################### 
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 10 # i meter
  nsize <- 10 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_10m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_10m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_10m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_10m) ## coef(m0_rnd_10m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_10m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_10m <- predict(m0_rnd_10m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_10m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_10 <- estimator
  summary(rnd_10[,2])

   #################################### 
  #### RANDOM PLOT SAMPLING LOOP 15 ####
   #################################### 
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 15 # i meter
  nsize <- 10 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_15m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_15m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_15m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_15m) ## coef(m0_rnd_15m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_15m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_15m <- predict(m0_rnd_15m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_15m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_15 <- estimator
  summary(rnd_15[,2])
  
   ###################################
  #### OVERVIEW OF RANDOM SAMPLING #### 
   ###################################  
  
  # 400 runs and 10 plots
  
  # Summary of the "H" coefficient across all repeated samplings
  dat_rnd_10ps <- matrix(0,400,8)
  dat_rnd_10ps[,1] <- rnd_1.5[,2]
  dat_rnd_10ps[,2] <- rnd_2.5[,2]
  dat_rnd_10ps[,3] <- rnd_3.5[,2]
  dat_rnd_10ps[,4] <- rnd_4.5[,2]
  dat_rnd_10ps[,5] <- rnd_5.5[,2]
  dat_rnd_10ps[,6] <- rnd_7.5[,2]
  dat_rnd_10ps[,7] <- rnd_10[,2]
  dat_rnd_10ps[,8] <- rnd_15[,2]
  colnames(dat_rnd_10ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  
  # RMSE - the spread of the prediction error
  RMSE_rnd_10ps <- matrix(0,8,1)
  RMSE_rnd_10ps[1,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_1.5m)^2))
  RMSE_rnd_10ps[2,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_2.5m)^2))
  RMSE_rnd_10ps[3,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_3.5m)^2))
  RMSE_rnd_10ps[4,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_4.5m)^2))
  RMSE_rnd_10ps[5,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_5.5m)^2))
  RMSE_rnd_10ps[6,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_7.5m)^2))
  RMSE_rnd_10ps[7,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_10m)^2))
  RMSE_rnd_10ps[8,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_15m)^2))
  rownames(RMSE_rnd_10ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  colnames(RMSE_rnd_10ps) <- "RootMeanSquaredError"
  plot(RMSE_rnd_10ps)
  
  # P value for "H"
  pforH_rnd_10ps <- matrix(0,8,1)
  pforH_rnd_10ps[1,] <- summary(m0_rnd_1.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_10ps[2,] <- summary(m0_rnd_2.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_10ps[3,] <- summary(m0_rnd_3.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_10ps[4,] <- summary(m0_rnd_4.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_10ps[5,] <- summary(m0_rnd_5.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_10ps[6,] <- summary(m0_rnd_7.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_10ps[7,] <- summary(m0_rnd_10m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_10ps[8,] <- summary(m0_rnd_15m)$parametric.anova["H","Pr(>F)"]
  rownames(pforH_rnd_10ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  colnames(pforH_rnd_10ps) <- "The P value for the Height coefficient"

  # PLOTTING RESULTS
  boxplot_rnd_10ps <- boxplot(dat_rnd_10ps, use.cols=TRUE, title("Spread of the H coeff."))
  plot(RMSE_rnd_10ps) # RMSE of prediction
  plot(pforH_rnd_10ps) # P value for H coeff.
 
    ##### 7 plots random sampling #####
  
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 1.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 1.5 # i meter
  nsize <- 7 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_1.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_1.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_1.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_1.5m) ## coef(m0_rnd_1.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_1.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_1.5m <- predict(m0_rnd_1.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_1.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_1.5 <- estimator
  summary(rnd_1.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 2.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 2.5 # i meter
  nsize <- 7 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_2.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_2.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_2.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_2.5m) ## coef(m0_rnd_2.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_2.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_2.5m <- predict(m0_rnd_2.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_2.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_2.5 <- estimator
  summary(rnd_2.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 3.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 3.5 # i meter
  nsize <- 7 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_3.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_3.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_3.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_3.5m) ## coef(m0_rnd_3.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_3.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_3.5m <- predict(m0_rnd_3.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_3.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_3.5 <- estimator
  summary(rnd_3.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 4.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 4.5 # i meter
  nsize <- 7 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_4.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_4.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_4.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_4.5m) ## coef(m0_rnd_4.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_4.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_4.5m <- predict(m0_rnd_4.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_4.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_4.5 <- estimator
  summary(rnd_4.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 5.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 5.5 # i meter
  nsize <- 7 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_5.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_5.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_5.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_5.5m) ## coef(m0_rnd_5.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_5.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_5.5m <- predict(m0_rnd_5.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_5.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_5.5 <- estimator
  summary(rnd_5.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 7.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 7.5 # i meter
  nsize <- 7 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_7.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_7.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_7.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_7.5m) ## coef(m0_rnd_7.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_7.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_7.5m <- predict(m0_rnd_7.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_7.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_7.5 <- estimator
  summary(rnd_7.5[,2])
  
   #################################### 
  #### RANDOM PLOT SAMPLING LOOP 10 ####
   #################################### 
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 10 # i meter
  nsize <- 7 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_10m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_10m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_10m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_10m) ## coef(m0_rnd_10m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_10m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_10m <- predict(m0_rnd_10m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_10m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_10 <- estimator
  summary(rnd_10[,2])
  
   #################################### 
  #### RANDOM PLOT SAMPLING LOOP 15 ####
   #################################### 
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 15 # i meter
  nsize <- 7 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_15m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_15m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_15m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_15m) ## coef(m0_rnd_15m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_15m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_15m <- predict(m0_rnd_15m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_15m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_15 <- estimator
  summary(rnd_15[,2])
  
   ###################################
  #### OVERVIEW OF RANDOM SAMPLING #### 
   ###################################  
  
  # 400 runs and 7 plots
  
  # Summary of the "H" coefficient across all repeated samplings
  dat_rnd_7ps <- matrix(0,400,8)
  dat_rnd_7ps[,1] <- rnd_1.5[,2]
  dat_rnd_7ps[,2] <- rnd_2.5[,2]
  dat_rnd_7ps[,3] <- rnd_3.5[,2]
  dat_rnd_7ps[,4] <- rnd_4.5[,2]
  dat_rnd_7ps[,5] <- rnd_5.5[,2]
  dat_rnd_7ps[,6] <- rnd_7.5[,2]
  dat_rnd_7ps[,7] <- rnd_10[,2]
  dat_rnd_7ps[,8] <- rnd_15[,2]
  colnames(dat_rnd_7ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  
  # RMSE - the spread of the prediction error
  RMSE_rnd_7ps <- matrix(0,8,1)
  RMSE_rnd_7ps[1,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_1.5m)^2))
  RMSE_rnd_7ps[2,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_2.5m)^2))
  RMSE_rnd_7ps[3,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_3.5m)^2))
  RMSE_rnd_7ps[4,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_4.5m)^2))
  RMSE_rnd_7ps[5,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_5.5m)^2))
  RMSE_rnd_7ps[6,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_7.5m)^2))
  RMSE_rnd_7ps[7,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_10m)^2))
  RMSE_rnd_7ps[8,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_15m)^2))
  rownames(RMSE_rnd_7ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  colnames(RMSE_rnd_7ps) <- "RootMeanSquaredError"
  plot(RMSE_rnd_7ps)
  
  # P value for "H"
  pforH_rnd_7ps <- matrix(0,8,1)
  pforH_rnd_7ps[1,] <- summary(m0_rnd_1.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_7ps[2,] <- summary(m0_rnd_2.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_7ps[3,] <- summary(m0_rnd_3.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_7ps[4,] <- summary(m0_rnd_4.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_7ps[5,] <- summary(m0_rnd_5.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_7ps[6,] <- summary(m0_rnd_7.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_7ps[7,] <- summary(m0_rnd_10m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_7ps[8,] <- summary(m0_rnd_15m)$parametric.anova["H","Pr(>F)"]
  rownames(pforH_rnd_7ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  colnames(pforH_rnd_7ps) <- "The P value for the Height coefficient"
  
  # PLOTTING RESULTS
  boxplot_rnd_7ps <- boxplot(dat_rnd_7ps, use.cols=TRUE, title("Spread of the H coeff."))
  plot(RMSE_rnd_7ps) # RMSE of prediction
  plot(pforH_rnd_7ps) # P value for H coeff.
  
  ##### 5 plots random sampling #####
  
  
  #####################################
  #### RANDOM PLOT SAMPLING LOOP 1.5 ####
  #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 1.5 # i meter
  nsize <- 5 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_1.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_1.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_1.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_1.5m) ## coef(m0_rnd_1.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_1.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_1.5m <- predict(m0_rnd_1.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_1.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_1.5 <- estimator
  summary(rnd_1.5[,2])
  
  #####################################
  #### RANDOM PLOT SAMPLING LOOP 2.5 ####
  #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 2.5 # i meter
  nsize <- 5 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_2.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_2.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_2.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_2.5m) ## coef(m0_rnd_2.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_2.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_2.5m <- predict(m0_rnd_2.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_2.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_2.5 <- estimator
  summary(rnd_2.5[,2])
  
  #####################################
  #### RANDOM PLOT SAMPLING LOOP 3.5 ####
  #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 3.5 # i meter
  nsize <- 5 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_3.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_3.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_3.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_3.5m) ## coef(m0_rnd_3.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_3.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_3.5m <- predict(m0_rnd_3.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_3.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_3.5 <- estimator
  summary(rnd_3.5[,2])
  
  #####################################
  #### RANDOM PLOT SAMPLING LOOP 4.5 ####
  #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 4.5 # i meter
  nsize <- 5 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_4.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_4.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_4.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_4.5m) ## coef(m0_rnd_4.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_4.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_4.5m <- predict(m0_rnd_4.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_4.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_4.5 <- estimator
  summary(rnd_4.5[,2])
  
  #####################################
  #### RANDOM PLOT SAMPLING LOOP 5.5 ####
  #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 5.5 # i meter
  nsize <- 5 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_5.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_5.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_5.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_5.5m) ## coef(m0_rnd_5.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_5.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_5.5m <- predict(m0_rnd_5.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_5.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_5.5 <- estimator
  summary(rnd_5.5[,2])
  
  #####################################
  #### RANDOM PLOT SAMPLING LOOP 7.5 ####
  #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 7.5 # i meter
  nsize <- 5 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_7.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_7.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_7.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_7.5m) ## coef(m0_rnd_7.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_7.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_7.5m <- predict(m0_rnd_7.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_7.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_7.5 <- estimator
  summary(rnd_7.5[,2])
  
  #################################### 
  #### RANDOM PLOT SAMPLING LOOP 10 ####
  #################################### 
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 10 # i meter
  nsize <- 5 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_10m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_10m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_10m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_10m) ## coef(m0_rnd_10m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_10m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_10m <- predict(m0_rnd_10m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_10m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_10 <- estimator
  summary(rnd_10[,2])
  
  #################################### 
  #### RANDOM PLOT SAMPLING LOOP 15 ####
  #################################### 
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 15 # i meter
  nsize <- 5 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_15m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_15m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_15m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_15m) ## coef(m0_rnd_15m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_15m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_15m <- predict(m0_rnd_15m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_15m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_15 <- estimator
  summary(rnd_15[,2])
  
   ###################################
  #### OVERVIEW OF RANDOM SAMPLING #### 
   ###################################  
  
  # 400 runs and 5 plots
  
  # Summary of the "H" coefficient across all repeated samplings
  dat_rnd_5ps <- matrix(0,400,8)
  dat_rnd_5ps[,1] <- rnd_1.5[,2]
  dat_rnd_5ps[,2] <- rnd_2.5[,2]
  dat_rnd_5ps[,3] <- rnd_3.5[,2]
  dat_rnd_5ps[,4] <- rnd_4.5[,2]
  dat_rnd_5ps[,5] <- rnd_5.5[,2]
  dat_rnd_5ps[,6] <- rnd_7.5[,2]
  dat_rnd_5ps[,7] <- rnd_10[,2]
  dat_rnd_5ps[,8] <- rnd_15[,2]
  colnames(dat_rnd_5ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  
  # RMSE - the spread of the prediction error
  RMSE_rnd_5ps <- matrix(0,8,1)
  RMSE_rnd_5ps[1,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_1.5m)^2))
  RMSE_rnd_5ps[2,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_2.5m)^2))
  RMSE_rnd_5ps[3,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_3.5m)^2))
  RMSE_rnd_5ps[4,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_4.5m)^2))
  RMSE_rnd_5ps[5,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_5.5m)^2))
  RMSE_rnd_5ps[6,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_7.5m)^2))
  RMSE_rnd_5ps[7,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_10m)^2))
  RMSE_rnd_5ps[8,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_15m)^2))
  rownames(RMSE_rnd_5ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  colnames(RMSE_rnd_5ps) <- "RootMeanSquaredError"
  plot(RMSE_rnd_5ps)
  
  # P value for "H"
  pforH_rnd_5ps <- matrix(0,8,1)
  pforH_rnd_5ps[1,] <- summary(m0_rnd_1.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_5ps[2,] <- summary(m0_rnd_2.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_5ps[3,] <- summary(m0_rnd_3.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_5ps[4,] <- summary(m0_rnd_4.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_5ps[5,] <- summary(m0_rnd_5.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_5ps[6,] <- summary(m0_rnd_7.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_5ps[7,] <- summary(m0_rnd_10m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_5ps[8,] <- summary(m0_rnd_15m)$parametric.anova["H","Pr(>F)"]
  rownames(pforH_rnd_5ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  colnames(pforH_rnd_5ps) <- "The P value for the Height coefficient"
  
  # PLOTTING RESULTS
  boxplot_rnd_5ps <- boxplot(dat_rnd_5ps, use.cols=TRUE, title("Spread of the H coeff."))
  plot(RMSE_rnd_5ps) # RMSE of prediction
  plot(pforH_rnd_5ps) # P value for H coeff.
 
  ##### 3 plots random sampling #####
  
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 1.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 1.5 # i meter
  nsize <- 3 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_1.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_1.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_1.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_1.5m) ## coef(m0_rnd_1.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_1.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_1.5m <- predict(m0_rnd_1.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_1.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_1.5 <- estimator
  summary(rnd_1.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 2.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 2.5 # i meter
  nsize <- 3 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_2.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_2.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_2.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_2.5m) ## coef(m0_rnd_2.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_2.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_2.5m <- predict(m0_rnd_2.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_2.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_2.5 <- estimator
  summary(rnd_2.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 3.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 3.5 # i meter
  nsize <- 3 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_3.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_3.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_3.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_3.5m) ## coef(m0_rnd_3.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_3.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_3.5m <- predict(m0_rnd_3.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_3.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_3.5 <- estimator
  summary(rnd_3.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 4.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 4.5 # i meter
  nsize <- 3 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_4.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_4.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_4.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_4.5m) ## coef(m0_rnd_4.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_4.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_4.5m <- predict(m0_rnd_4.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_4.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_4.5 <- estimator
  summary(rnd_4.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 5.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 5.5 # i meter
  nsize <- 3 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_5.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_5.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_5.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_5.5m) ## coef(m0_rnd_5.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_5.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_5.5m <- predict(m0_rnd_5.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_5.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_5.5 <- estimator
  summary(rnd_5.5[,2])
  
   #####################################
  #### RANDOM PLOT SAMPLING LOOP 7.5 ####
   #####################################
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 7.5 # i meter
  nsize <- 3 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_7.5m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_7.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_7.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_7.5m) ## coef(m0_rnd_7.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_7.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_7.5m <- predict(m0_rnd_7.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_7.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_7.5 <- estimator
  summary(rnd_7.5[,2])
  
   #################################### 
  #### RANDOM PLOT SAMPLING LOOP 10 ####
   #################################### 
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 10 # i meter
  nsize <- 3 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_10m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_10m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_10m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_10m) ## coef(m0_rnd_10m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_10m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_10m <- predict(m0_rnd_10m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_10m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_10 <- estimator
  summary(rnd_10[,2])
  
   #################################### 
  #### RANDOM PLOT SAMPLING LOOP 15 ####
   #################################### 
  
  estimator <- matrix(0,400,4)
  runs <- 400
  psize <- 15 # i meter
  nsize <- 3 # n of samples
  
  for (i in 1:runs) {
    x <- runif(nsize,min=min(spruce$X),max=max(spruce$X))
    y <- runif(nsize,min=min(spruce$Y),max=max(spruce$Y))
    
    ii <- sapply(1:nrow(spruce),
                 function(i){(min(abs(spruce$X[i]-x))<psize)&(min(abs(spruce$Y[i]-y))<psize)})
    smalldata_rnd_15m <- spruce[ii,]
    
    # Model from samples
    m0_rnd_15m <- gam(Dmg~H+lo(X,Y),data=smalldata_rnd_15m)
    
    # Saving results
    estimator[i,] <- coef(m0_rnd_15m) ## coef(m0_rnd_15m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_rnd_15m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_rnd_15m <- predict(m0_rnd_15m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_rnd_15m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  rnd_15 <- estimator
  summary(rnd_15[,2])
  
   ###################################
  #### OVERVIEW OF RANDOM SAMPLING #### 
   ###################################  
  
  # 400 runs and 3 plots
  
  # Summary of the "H" coefficient across all repeated samplings
  dat_rnd_3ps <- matrix(0,400,8)
  dat_rnd_3ps[,1] <- rnd_1.5[,2]
  dat_rnd_3ps[,2] <- rnd_2.5[,2]
  dat_rnd_3ps[,3] <- rnd_3.5[,2]
  dat_rnd_3ps[,4] <- rnd_4.5[,2]
  dat_rnd_3ps[,5] <- rnd_5.5[,2]
  dat_rnd_3ps[,6] <- rnd_7.5[,2]
  dat_rnd_3ps[,7] <- rnd_10[,2]
  dat_rnd_3ps[,8] <- rnd_15[,2]
  colnames(dat_rnd_3ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  
  # RMSE - the spread of the prediction error
  RMSE_rnd_3ps <- matrix(0,8,1)
  RMSE_rnd_3ps[1,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_1.5m)^2))
  RMSE_rnd_3ps[2,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_2.5m)^2))
  RMSE_rnd_3ps[3,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_3.5m)^2))
  RMSE_rnd_3ps[4,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_4.5m)^2))
  RMSE_rnd_3ps[5,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_5.5m)^2))
  RMSE_rnd_3ps[6,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_7.5m)^2))
  RMSE_rnd_3ps[7,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_10m)^2))
  RMSE_rnd_3ps[8,] <- sqrt(mean((hat.Dmg-hat.Dmg_rnd_15m)^2))
  rownames(RMSE_rnd_3ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  colnames(RMSE_rnd_3ps) <- "RootMeanSquaredError"
  plot(RMSE_rnd_3ps)
  
  # P value for "H"
  pforH_rnd_3ps <- matrix(0,8,1)
  pforH_rnd_3ps[1,] <- summary(m0_rnd_1.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_3ps[2,] <- summary(m0_rnd_2.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_3ps[3,] <- summary(m0_rnd_3.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_3ps[4,] <- summary(m0_rnd_4.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_3ps[5,] <- summary(m0_rnd_5.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_3ps[6,] <- summary(m0_rnd_7.5m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_3ps[7,] <- summary(m0_rnd_10m)$parametric.anova["H","Pr(>F)"]
  pforH_rnd_3ps[8,] <- summary(m0_rnd_15m)$parametric.anova["H","Pr(>F)"]
  rownames(pforH_rnd_3ps) <- c("1.5", "2.5", "3.5", "4.5", "5.5", "7.5", "10", "15")
  colnames(pforH_rnd_3ps) <- "The P value for the Height coefficient"
  
  # PLOTTING RESULTS
  # EXTREME OUTLIERS ON 1.5m PLOTS
  boxplot_rnd_3ps <- boxplot(dat_rnd_3ps, use.cols=TRUE, title("Spread of the H coeff."))
  plot(RMSE_rnd_3ps) # RMSE of prediction
  plot(pforH_rnd_3ps) # P value for H coeff.
  
  # ggplot2 - not working -> not working !!
  ggplot(data = data.frame(dat_rnd_3ps), aes(x = dat_rnd_3ps, y = "")) + 
    geom_boxplot()
  
  #### plot all plots in order to gain overview:
  # BOXPLOT - without outliers!!!
  boxplot_rnd_10ps <- boxplot(dat_rnd_10ps, use.cols=TRUE, main="H coeff. 10 plots, 400 reps", outline=FALSE, xlab="Sample size")
  boxplot_rnd_7ps <- boxplot(dat_rnd_7ps, use.cols=TRUE, main="H coeff. 7 plots, 400 reps", outline=FALSE, xlab="Sample size")
  boxplot_rnd_5ps <- boxplot(dat_rnd_5ps, use.cols=TRUE, main="H coeff. 5 plots, 400 reps", outline=FALSE, xlab="Sample size")
  boxplot_rnd_3ps <- boxplot(dat_rnd_3ps, use.cols=TRUE, main="H coeff. 3 plots, 400 reps", outline=FALSE, xlab="Sample size")
  
  # RMSE
  plot(RMSE_rnd_10ps)
  plot(RMSE_rnd_7ps)
  plot(RMSE_rnd_5ps)
  # /wo 1.5 and 3.5
  ezy <- RMSE_rnd_5ps[-1,]
  ezy2 <- ezy[-2,]
  plot(tmp.plot)
  plot(RMSE_rnd_3ps)
  
  # P value for H coeff.
  plot(pforH_rnd_10ps)
  plot(pforH_rnd_10ps[-1,])
  plot(pforH_rnd_7ps)
  plot(pforH_rnd_7ps[-1,])
  plot(pforH_rnd_5ps)
  plot(pforH_rnd_5ps[-1,])
  plot(pforH_rnd_3ps)
  
  
            #### T R A N S E C T ####

     ########################################
    #### STRIP TRANSECT SAMPLING 1.5width ####
     ########################################  
  
  # Setting up 'container' for data
  estimator <- matrix(0,6,4)
  
  # How many times does a 1.5m transect fit along the X-axis - "the shortest axis"
  possible_transects_1.5 <- (max(spruce$X)-min(spruce$X))/1.5
  
  # To cover the whole area we need to do ~59 transects of 1.5 metres
  # Lets say we do 6 transects over ten times
  # We want to know the coorinates on the xaxis if divided into ~59 intervals ...
  x <- seq(min(spruce$X),max(spruce$X),length.out=59)
  
  # Width of transect / 2
  width <- 0.75
  
  # Now we want to select a coordinate and skip 6 and select that one, skip 6 etc..
  for (i in 1:6) {
    xx <- x[seq(i, length(x), 6)]      # spiller rollen for a,b,...
    
    # Next take the coordinates and select every observation 
    # within 0.75 metres of them along the Y axis
    ii <- sapply(1:nrow(spruce),function(i){(min(abs(spruce$X[i]-xx))<width)})
    
    
    # small dataset
    smalldata_transect_1.5m <- spruce[ii,]
    
    # refit model on small data and compare results with full data
    m0_transect_1.5m <- gam(Dmg~H+lo(X,Y),data=smalldata_transect_1.5m)
    
    # Saving results
    estimator[i,] <- coef(m0_transect_1.5m) ## coef(m0_transect_1.5m)["H"]
    
  }
  
  # P value for "H"
  summary(m0_transect_1.5m)$parametric.anova["H","Pr(>F)"]
  
  # lav pr??diktion med subsampled data
  hat.Dmg_transect_1.5m <- predict(m0_transect_1.5m, newdata=test.grid)
  
  # RMSE (root mean squared error)
  sqrt(mean((hat.Dmg-hat.Dmg_transect_1.5m)^2))
  
  # Summary of the "H" coefficient across all repeated samplings
  transect_1.5 <- estimator
  
  
  # Boxplot 
  colnames(transect_1.5) <- c("","1.5","","")
  
  boxplot_transect_1.5 <- boxplot(transect_1.5[,2], use.cols=TRUE, main="Spread of the H coeff.", xlab="Transect width")
  
  summary(transect_1.5[,2])

  
  
  
  
  