############################################################
## Roadside proximity extraction for species occurrence    #
## data.                                                   #
##                                                         #
## Author: Vaughn Shirey                                   #
############################################################

library(rgeos)
library(raster)

set.seed(11082017)

## read in occurrence data and raster files ##

o <- read.csv(file.choose(), header=TRUE)
r <- raster(file.choose())
road <- raster(file.choose())

m <- c(-Inf,0,0, 1,11,1, 12,15,2, 16,19,3, 20,21,4, 22,23,5, 24,26,6, 27,29,7, 30,31,8, 32,36,9, 37,39,10, 40,43,11, 44,45,12, 46,Inf,13)
rclmat <- matrix(m, ncol=3, byrow=TRUE)

r2 <- reclassify(r, rclmat, right=NA)

## truncate occurrence data to necessary variables and convert to spatial points ##

o.simp <- subset(o, select=c(Taxon.Scie, Gatherin26, Gatherin27, Gatherin28, Gatherin29))
o.points <- SpatialPoints(cbind((o.simp$Gatherin29+o.simp$Gatherin28)/2, (o.simp$Gatherin27+o.simp$Gatherin26)/2))

########################################################################################################
## process the first 50 occurrences by sampling the same habitat within a 1000 meter radius 100 times ##
########################################################################################################

prox <- vector()
prox.sim <- vector()
sim <- list()
j = 0

for(i in 1:50){
  
  p <- o.points[i,]
  prox[i] <- extract(road, p)
  
  while(j < 101){
    e <- extract(r2, p)
    
    x.f <- runif(1, -1000.0, 1000.0)
    y.f <- runif(1, -1000.0, 1000.0)
    
    df <- as.data.frame(p)
    
    df$coords.x1 <- df$coords.x1+x.f
    df$coords.x2 <- df$coords.x2+y.f
    
    p <- SpatialPoints(df)
    
    e2 <- extract(r2, p)
    
    print(paste("Iteration: ", i, ".", j, " | Original: ", e, " | Simulated: ", e2, sep=""))
    
    if(e == e2){
      
      prox.sim[j] <- extract(road, p)
      j = j+1
    }
  }
  sim[[i]] <- prox.sim
  j = 0
}

## plot the two datasets for comparison ##

par(mfrow=c(1,2))
hist(prox, main="Original Proximity to Roadsides")
hist(unlist(sim), main="Proximity of Simulated Data to Roadsides")

#######################################################################################
## process the first 50 occurrences that do not occur in an urban classified habitat ##
#######################################################################################

prox <- vector()
prox.sim <- vector()
sim <- list()
j = 0

for(i in 1:50){
  
  p <- o.points[i,]
  prox[i] <- extract(road, p)
  e <- extract(r2, p)
  
  if(e > 4){
    while(j < 101){
      
      
      x.f <- runif(1, -1000.0, 1000.0)
      y.f <- runif(1, -1000.0, 1000.0)
      
      df <- as.data.frame(p)
      
      df$coords.x1 <- df$coords.x1+x.f
      df$coords.x2 <- df$coords.x2+y.f
      
      p <- SpatialPoints(df)
      
      e2 <- extract(r2, p)
      
      print(paste("Iteration: ", i, ".", j, " | Original: ", e, " | Simulated: ", e2, sep=""))
      
      if(e == e2){
        prox.sim[j] <- extract(road, p)
        j = j+1
      }
    }
    sim[[i]] <- prox.sim
    j = 0
  }
}

par(mfrow=c(1,2))
hist(prox, main="Original Proximity to Roadsides")
hist(unlist(sim), main="Proximity of Simulated Data to Roadsides (Land Class Restricted)")
