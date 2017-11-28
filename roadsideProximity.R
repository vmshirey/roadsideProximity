############################################################
## Roadside proximity extraction for species occurrence    #
## data.                                                   #
##                                                         #
## Author: Vaughn Shirey                                   #
############################################################

library(rgeos)
library(raster)

set.seed(11082017)
iterations <- 5

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

rand1 <- sample(1:nrow(o), iterations, replace=FALSE)
rand2 <- sample(1:nrow(o), iterations, replace=FALSE)

########################################################################################################
## process the occurrences by sampling the same habitat within a 1000 meter radius 100 times ##
########################################################################################################

prox <- vector()
prox.sim <- vector()
sim <- list()
j = 0

numGreater <- 0
numCloser <- 0
numEqual <- 0

for(i in 1:iterations){
  
  p <- o.points[rand1[i],]
  prox[i] <- extract(road, p)
  
  e <- extract(r2, p)
  
  while(j < 101){

    x.f <- runif(1, -1000.0, 1000.0)
    y.f <- runif(1, -1000.0, 1000.0)
    
    df <- as.data.frame(p)
    
    df$coords.x1 <- df$coords.x1+x.f
    df$coords.x2 <- df$coords.x2+y.f
    
    p <- SpatialPoints(df)
    
    e2 <- extract(r2, p)
    
    print(paste("Iteration: ", i, ".", j, " | Original: ", e, " | Simulated: ", e2, sep=""))
    
    if(e == e2 && !is.na(e2)){
      
      prox.sim[j] <- extract(road, p)
      j = j+1
    }
    if(is.na(e2)){
      print("BAD SAMPLE: Simulation fell off the Earth! Skipping...")
      j = 999
    }
  }
  sim[[i]] <- prox.sim
  j = 0
}

## plot the two datasets for comparison ##

par(mfrow=c(1,2))
hist(prox, main="Original Proximity to Roadsides", col="gray", border="white", breaks=20, ylim=c(0,20), xlim=c(0,8000))
hist(unlist(sim), main="Proximity of Simulated Data to Roadsides", col="gray", border="white", breaks=20, xlim=c(0,8000))

t <- ks.test(prox, unlist(sim))

par(mfrow=c(1,1))
plot(density(prox), main=paste("Density of Occurrence Record Distances to Roadsides (", o[1,]$Taxon.Scie, ")", sep=""), xlab="Distance from Road (meters)", sub=paste("Two Sample Kolmogorov-Smirnov Test:", t$p.value))
lines(density(unlist(sim)), col="red")
mtext(paste("Sample Mean:", round(mean(prox), digits=2), "  |    Simulation Mean:", round(mean(unlist(sim)), digits=2)), side=4)

## compute summary statistics for this run ##

for(k in 1:length(sim)){numGreater = numGreater + sum(sim[[k]] > prox[k])}
for(k in 1:length(sim)){numCloser = numCloser + sum(sim[[k]] < prox[k])}
for(k in 1:length(sim)){numEqual = numEqual + sum(sim[[k]] == prox[k])}

paste("For this simulation", numGreater, "simulations of", i*iterations, "simulations were farther from roadsides than the original data.")
paste("For this simulation", numCloser, "simulations of", i*iterations, "simulations were closer to roadsides than the original data.")
paste("For this simulation", numEqual, "simulations of", i*iterations, "simulations were equal distance to roadsides than the original data.")

#######################################################################################
## process the occurrences that do not occur in an urban classified habitat ##
#######################################################################################

prox <- vector()
prox.sim <- vector()
sim <- list()
j = 0

numGreater <- 0
numCloser <- 0
numEqual <- 0

for(i in 1:iterations){
  
  p <- o.points[rand2[i],]
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
      
      if(e == e2  && !is.na(e2)){
        prox.sim[j] <- extract(road, p)
        j = j+1
      }
      if(is.na(e2)){
        print("BAD SAMPLE: Simulation fell off the Earth! Skipping...")
        j = 999
      }
    }
    sim[[i]] <- prox.sim
    j = 0
  }
}

par(mfrow=c(1,2))
hist(prox, main="Original Proximity to Roadsides", col="gray", border="white", breaks=20, ylim=c(0,20), xlim=c(0,8000))
hist(unlist(sim), main="Proximity of Simulated Data to Roadsides (Land Class Restricted)", col="gray", border="white", breaks=20, xlim=c(0,8000))

t <- ks.test(prox, unlist(sim))

par(mfrow=c(1,1))
plot(density(prox), main=paste("Density of Restricted Occurrence Record Distances to Roadsides (", o[1,]$Taxon.Scie, ")", sep=""), xlab="Distance from Road (meters)", sub=paste("Two Sample Kolmogorov-Smirnov Test:", t$p.value))
lines(density(unlist(sim)), col="red")
mtext(paste("Sample Mean:", round(mean(prox), digits=2), "  |    Simulation Mean:", round(mean(unlist(sim)), digits=2)), side=4)

## compute summary statistics for this run ##

for(k in 1:length(sim)){numGreater = numGreater + sum(sim[[k]] > prox[k])}
for(k in 1:length(sim)){numCloser = numCloser + sum(sim[[k]] < prox[k])}
for(k in 1:length(sim)){numEqual = numEqual + sum(sim[[k]] == prox[k])}

paste("For this simulation", numGreater, "simulations of", i*iterations, "simulations were farther from roadsides than the original data.")
paste("For this simulation", numCloser, "simulations of", i*iterations, "simulations were closer to roadsides than the original data.")
paste("For this simulation", numEqual, "simulations of", i*iterations, "simulations were equal distance to roadsides than the original data.")


