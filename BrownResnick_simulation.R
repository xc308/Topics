#====================
# Max-stable training
#====================

rm(list = ls())

install.packages("lattice")
install.packages("dqrng")
install.packages("SpatialExtremes")

library(lattice)
library(dqrng)
library(SpatialExtremes)

source("Algorithms.R")
rnorm <- dqrnorm

# smooth parameter: difficult to estimate
  # but the convergence quantile is easier
   

#==========================================
# Simulation of 1 Std. Brown-Resnick process
# on dense grid on [0, 1]-interval
#==========================================

coord <- seq(0,1,0.004) # 251 
vario <- function(h){abs(h)}

Z <- simu_extrfcts(no.simu=1, coord=coord, vario=vario, type="brownresnick") 
plot(coord,Z$res,type="l")


#-------------------
# 10 such processes
#-------------------

Z <- simu_extrfcts(no.simu = 10, coord = coord, vario = vario, type = "brownresnick")
matplot(coord, t(Z$res), type  = "l") # matrix plot
# value can only be positive
str(Z$res) # num [1:10, 1:251]

## from Frechet to Gumbel
matplot(coord, t(log(Z$res)), type = "l")
# value can be both positive and negative



#=====================
# 2-D grid simulation
#=====================

xseq <- yseq <- seq(0,5,0.2)
coord <- as.matrix(expand.grid(xseq, yseq))
vario <- function(h) {
  sqrt(sum(h^2))
}


#-------
# 1-D standard Brown Resnick process on grid [0, 5]*[0, 5]
#-------

Z <- simu_extrfcts(no.simu = 1, coord = coord, vario = vario, type = "brownresnick")
data <- data.frame(x = coord[, 1], y = coord[, 2], z = t(log(Z$res)))
levelplot(z ~ x * y, data)

data <- data.frame(x=coord[,1],y=coord[,2],z=t(log(Z$res)))
levelplot(z ~ x * y, data)


#----------
# 6 Standard Brown Resnick processes on grid on [0, 5]*[0, 5]
#----------

Z <- simu_extrfcts(no.simu = 6, coord = coord, vario = vario, type = "brownresnick")
str(Z$res) # num [1:6, 1:676]
dim(coord) # [1] 676   2

data <- data.frame(x = rep(coord[, 1], 6), y = rep(coord[, 2], 6),
           z = as.vector(t(log(Z$res))),
           no.simu.factor = rep(1:6, each = dim(coord)[1]))

levelplot(z ~ x * y | factor(no.simu.factor), data)


















