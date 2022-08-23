#################################################
############## Density Estimation ###############
##################  Test 1 ######################
#################################################

graphics.off()
rm(list=ls())

setwd("DensityEstimation/")
source("../utils.R")
source("setting.R")

sett = setting("ontario")

mesh = sett$mesh
FEMbasis = sett$FEMbasis
nnodes = sett$nnodes

spat.stat.linnet = sett$spat.stat.linnet

# point pattern #
#n = c(50, 100, 250, 500)
n = 50
set.seed(1234)

#sources = sample(nnodes, 4)

#sources = c(10, 500, 1000)
sources = c(10, 250, 350, 1)

auxiliary_test1 = function(x, y, seg, tp, sigma= 0.375, 
                           nodes.lpp = ppp(x = mesh$nodes[,1], y = mesh$nodes[,2], 
                                           window = owin(xrange = c(min(mesh$nodes[,1]),max(mesh$nodes[,1])),
                                                         yrange = c(min(mesh$nodes[,2]),max(mesh$nodes[,2])))),
                           L = spat.stat.linnet,
                           source = sources)
{ 
  PP = ppp(x = x, y = y, window = nodes.lpp$window)
  ND = crossdist.lpp(lpp(nodes.lpp, L), lpp(PP, L))
  
  return(   0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[1],]^2/(2*sigma^2)) + 
            0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[2],]^2/(2*sigma^2)) + 
            0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[3],]^2/(2*sigma^2)) + 
            0.25 * 1/sqrt(2*pi*sigma^2) * exp(-ND[source[4],]^2/(2*sigma^2))   )
  
  
}


DENSITY = linfun(auxiliary_test1, spat.stat.linnet)

Mass = CPP_get.FEM.Mass.Matrix(FEMbasis)
true.density = DENSITY(x=mesh$nodes[,1], y=mesh$nodes[,2])
true.density = true.density / sum( Mass %*% true.density)

PP = rlpp(n=n, f = DENSITY)
data = cbind(PP$data$x, PP$data$y)
nsim = 30

rmse = matrix(nrow=nsim, ncol = 5)
# DE-PDE #

for( i in 1:nsim){
  
lambda = 10^seq(from=-5, to=-3,length.out = 5)
DE_PDE = fdaPDE::DE.FEM(data = data, FEMbasis = FEMbasis,
                            lambda = lambda,
                            preprocess_method ="SimplifiedCV",
                            nfolds = 5) 

rmse[,1] = sqrt(mean((true.density - exp(DE_PDE$g))^2 ))

plot(FEM(true.density, FEMbasis))
plot(FEM(exp(DE_PDE$g), FEMbasis))

# spatstat returns the INTENSITY function. 
# The estimation of DENSITY or INTENSITY is equivalent, if n is fixed,
# INTENSITY(p) = n DENSITY(p) for all p. 

# See 
# McSwiggan, Greg, Adrian Baddeley, and Gopalan Nair. 
# "Kernel density estimation on a linear network." 
# Scandinavian Journal of Statistics 44.2 (2017): 324-345.

# KDE-PDE  
bw = bw.lppl(X = PP)
KDE_PDE = densityHeat(x = as.lpp(PP), sigma = as.numeric(bw)) 

rmse[,2] = sqrt(mean( (true.density - as.linfun(KDE_PDE)(mesh$nodes[,1], mesh$nodes[,2])/n)^2 ))

# KDE-ES 
bw = bw.lppl(X = PP)
KDE_ES = densityEqualSplit(x = PP, sigma = bw)

rmse[,3] = sqrt(mean( (true.density - as.linfun(KDE_ES)(mesh$nodes[,1], mesh$nodes[,2])/n)^2 ))

# KDE-2D
bw = bw.scott(X = PP)
KDE_2D = densityQuick.lpp(X = PP, sigma = bw) #, at = points)

rmse[,4] = sqrt(mean( (true.density - as.linfun(KDE_2D)(mesh$nodes[,1], mesh$nodes[,2])/n)^2 ) )

# KDE-VORONOI #
bw = bw.voronoi(X = PP) 
KDE_VORONOI = densityVoronoi(X = PP, sigma = bw)
rmse[,5] = sqrt(mean( (true.density - as.linfun(KDE_VORONOI)(mesh$nodes[,1], mesh$nodes[,2])/n)^2 ) )

}

nodes.lpp = ppp(x = mesh$nodes[,1], y = mesh$nodes[,2], 
                window = owin(xrange = c(min(mesh$nodes[,1]),max(mesh$nodes[,1])),
                              yrange = c(min(mesh$nodes[,2]),max(mesh$nodes[,2]))))

ND = crossdist(nodes.lpp, PP)
