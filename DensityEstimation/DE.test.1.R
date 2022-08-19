#################################################
############## Density Estimation ###############
##################  Test 1 ######################
#################################################

setwd("DensityEstimation/")
source("../utils.R")
source("setting.R")

sett = setting("simplenet")

mesh = sett$mesh
FEMbasis = sett$FEMbasis
nnodes = sett$nnodes

spat.stat.linnet = sett$spat.stat.linnet

# point pattern #
n = 100
set.seed(1234)
DENSITY = linfun(auxiliary, spat.stat.linnet)

Mass = CPP_get.FEM.Mass.Matrix(FEMbasis)
true.density = DENSITY(x=mesh$nodes[,1], y=mesh$nodes[,2])
true.density = true.density / sum( Mass %*% true.density)

PP = rlpp(n=n, f = DENSITY)
data = cbind(PP$data$x, PP$data$y)
nsim = 1

rmse = matrix(nrow=nsim, ncol = 5)
# DE-PDE #
lambda = 10^seq(from=-4, to=-3,length.out = 10)
DE_PDE = fdaPDE::DE.FEM(data = data, FEMbasis = FEMbasis,
                            lambda = lambda,
                            preprocess_method = "RightCV",
                            nfolds = 10) 

# Integral = sum( Mass %*% exp(DE_PDE$g) )

rmse[,1] = sqrt(mean((true.density - exp(DE_PDE$g))^2 ))

plot(FEM(true.density, FEMbasis))
plot(FEM(exp(DE_PDE$g), FEMbasis))

# KDE-PDE # DOES NOT WORK 
# bw = bw.lppl(X = PP)
# KDE_PDE = densityHeat.lpp(x = PP, sigma = bw) 

# rmse[,2] = sqrt(mean(true.density - ))

# KDE-ES # return the INTENSITY FUNCTION ...
bw = bw.lppl(X = PP)
KDE_ES = densityEqualSplit(x = PP, sigma = bw, )

ND = crossdist()
rmse[,3] = sqrt(mean( (true.density - as.linfun(KDE_ES)(mesh$nodes[,1], mesh$nodes[,2]))^2 ) )

# KDE-2D
bw = bw.scott(X = PP)
KDE_2D = densityQuick.lpp(X = PP, sigma = bw, at = points)

KDE_2D_density <- KDE_2D / sum(Mass %*% as.linfun(KDE_2D)(mesh$nodes[,1], mesh$nodes[,2])) # Così la normalizzi...

#pixarea <- with(KDE_2D, KDE_2D$xstep * KDE_2D$ystep)
#KDE_2D_density <- KDE_2D * 128 * pixarea

rmse[,4] = sqrt(mean( (true.density - as.linfun(KDE_2D_density)(mesh$nodes[,1], mesh$nodes[,2]))^2 ) )

plot(KDE_2D)
# KDE-VORONOI #
bw = bw.voronoi(X = PP) 
KDE_VORONOI = densityVoronoi(X = PP, sigma = bw)
rmse[,5] = sqrt(mean( (true.density - as.linfun(KDE_VORONOI)(mesh$nodes[,1], mesh$nodes[,2]))^2 ) )



nodes.lpp = ppp(x = mesh$nodes[,1], y = mesh$nodes[,2], 
                window = owin(xrange = c(min(mesh$nodes[,1]),max(mesh$nodes[,1])),
                              yrange = c(min(mesh$nodes[,2]),max(mesh$nodes[,2]))))

ND = crossdist(nodes.lpp, PP)

KDE_2D_density = 1 / n * KDE_2D