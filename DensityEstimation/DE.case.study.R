#################################################
############## Density Estimation ###############
##############   Chicago Crimes   ###############
#################################################
graphics.off()
rm(list=ls())
library(fdaPDE)
setwd("DensityEstimation/")
source("../utils.R")
source("utils.R")
source("setting.R")

data("chicago")

{
vertices = cbind(chicago$domain$vertices$x, chicago$domain$vertices$y)
edges = cbind( chicago$domain$from, chicago$domain$to)

mesh = create.mesh.1.5D(nodes = vertices, edges = edges)

#mesh = normalize_mesh(mesh)
res <- normalize_mesh_unit(mesh)
mesh <- res$mesh
LN = as.spatstat.linnet.fdaPDE(mesh)

#x.m = mean(chicago$domain$vertices$x)
#y.m = mean(chicago$domain$vertices$y)
#x.sd = sd(chicago$domain$vertices$x)
#y.sd = sd(chicago$domain$vertices$y)
#x.norm = (chicago$data$x - x.m)/x.sd
#y.norm = (chicago$data$y - y.m)/y.sd

x.norm = (chicago$data$x - res$x.min)/(res$x.max-res$x.min)
y.norm = (chicago$data$y - res$y.min)/(res$y.max-res$y.min)

chicago.norm = spatstat.linnet::lpp(X= ppp(x.norm, y=y.norm, 
                                        window= LN$window), 
                                 L= LN)
}

#spat.stat.linnet = chicago$domain
spat.stat.linnet = chicago.norm
mesh = as.fdaPDE.spatstat.linnet(spat.stat.linnet$domain)
FEMbasis = create.FEM.basis(mesh)
Mass = CPP_get.FEM.Mass.Matrix(FEMbasis)

K = 10

dataList = set_Kfold_data(spat.stat.linnet$data, seed = 0)
n = nrow(spat.stat.linnet$data)

date_ = gsub(":","_",gsub(" ","-",Sys.time()))
if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists(paste("data/chicago/",sep=""))){
  dir.create(paste("data/chicago/",sep=""))
}

folder.name = paste("data/chicago/", date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}


CV_errors = matrix(0, nrow = K, ncol = 4)

for(i in 1:K){
  
  tmp = get_Kfold_data(dataList, iter = i)
  train_data = tmp$train_data
  test_data = tmp$test_data

  # DE-PDE 
  #lambda = 10^seq(from=4.45, to=4.75,length.out = 20)
  lambda = 10^seq(from=-3.5, to=-2.5,length.out = 20)
  DE_PDE = fdaPDE::DE.FEM(data = cbind(train_data$x, train_data$y), FEMbasis = FEMbasis,
                          lambda = lambda,
                          preprocess_method ="RightCV",
                          nfolds = 10)
  
  CV_errors[i,1] = cv_error(FEM = FEM(coeff= exp(DE_PDE$g), FEMbasis), 
                            R0 = Mass, data.k = test_data)
    
  # Training Point Pattern over spat.stat.linnet Road Network
  PP_train = lpp(X = ppp(x = train_data$x, y = train_data$y, window = spat.stat.linnet$domain$window),
                 L = spat.stat.linnet$domain)
  # KDE-PDE
  bw = bw.lppl(X = PP_train)
  KDE_PDE = densityHeat(x = PP_train, sigma = as.numeric(bw)) 
  
  CV_errors[i,2] = cv_error(FEM = FEM(coeff=as.linfun(KDE_PDE)(mesh$nodes[,1], mesh$nodes[,2])/n, FEMbasis), 
                            R0 = Mass, data.k = test_data)
  
  # KDE-2D
  bw = bw.scott(X = PP_train)
  KDE_2D = densityQuick.lpp(X = PP_train, sigma = bw) #, at = points)
  
  CV_errors[i,3] = cv_error(FEM = FEM(coeff=as.linfun(KDE_2D)(mesh$nodes[,1], mesh$nodes[,2])/n, FEMbasis), 
                            R0 = Mass, data.k = test_data)
  
  # KDE-VORONOI #
  bw = bw.voronoi(X = PP_train) 
  KDE_VORONOI = densityVoronoi(X = PP_train, sigma = bw)
  
  CV_errors[i,4] = cv_error(FEM = FEM(coeff=as.linfun(KDE_VORONOI)(mesh$nodes[,1], mesh$nodes[,2])/n, FEMbasis), 
                            R0 = Mass, data.k = test_data)
  
}

date_ = gsub(":","_",gsub(" ","-",Sys.time()))

save(CV_errors, date_, folder.name,
     file = paste(folder.name, "CV_error.RData", sep=""))

# Competing methods
# DE-PDE
DE_PDE = fdaPDE::DE.FEM(data = cbind(spat.stat.linnet$data$x, spat.stat.linnet$data$y), FEMbasis = FEMbasis,
                        lambda = lambda,
                        preprocess_method ="RightCV",
                        nfolds = 10)

DE_PDE.FEM = FEM(coeff= exp(DE_PDE$g), FEMbasis)

# Training Point Pattern over Chicago Road Network
PP = lpp(X = ppp(x = spat.stat.linnet$data$x, y = spat.stat.linnet$data$y, window = spat.stat.linnet$domain$window),
               L = spat.stat.linnet$domain)
# KDE-PDE
bw = bw.lppl(X = PP)
KDE_PDE = densityHeat(x = PP, sigma = as.numeric(bw)) 

KDE_PDE.FEM = FEM(coeff=as.linfun(KDE_PDE)(mesh$nodes[,1], mesh$nodes[,2])/n, FEMbasis)

# KDE-2D
bw = bw.scott(X = PP)
KDE_2D = densityQuick.lpp(X = PP, sigma = bw) #, at = points)

KDE_2D.FEM = FEM(coeff=as.linfun(KDE_2D)(mesh$nodes[,1], mesh$nodes[,2])/n, FEMbasis)

# KDE-VORONOI #
bw = bw.voronoi(X = PP) 
KDE_VORONOI = densityVoronoi(X = PP, sigma = bw)

KDE_VORONOI.FEM = FEM(coeff=as.linfun(KDE_VORONOI)(mesh$nodes[,1], mesh$nodes[,2])/n, FEMbasis)

save(DE_PDE.FEM, KDE_PDE.FEM, KDE_2D.FEM, KDE_VORONOI.FEM,
     spat.stat.linnet,
     file = paste(folder.name,"estimates.RData",sep=""))

source("DE.case.study.post.processing.R")


