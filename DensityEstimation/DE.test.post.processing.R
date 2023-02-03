#################################################
############## Density Estimation ###############
##############       Test 1       ###############
#############   Post - Processing    ############
#################################################

folder.imgs = paste(folder.name,"img/",sep="")
if(!dir.exists(folder.imgs)) {
  dir.create(folder.imgs)
}
nsim = nrow(rmse.DE_PDE)
RMSE = matrix(nrow = nsim * length(n), ncol=length(methods.names)) 
RMSE[,1] = as.vector(rmse.DE_PDE)
RMSE[,2] = as.vector(rmse.KDE_PDE)
RMSE[,3] = as.vector(rmse.KDE_ES)
RMSE[,4] = as.vector(rmse.KDE_2D)
RMSE[,5] = as.vector(rmse.VORONOI)


pdf(paste(folder.imgs,"RMSE",".pdf",sep=""), width=12)

print(boxplot_RMSE(RMSE, 
             methods = methods,
             methods.names = methods.names, 
             nsim = nsim,
             title.size=26,
             begin=1, #color
             end=0.25,   #color
             width =0.75,
             n = n,
             title="RMSE"))
dev.off()

line.size = 1.
pdf(paste(folder.imgs,"estimates-",line.size,".pdf",sep=""))
estimates = list()
estimates[[1]] = FEM(DE_PDE.FEM, FEMbasis)
estimates[[2]] = FEM(KDE_PDE.FEM, FEMbasis)
estimates[[3]] = FEM(KDE_ES.FEM, FEMbasis)
estimates[[4]] = FEM(KDE_2D.FEM, FEMbasis)
estimates[[5]] = FEM(VORONOI.FEM, FEMbasis)

PLOTS <- plot_estimates(estimates,
                        true.density = FEM(true.density.FEM,FEMbasis),
                        methods = methods,
                        methods.names = methods.names,
                        true.density.main = "Density",
                        palette="viridis",
                        line.size = line.size)

print(PLOTS$density.plot)
for(i in 1:length(estimates)){
  if(methods[i]) print(PLOTS$estimates.plot[[i]])
}

PLOTS <- plot_estimates(estimates, 
                        true.density = FEM(true.density.FEM,FEMbasis),
                        methods = methods,
                        methods.names = methods.names,
                        true.density.main = "Density",
                        palette="magma",
                        line.size = line.size)
print(PLOTS$density.plot)
for(i in 1:length(estimates)){
  if(methods[i])  print(PLOTS$estimates.plot[[i]])
}

PLOTS <- plot_estimates(estimates,
                        true.density = FEM(true.density.FEM,FEMbasis),
                        methods = methods,
                        methods.names = methods.names,
                        true.density.main = "Density",
                        line.size = line.size)

print(PLOTS$density.plot)
for(i in 1:length(estimates)){
  if(methods[i])  print(PLOTS$estimates.plot[[i]])
}

dev.off()

pdf(paste(folder.imgs,"point_pattern",".pdf",sep=""))
for(i in n){
PP = rlpp(i, DENSITY)
plot(mesh, pch=".")
points(PP$data$x, PP$data$y, pch=16, col ="red3")
}
dev.off()

estimates = list()
estimates[[1]] = FEM(DE_PDE.FEM, FEMbasis)
estimates[[2]] = FEM(KDE_PDE.FEM, FEMbasis)
estimates[[3]] = FEM(KDE_ES.FEM, FEMbasis)
estimates[[4]] = FEM(KDE_2D.FEM, FEMbasis)
estimates[[5]] = FEM(VORONOI.FEM, FEMbasis)

PLOTS <- plot_estimates(estimates,
                        true.density = FEM(true.density.FEM,FEMbasis),
                        methods = methods,
                        methods.names = methods.names,
                        true.density.main = "Density",
                        palette="inferno",
                        line.size = line.size)

print(PLOTS$density.plot)
for(i in 1:length(estimates)){
  if(methods[i]) print(PLOTS$estimates.plot[[i]])
}



