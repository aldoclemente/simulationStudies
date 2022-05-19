library(grid)
library(gridExtra)
setwd(paste(getwd(),"/SpaceTimeRegression",sep=""))
load("data/SpaceTime-ontario-2022-05-14-02_03_58.RData")
source("SpatialTimePlot.R")

time_locations = seq(from=0.,to=pi/2,length.out=6)
library(fdaPDE)
FEM_time = FEM.time(coeff=mean.field.fdaDPE[,4],
                    time_mesh=time_locations,
                    FEMbasis = FEMbasis,
                    FLAG_PARABOLIC = T) 
time_evals = seq(from=0, to=pi/2, length.out=50)

tmp = eval.FEM.time(FEM.time=FEM_time, locations = FEMbasis$mesh$nodes, time.instants = time_evals)


coeff.evals = matrix(tmp, nrow=FEMbasis$nbasis, ncol=length(time_evals))

FEM_time = FEM.time(coeff=field,
                    time_mesh=time_locations,
                    FEMbasis = FEMbasis,
                    FLAG_PARABOLIC = T)
tmp = eval.FEM.time(FEM.time=FEM_time, locations = FEMbasis$mesh$nodes, time.instants = time_evals)
coeff.ex.evals = matrix(tmp, nrow=FEMbasis$nbasis, ncol=length(time_evals))

dir = "img/esempio/prova-"


col.max = max(coeff.evals, coeff.ex.evals)
col.min = min(coeff.evals, coeff.ex.evals)
line.size = 1
estimate.field = list()
true.field = list()
loop.list = list()
for( i in 1:length(time_evals)){
  title_ = TeX(sprintf("$t = %f", round(time_evals[i],digits = 2)))
  
   
  estimate.field[[i]] <- R_plot_graph.ggplot2.2(FEM(coeff.evals[,i], FEMbasis),
                                                    line.size = line.size,
                                                    color.max = col.max,
                                                    color.min = col.min,
                                                    title = "estimate field",
                                                    return.ggplot.object = T,
                                                    legend.pos = "right",
                                                    title.size = 16)
  
  true.field[[i]] <- R_plot_graph.ggplot2.2(FEM(coeff.ex.evals[,i], FEMbasis),
                                                line.size = line.size,
                                                color.max = col.max,
                                                color.min = col.min,
                                                title = "true field",
                                                return.ggplot.object = T,
                                                legend.pos = "right",
                                                title.size = 16)

  png(filename = paste(dir,i-1,".png",sep=""))
  grid.arrange(true.field[[i]], estimate.field[[i]], 
               ncol=2, nrow=1, top=textGrob(title_, gp=gpar(fontsize=26,font=8)) )
  dev.off()
  
  #png(filename = paste(dir,i-1,".png",sep=""))
  #print(estimate.field[[i]])
  #dev.off()  
}


pdf("img/esempio/prova.pdf")
for(i in 1:length(time_evals)){
  title_ = TeX(sprintf("$t = %f", round(time_evals[i],digits = 2)))
  grid.arrange(true.field[[i]], estimate.field[[i]], 
               ncol=2, nrow=1, top=textGrob(title_, gp=gpar(fontsize=26,font=8)) )
  
}
dev.off()
