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


begin=0.25
end=0.75

res = boxplot_RMSE(RMSE,n_data = n_data, model_=c(T,T,F,F), palette = "viridis", begin=0.35, end=0.65)

tmp = gray.colors(2, start = begin,end=end)  
res + scale_fill_grey(start=begin, end=end) + scale_color_manual(values=darken(tmp, amount = 0.25))


### no eval.FEM.time per sol esatta ###
library(grid)
library(gridExtra)

load("data/SpaceTime-ontario-2022-05-14-02_03_58.RData")
source("SpatialTimePlot.R")

time_locations = seq(from=0.,to=pi/2,length.out=6)
library(fdaPDE)
FEM_time = FEM.time(coeff=mean.field.fdaDPE[,4],
                    time_mesh=time_locations,
                    FEMbasis = FEMbasis,
                    FLAG_PARABOLIC = T) 
time_evals = seq(from=0, to=pi/2, length.out=20)

tmp = eval.FEM.time(FEM.time=FEM_time, locations = FEMbasis$mesh$nodes, time.instants = time_evals)
coeff.evals = matrix(tmp, nrow=FEMbasis$nbasis, ncol=length(time_evals))

# fino a riga 50 di SpaceTimeRegression prima di proseguire
x.loc = rep(mesh$nodes[,1], times=length(time_evals))
y.loc = rep(mesh$nodes[,2], times=length(time_evals))
space.locations = cbind(x.loc, y.loc)
time.locations = rep(time_evals, each=nrow(mesh$nodes))
n_time_locs= length(time_evals)
ND.space.only = compute_NetworkDist_matrix(mesh, 1:nnodes)
ND.space.time = matrix(nrow=nnodes*n_time_locs, ncol=nnodes*n_time_locs)

for(i in 1:n_time_locs){
  for( j in 1:n_time_locs){
    ND.space.time[(1 +(i-1)*nnodes):(i*nnodes), 
                  (1 +(j-1)*nnodes):(j*nnodes)] = ND.space.only
  }
}

field = fun.time(func,func.time)(ND.space.time, time_locations, 
                                 source = 63,sigma = 2.5)
field.2 = fun.time(func,func.time)(ND.space.time, time_locations,
                                   source = 1308, sigma = 1.25)
field.3 = fun.time(func,func.time)(ND.space.time, time_locations,
                                   source = 1018,sigma = 2)
field.4 = fun.time(func,func.time)(ND.space.time, time_locations,
                                   source = 1125,sigma = 1.25)
field = field + field.2 + field.3 + field.4

rm(ND.space.only, ND.space.time)
coeff.ex.evals = matrix(field, nrow=FEMbasis$nbasis, ncol=length(time_evals))

dir = "img/esempio/"


col.max = max(coeff.evals, coeff.ex.evals)
col.min = min(coeff.evals, coeff.ex.evals)
line.size = 1
estimate.field = list()
true.field = list()
loop.list = list()

pdf("img/esempio/provaAAA-magma.pdf")
for(i in 1:length(time_evals)){
  title_ = TeX(sprintf("$t = %f", round(time_evals[i],digits = 2)))
  
  estimate.field[[i]] <- R_plot_graph.ggplot2.2(FEM(coeff.evals[,i], FEMbasis),
                                                line.size = line.size,
                                                color.max = col.max,
                                                color.min = col.min,
                                                title = "estimate field",
                                                return.ggplot.object = T,
                                                legend.pos = "right",
                                                palette = magma,
                                                title.size = 16)
  
  true.field[[i]] <- R_plot_graph.ggplot2.2(FEM(coeff.ex.evals[,i], FEMbasis),
                                            line.size = line.size,
                                            color.max = col.max,
                                            color.min = col.min,
                                            title = "true field",
                                            return.ggplot.object = T,
                                            legend.pos = "right",
                                            palette = magma,
                                            title.size = 16)
  
  grid.arrange(true.field[[i]], estimate.field[[i]], 
               ncol=2, nrow=1, top=textGrob(title_, gp=gpar(fontsize=26,font=8)) )
  
}
dev.off()
