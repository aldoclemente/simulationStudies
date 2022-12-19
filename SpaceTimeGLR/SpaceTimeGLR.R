############################
###   Space-Time GLR     ###
###                      ###
###  locations at nodes  ###
############################


setwd("SpaceTimeGLR/")
library(plotrix)
source("../utils.R")
source("../SpatialRegression-NoCovariate/NoCovariatesCore.R")

source("settings.R")
source("SpaceTimeGLRCore.R")
FAMILY = "poisson"
if(FAMILY=="poisson"){
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
}

func = "exp" # 'exp' 'sin', see settings.R
func.time = "cos"
domain = "ontario" # 'estevan' 'ontario', see settings.R

set = setting(network = domain)

mesh = set$mesh
FEMbasis = set$FEMbasis
nnodes = set$nnodes

# from =0 to= pi/den
den = 2 
time_locations = seq(from=0.,to=pi/den,length.out=6) # pi/3 ok
n_time_locs = length(time_locations)

ND.space.only = compute_NetworkDist_matrix(mesh, 1:nnodes)
ND.space.time = matrix(nrow=nnodes*n_time_locs, ncol=nnodes*n_time_locs)

for(i in 1:n_time_locs){
  for( j in 1:n_time_locs){
    ND.space.time[(1 +(i-1)*nnodes):(i*nnodes), 
                  (1 +(j-1)*nnodes):(j*nnodes)] = ND.space.only
  }
}
# data 
x.loc = rep(mesh$nodes[,1], times=length(time_locations))
y.loc = rep(mesh$nodes[,2], times=length(time_locations))
space.locations = cbind(x.loc, y.loc)
time.locations = rep(time_locations, each=nrow(mesh$nodes))

space.time.locations = cbind(space.locations, time.locations)

field = fun.time(func,func.time)(ND.space.time, time_locations, 
                                 source = 63,sigma = 2.5)
field.2 = fun.time(func,func.time)(ND.space.time, time_locations,
                                   source = 1308, sigma = 1.25)
field.3 = fun.time(func,func.time)(ND.space.time, time_locations,
                                   source = 1018,sigma = 2)
field.4 = fun.time(func,func.time)(ND.space.time, time_locations,
                                   source = 1125,sigma = 1.25)
field = field + field.2 + field.3 + field.4

# Test 1: ONTARIO (func -> EXP) 
#betas= c(1,1)
#W = matrix(nrow=nnodes*length(time_locations), ncol=2)
#W[,1] = rnorm(nnodes*length(time_locations), mean=0.5,sd=0.25)
#W[,2] = 1/4 * fun("sin")(ND.space.only, source=63, sigma=10)
W = NULL
signal = field  #+ W%*%betas
param = signal 

n_sim=30
lambdaS = 10^seq(from=-0.5,to=0.5, length.out=4)
lambdaT = 10^seq(from=-1.5,to=0.5, length.out=4)  

mu<-inv.link(param)
range(param)

response <- as.numeric(rpois(nnodes*n_time_locs, lambda = mu)) 
range(response)

ED_=NULL
n_data = c(50, 100, 250, 500)

RMSE.fdaPDE  = matrix(0,nrow=n_sim, ncol=length(n_data))
mean.field.fdaPDE = matrix(0,nrow=nnodes*length(time_locations), ncol=length(n_data))

start.time = Sys.time()
for(j in 1:length(n_data)){  
  for(i in 1:n_sim){
    
    tmp.sample_ = sample(x=(1:nnodes), size=n_data[j])
    sample_ = vector(mode="integer")
    for(k in 1:length(time_locations)){
      sample_ = cbind(sample_, tmp.sample_ + (k-1)*nnodes)
    }
    sample_ = as.vector(sample_)
    
    locations = mesh$nodes[tmp.sample_,]
    
    observations_ = response[sample_]
    
    if(is.null(W)){
      ### fdaPDE ### 
      output_CPP = smooth.FEM.time(observations = matrix(observations_,
                                                         nrow = n_data[j],
                                                         ncol = length(time_locations)),
                                   time_locations = time_locations,                           
                                   locations = locations,
                                   FEMbasis = FEMbasis,
                                   lambdaS = lambdaS,
                                   lambdaT = lambdaT,
                                   lambda.selection.criterion = "grid",
                                   lambda.selection.lossfunction = "GCV",
                                   DOF.evaluation = "stochastic",
                                   family = FAMILY) 
      coeff = output_CPP$fit.FEM.time$coeff[,output_CPP$bestlambda[1],output_CPP$bestlambda[2]]
      f = FEM.time(coeff=array(coeff),
                   time_mesh=output_CPP$fit.FEM.time$mesh_time, 
                   FEMbasis=output_CPP$fit.FEM.time$FEMbasis, 
                   FLAG_PARABOLIC=output_CPP$fit.FEM.time$FLAG_PARABOLIC)
      prediction = eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                                 locations=mesh$nodes,
                                 time.instants=time_locations,
                                 lambdaS = output_CPP$bestlambda[1], 
                                 lambdaT = output_CPP$bestlambda[2])
      RMSE.fdaPDE[i,j] = rmse(field, prediction)
      
      mean.field.fdaPDE[,j] = mean.field.fdaPDE[,j] + eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                                                                    locations=mesh$nodes,
                                                                    time.instants=time_locations,
                                                                    lambdaS = output_CPP$bestlambda[1],
                                                                    lambdaT = output_CPP$bestlambda[2])/n_sim
    }
    else{
      output_CPP = smooth.FEM.time(observations = matrix(observations_,
                                                         nrow = n_data[j],
                                                         ncol = length(time_locations)),
                                   time_locations = time_locations,
                                   covariates = W[sample_,],
                                   locations = locations,
                                   FEMbasis = FEMbasis,
                                   lambdaS = lambdaS,
                                   lambdaT = lambdaT,
                                   lambda.selection.criterion = "grid",
                                   lambda.selection.lossfunction = "GCV",
                                   DOF.evaluation = "stochastic",
                                   family = FAMILY)
      prediction = 
        eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                      locations=locations,
                      time.instants=time_locations,
                      lambdaS = output_CPP$optimization$lambda_solution[1],
                      lambdaT = output_CPP$optimization$lambda_solution[2]) +
        W[sample_,]%*%as.vector(output_CPP$beta)
      RMSE.fdaPDE[i,j] = rmse(true.signal[sample_], prediction)
      
      mean.field.fdaPDE[,j] = mean.field.fdaPDE[,j] + eval.FEM.time(FEM.time = output_CPP$fit.FEM.time, 
                                                                    locations=mesh$nodes,
                                                                    time.instants=time_locations,
                                                                    lambdaS = output_CPP$optimization$lambda_solution[1],
                                                                    lambdaT = output_CPP$optimization$lambda_solution[2])/n_sim
      
    }
    
  }
}
tot.time = difftime(Sys.time(), start.time)

RMSE = list( RMSE.fdaPDE  = RMSE.fdaPDE)

results = list(RMSE = RMSE, tot.time=tot.time, mean.field.fdaPDE=mean.field.fdaPDE)

results$tot.time
RMSE = results$RMSE
mean.field.fdaPDE = results$mean.field.fdaPDE

#boxplot_RMSE(RMSE, n_data, model_ = c(T,F,F,F),
#             names_ = c("fdaPDE","GWR","",""))

date_ = gsub(":","_",gsub(" ","-",Sys.time()))
head = paste("SpaceTimeGSR",domain, date_, sep="-")

filename_ = paste("data/",head,".RData", sep="")
imgfile_  = paste("img/" ,head,".pdf",sep="")

if(!dir.exists("data/")) {
  dir.create("data/")
}

save(RMSE, 
     n_data,
     time_locations,
     response, 
     signal, 
     field, 
     imgfile_,
     mean.field.fdaPDE,
     FEMbasis, 
     W, betas, file = filename_)

imgfile = paste(paste("plots",domain,func,date_, sep="-"), ".pdf", sep="")

##############################################

setwd(paste("../SpaceTimeGLR",sep=""))

if(!dir.exists("img/")) {
  dir.create("img/")
}

palette = "viridis" # "viridis" "magma
imgfile_ = paste("img/SpaceTimeGLR-",palette,".pdf",sep="")

if(palette == "ggplot")
  palette=NULL


rmse <- boxplot_RMSE(RMSE=RMSE,
             n_data=n_data,
             model=c(T,F,F,F),names_ = c("GST-PDE","","",""),
             palette=palette,begin=0.25,end=0.25)

pdf(imgfile_)
print(rmse)
dev.off()

source("../SpaceTimeRegression/SpaceTimePlotLoop.R")
imgfile.loop = paste("img/SpaceTimeGLR-Loop-",palette,".pdf",sep="")
SpaceTimeLoop(imgfile = imgfile.loop,
              FEMbasis = FEMbasis,
              time_locations = time_locations,
              field = field,
              mean.field.fdaPDE = mean.field.fdaPDE,
              palette =palette,
              line.size=0.5)
colors = viridis(n=1, begin=0.25, end=0.25)
pdf("img/RMSE-ST-PDE.pdf")
boxplot_RMSE(RMSE, n_data, model_ = c(T,F,F,F), 
             names_ = c("SR-PDE","","",""),
             legend.pos = c(0.825,0.8625), palette=palette,
             colors=colors )
dev.off()

imgfile.6plots = paste("img/SpaceTimeGLR-6plots-", palette,".pdf",sep="")
SpaceTime6Plots(imgfile = imgfile.6plots,
                FEMbasis = FEMbasis,
                time_locations = time_locations,
                field = field,
                mean.field.fdaPDE = mean.field.fdaPDE,
                palette =palette,
                line.size=0.5)



#############################################################

den = 2 
time_locations = seq(from=0.,to=pi/den,length.out=12) # pi/3 ok
n_time_locs = length(time_locations)

ND.space.only = compute_NetworkDist_matrix(mesh, 1:nnodes)
ND.space.time = matrix(nrow=nnodes*n_time_locs, ncol=nnodes*n_time_locs)

for(i in 1:n_time_locs){
  for( j in 1:n_time_locs){
    ND.space.time[(1 +(i-1)*nnodes):(i*nnodes), 
                  (1 +(j-1)*nnodes):(j*nnodes)] = ND.space.only
  }
}
# data 
x.loc = rep(mesh$nodes[,1], times=length(time_locations))
y.loc = rep(mesh$nodes[,2], times=length(time_locations))
space.locations = cbind(x.loc, y.loc)
time.locations = rep(time_locations, each=nrow(mesh$nodes))

space.time.locations = cbind(space.locations, time.locations)

field = fun.time(func,func.time)(ND.space.time, time_locations, 
                                 source = 63,sigma = 2.5)
field.2 = fun.time(func,func.time)(ND.space.time, time_locations,
                                   source = 1308, sigma = 1.25)
field.3 = fun.time(func,func.time)(ND.space.time, time_locations,
                                   source = 1018,sigma = 2)
field.4 = fun.time(func,func.time)(ND.space.time, time_locations,
                                   source = 1125,sigma = 1.25)
field = field + field.2 + field.3 + field.4

# Test 1: ONTARIO (func -> EXP) 
#betas= c(1,1)
#W = matrix(nrow=nnodes*length(time_locations), ncol=2)
#W[,1] = rnorm(nnodes*length(time_locations), mean=0.5,sd=0.25)
#W[,2] = 1/4 * fun("sin")(ND.space.only, source=63, sigma=10)
W = NULL
signal = field  #+ W%*%betas
param = signal 

mu<-inv.link(param)
range(param)

response <- as.numeric(rpois(nnodes*n_time_locs, lambda = mu)) 
range(response)
n_data=100
tmp.sample_ = sample(x=(1:nnodes), size=n_data)
sample_ = vector(mode="integer")
for(k in 1:length(time_locations)){
  sample_ = cbind(sample_, tmp.sample_ + (k-1)*nnodes)
}
sample_ = as.vector(sample_)

locations = mesh$nodes[tmp.sample_,]

observations_ = response[sample_]

estimates = list()
estimates$locations = locations
estimates$response = observations_

estimates$param = param[sample_]
estimates$param = matrix(estimates$param, nrow=n_data, ncol=length(time_locations))
estimates$response = matrix(estimates$response, nrow=n_data, ncol=length(time_locations))

min.col = min(estimates$param)
max.col = max(estimates$param)

library(viridis)
library(ggplot2)
library(latex2exp)
library(grid)
library(gridExtra)
source("../Auxiliary/R_plot_graph.ggplot2.R")
source("../SpatialRegression-NoCovariate/NoCovariatesCore.R")
p=viridis
estimates.obs = list()
estimates.response = list()
pdf("Space-Time-GLR-observations-time.pdf")

for(i in 1:length(time_locations)){

  estimates.obs[[i]] <- R_plot_mesh.ggplot(mesh = mesh,
                                    points_ = estimates$locations,
                                    mu = estimates$param[,i], 
                                    color.min = min.col,
                                    color.max = max.col,
                                    line.size=1,
                                    palette=p,
                                    title = bquote(g(mu[ij]) == f(bold(p)[i], t_j)) )

  print(estimates.obs[[i]])
}
dev.off()


min.col = min(estimates$response)
max.col = max(estimates$response)

pdf("Space-Time-GLR-response-time.pdf")
for(i in 1:length(time_locations)){
estimates.response[[i]] <- R_plot_mesh.ggplot(mesh = mesh,
                                              points_ = estimates$locations,
                                              mu = estimates$response[,i], 
                                              color.min = min.col,
                                              color.max = max.col,
                                              line.size=1,
                                              palette=p,
                                              title = bquote(y[ij]) )


print(estimates.response[[i]])
}
dev.off()
