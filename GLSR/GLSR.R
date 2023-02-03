############################
###         GAM          ###
###                      ###
###  locations at nodes  ###
############################

setwd("../GLSR")

library(purrr)
library(fdaPDE)
source("../utils.R")
source("../settings.R")
source("../SpatialRegression-NoCovariate/NoCovariatesCore.R")
source("GLSRCore.R")
#source("plots.R")
source("../Auxiliary/R_plot_graph.ggplot2.R")
#source("settings.R")
{
 pdf("img/networks.pdf")  
 set = setting(network = "estevan")
 #estevan <- R_plot_mesh.ggplot(set$mesh, line.size = 0.65)
 plot(set$mesh, type="n")
 set = setting(network = "ontario")
 #ontario <-R_plot_mesh.ggplot(set$mesh, line.size = 0.65)
 plot(set$mesh, type="n")
 #estevan
 #ontario
 dev.off()
}
func = "exp" # 'exp' 'sin', see settings.R
domain = "ontario" # 'estevan' 'ontario', see settings.R

set = setting(network = domain)

mesh = set$mesh
FEMbasis = set$FEMbasis
nnodes = set$nnodes

mesh.2D = set$mesh.2D
FEMbasis.2D = set$FEMbasis.2D

# family
FAMILY = "poisson" #"binomial"
if(FAMILY=="poisson"){
l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv
}else if(FAMILY == "binomial"){
  logit <- function(x){qlogis(x)}
  inv.logit <- function(x){plogis(x)}
  link = logit
  inv.link = inv.logit
}
# in ../utils.R 
# ND_ = compute_NetworkDist_matrix(mesh, 1:nnodes)
# ED_ = computed_EuclideanDist_matrix(mesh, 1:nnodes)
spat.stat.linnet <- as.spatstat.linnet.fdaPDE(mesh)
ND_ <- spatstat.linnet::crossdist.lpp(X= as.LPP(mesh$nodes, spat.stat.linnet),
                                        Y= as.LPP(mesh$nodes,spat.stat.linnet))
ED_ <- NULL 
set.seed(31415)

# Test 1: ONTARIO, POISSON (func -> EXP)
#sources 63, 1308, 1018, 1125
sample.idx =c(250, 1000, 950, 63)
field = fun(func)(ND_,source = sample.idx[1],sigma = 0.175)
field.2 = fun(func)(ND_,source = sample.idx[2],sigma = 0.175)
field.3 = fun(func)(ND_,source = sample.idx[3],sigma = 0.35)
field.4 = fun(func)(ND_,source = sample.idx[4],sigma = 0.35)

field = field + field.2 + field.3 + field.4
plot(FEM(field, FEMbasis))

W = matrix(nrow=nnodes, ncol=2)
W[,1] = rnorm(nnodes, mean=0.5,sd=0.25)
W[,2] = 1/4 * fun("sin")(ND_, source=63, sigma=1.5)
plot(FEM(W[,2],FEMbasis))
betas = c(1.,1.)

signal = field + W%*%betas
param = signal + rnorm(nnodes, mean=0, sd=0.05*diff(range(signal)))
param = signal 
lambda = 10^seq(from=-3,to=0,length.out=25)
lambda.2D = NULL
####

mu<-inv.link(param)
range(param)
if(FAMILY=="poisson"){
  response <- as.numeric(rpois(nnodes, lambda = mu)) 
  range(response)
}else{
  response <- as.numeric(rbernoulli(nnodes, p = mu))
  range(response)
}

n_data = c(100, 250, 500, 1000)

n_sim=30
invisible(capture.output(
  results<- gamCore(ND_ = ND_, ED_ = ED_, 
                               observations=response,
                               n_sim=n_sim, 
                               n_data= n_data, 
                               lambda=lambda,
                               mesh= mesh,
                               FEMbasis=FEMbasis,
                               FEMbasis.2D=FEMbasis.2D, lambda.2D = lambda,
                               true.signal=signal, 
                               model_=c(T,T,F,F), W = W, betas=betas, 
                               FAMILY=FAMILY) ))


results$tot.time
RMSE = results$RMSE
mean.field.fdaPDE=results$mean.field.fdaPDE
estimates = results$estimates
boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F),
             names_ = c("fdaPDE","GWR","",""))

date_ = gsub(":","_",gsub(" ","-",Sys.time()))

if(domain=="estevan"){
  head = "test-1/"
}else{
  head = "test-2/"
}
date_ = gsub(":","_",gsub(" ","-",Sys.time()))

if(!dir.exists("data/"))
  dir.create("data/")

foldername_ = paste("data/",head, sep = "")
if(!dir.exists(foldername_))
  dir.create(foldername_)

foldername_ = paste(foldername_, date_, "/",sep = "")
dir.create(foldername_)

save(RMSE, 
     n_data, 
     field,  # f  
     signal, # f + W beta
     param,  # f + W beta + eps 
     response, # observations
     mean.field.fdaPDE,
     FEMbasis,
     estimates,
     W,
     betas, 
     foldername_,
     file = paste(foldername_, "data.RData",sep=""))

palette = "magma" # "viridis" "magma
imgfile_ = paste(foldername_,"imgs.pdf",sep="")

if(palette == "ggplot")
  palette=NULL

# # # img with cov # # # 
foldername_ = "data/test-2/2023-02-01-18_18_38/"
load("data/test-2/2023-02-01-18_18_38/data.RData")
source("../SpatialRegression-NoCovariate/NoCovariatesCore.R")
library(fdaPDE)
library(plotrix)
source("GLSRWithCovPlots.R") 
palette = "viridis"

imgfolder_ = paste(foldername_,"imgs/",sep="")
if(!dir.exists(imgfolder_))
  dir.create(imgfolder_)

line.size = 1.
GLRWithCovPlots(imgfile=paste(imgfolder_,"estimates-",line.size,".pdf",sep=""),
              true.field=field,
              true.signal = signal,
              mean.field.fdaPDE = mean.field.fdaPDE,
              param = param,
              observations = response,
              FEMbasis = FEMbasis,
              n_data = n_data,
              W=W, betas=betas,
              RMSE,legend.pos.RMSE = "right",
              palette=palette,
              line.size = line.size)

colors = viridis(n=2, begin=0.95, end=0.25)
pdf(paste(imgfolder_,"RMSE.pdf",sep=""))
boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F), 
             names_ = c("SR-PDE","GWR","",""),
             legend.pos = c(0.825,0.8625), palette=palette,
             colors=colors)
dev.off()