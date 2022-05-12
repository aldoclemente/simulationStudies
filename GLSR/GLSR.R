############################
###         GAM          ###
###                      ###
###  locations at nodes  ###
############################

setwd("C:/Users/Aldo/Documents/SimulationStudies/GLSR")
library(plotrix)
library(purrr)
source("../utils.R")
source("../SpatialRegression-NoCovariate/NoCovariatesCore.R")
source("GLSRCore.R")
#source("plots.R")
source("../Auxiliary/R_plot_graph.ggplot2.R")
source("settings.R")

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
ND_ = compute_NetworkDist_matrix(mesh, 1:nnodes)
ED_ = computed_EuclideanDist_matrix(mesh, 1:nnodes)
set.seed(31415)

# Test 1: ONTARIO, POISSON (func -> EXP)
#sources 63, 1308, 1018, 1125

field = fun(func)(ND_,source = 63,sigma = 2)
field.2 = fun(func)(ND_,source = 1308,sigma = 1.25)
field.3 = fun(func)(ND_,source = 1018,sigma = 2)
field.4 = fun(func)(ND_,source = 1125,sigma = 1.25)

field = field + field.2 + field.3 + field.4
plot(FEM(field, FEMbasis))

W = matrix(nrow=nnodes, ncol=2)
W[,1] = rnorm(nnodes, mean=0.5,sd=0.25)
W[,2] = 1/4 * fun("sin")(ND_, source=63, sigma=10)
plot(FEM(W[,2],FEMbasis))
betas = c(1.,1.)

signal = field + W%*%betas
param = signal + rnorm(nnodes, mean=0, sd=0.05*diff(range(signal)))
lambda = 10^seq(from=0,to=2,length.out=25)
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
boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F),
             names_ = c("fdaPDE","GWR","",""))


date_ = gsub(":","_",gsub(" ","-",Sys.time()))

if(domain=="estevan"){
  head = "test-1"
}else{
  head = "test-2"
}
date_ = gsub(":","_",gsub(" ","-",Sys.time()))
tail_ = paste(domain, func, date_, sep="-")

filename_ = paste(paste(paste("data/",head,sep=""),tail_, sep="-"), ".RData", sep="")
imgfile_  = paste(paste(paste(paste("img/" ,head,sep=""),"plots",sep="-"),tail_,sep="-"),".pdf",sep="")

save(RMSE, 
     n_data, 
     field,  # f  
     signal, # f + W beta
     param,  # f + W beta + eps 
     response, # observations
     mean.field.fdaPDE,
     W,
     betas, 
     file = filename_)

if(is.null(W)){ 
source("GLRNoCovPlots.R")
GLRNoCovPlots(imgfile,
             field=param, line.size.field = 0.5,    
             response, 
             RMSE,legend.pos.RMSE = c(0.85, 0.85))
}else{
# # # img with cov # # # 
source("GLSRWithCovPlots.R") # è esattamente la stessa... interessant
GLRWithCovPlots(imgfile=imgfile_,
              true.field=field,
              true.signal = signal,
              mean.field.fdaPDE = mean.field.fdaPDE,
              param = param,
              observations = response,
              FEMbasis = FEMbasis,
              n_data = n_data,
              W=W, betas=betas,
              RMSE,legend.pos.RMSE = "right",
              line.size = 1)
}

source("../SpatialRegression-Covariates/RegressionWithCovPlots.R")

RegressionWithCovPlots(imgfile="img/prova.pdf",
                true.field=field,
                true.signal = signal,
                mean.field.fdaPDE = mean.field.fdaPDE,
                #param = param,
                observations = response,
                FEMbasis = FEMbasis,
                n_data = n_data,
                W=W, betas=betas,
                RMSE,legend.pos.RMSE = "right",
                line.size = 1)
