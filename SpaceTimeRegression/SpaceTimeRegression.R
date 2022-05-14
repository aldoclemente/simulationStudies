############################
###        Space-Time    ###
###                      ###
###  locations at nodes  ###
############################


setwd(paste(getwd(),"/SpaceTimeRegression",sep=""))
library(plotrix)
source("../utils.R")
source("../SpatialRegression-NoCovariate/NoCovariatesCore.R")

source("settings.R")
source("SpaceTimeCore.R")

func = "exp" # 'exp' 'sin', see settings.R
func.time = "cos"
domain = "ontario" # 'estevan' 'ontario', see settings.R

set = setting(network = domain)

mesh = set$mesh
FEMbasis = set$FEMbasis
nnodes = set$nnodes

mesh.2D = set$mesh.2D
FEMbasis.2D = set$FEMbasis.2D
# from =0 to= pi/den
den = 2 
time_locations = seq(from=0.,to=pi/den,length.out=6) # pi/3 ok
n_time_locs = length(time_locations)
# in ../utils.R
#ND_ = compute_NetworkDist_matrix.time(mesh, 1:nnodes, time_locations)
ED_ = NULL

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
betas= c(1,1)
W = matrix(nrow=nnodes*length(time_locations), ncol=2)
W[,1] = rnorm(nnodes*length(time_locations), mean=0.5,sd=0.25)
W[,2] = 1/4 * fun("sin")(ND.space.only, source=63, sigma=10)

signal = field + W%*%betas
observations = signal + rnorm(length(field), mean=0, sd = 0.05*abs(diff(range(signal))) )

true.field  = matrix(field,
                    nrow = nnodes,
                    ncol = length(time_locations))
true.signal = matrix(signal,
                     nrow = nnodes,
                     ncol = length(time_locations))
                     
{
  library(gridExtra)
  library(latex2exp)
source("../Auxiliary/R_plot_graph.ggplot2.R")  
col.min = min(field)
col.max = max(field)  
ggplot_list = list()
ggplot_list.signal = list()
n_istants = length(time_locations)
for( i in 1:n_time_locs){
  if(i == 1)
    title_ = TeX(sprintf("$t = 0"))
  else if(i==2)
     title_ = TeX(sprintf("$t = \\pi / %d$",den*(n_istants-1)))
  else if(i != n_time_locs)
     title_ = TeX(sprintf("$t =  %d \\pi / %d$",i-1,den*(n_istants-1)))
  else{
     if(den==1){
      title_ = TeX(sprintf("$t =  \\pi$"))
     }else{
       title_ = TeX(sprintf("$t =  \\pi / %d$", den))  
      }
   }
   
    
  #print(
  ggplot_list[[i]] <- R_plot_graph.ggplot2.2(FEM(true.field[,i], FEMbasis),
                                             line.size = 1,
                                             color.max = col.max,
                                             color.min = col.min,
                                             title = title_,
                                             return.ggplot.object = T,
                                             legend.pos = "right") 
  ggplot_list.signal[[i]] <- R_plot_graph.ggplot2.2(FEM(true.signal[,i], FEMbasis),
                                             line.size = 1,
                                             title = title_,
                                             return.ggplot.object = T,
                                             legend.pos = "right") 
  
}  

pdf("prova.pdf", height=9,width = 9)
grid.arrange(grobs = ggplot_list, 
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)))
grid.arrange(grobs = ggplot_list.signal, 
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)))


dev.off()
}
n_sim=30
lambdaS = 10^seq(from=-0.5,to=0.5, length.out=4)
lambdaT = 10^seq(from=-1,to=1, length.out=4)  

ED_=NULL
n_data = c(50, 100, 250, 500)

invisible(capture.output(
  results<- SpaceTimeRegressionCore(
                    ND_ = ND.space.time, 
                    ED_ = ED_, 
                    observations = observations,
                    n_sim = n_sim, 
                    n_data = n_data, 
                    lambdaS = lambdaS,
                    lambdaT = lambdaT,
                    mesh= mesh,
                    FEMbasis = FEMbasis,
                    FEMbasis.2D = FEMbasis.2D, 
                    lambdaS.2D = lambdaS,
                    lambdaT.2D = lambdaT,
                    true.signal = signal,
                    model_=c(T,T,F,F),
                    W = W, betas=betas,
                    time_locations = time_locations)) )

results$tot.time
RMSE = results$RMSE
mean.field.fdaPDE = results$mean.field.fdaPDE

boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F),
             names_ = c("fdaPDE","GWR","",""))

date_ = gsub(":","_",gsub(" ","-",Sys.time()))
head = paste("SpaceTime",domain, date_, sep="-")

filename_ = paste("data/",head,".RData", sep="")
imgfile_  = paste("img/" ,head,".pdf",sep="")

if(!dir.exists("data/")) {
  dir.create("data/")
}

save(RMSE, 
     n_data,
     time_locations,
     observations, 
     signal, 
     field, 
     imgfile_,
     mean.field.fdaPDE,
     FEMbasis, 
     W, betas, file = filename_)

imgfile = paste(paste("plots",domain,func,date_, sep="-"), ".pdf", sep="")

##############################################

setwd(paste(getwd(),"/SpaceTimeRegression",sep=""))
load("data/SpaceTime-ontario-2022-05-14-02_03_58.RData")
source("SpatialTimePlot.R")

if(!dir.exists("img/")) {
  dir.create("img/")
}

SpaceTimePlots(imgfile = imgfile_, 
               time_locations = time_locations,
               true.field = field,            # f 
               true.signal = signal,           # f + W beta 
               mean.field.fdaPDE = mean.field.fdaDPE,
               observations = observations,          # f + W beta + eps
               FEMbasis = FEMbasis,
               n_data = n_data,
               W =W, betas=betas,
               RMSE=RMSE,
               legend.pos.RMSE = "right",
               line.size=1)
