rm(list=ls())
library(fdaPDE)
library(spatstat)
setwd("SR-chicago-case-study/")

source("refine1D.R")
source("../utils.R")
source("../DensityEstimation/utils.R")
data("chicago")

vertices = cbind(chicago$domain$vertices$x, chicago$domain$vertices$y)
edges = cbind( chicago$domain$from, chicago$domain$to)

mesh = create.mesh.1.5D(nodes = vertices, edges = edges)

mesh = normalize_mesh(mesh)
LN = as.spatstat.linnet.fdaPDE(mesh)

x.m = mean(chicago$domain$vertices$x)
y.m = mean(chicago$domain$vertices$y)

x.sd = sd(chicago$domain$vertices$x)
y.sd = sd(chicago$domain$vertices$y)

x.norm = (chicago$data$x - x.m)/x.sd
y.norm = (chicago$data$y - y.m)/y.sd

chicago.norm = spatstat.linnet::lpp(X= ppp(x.norm, y=y.norm, 
                                        window= LN$window), 
                                 L= LN)
#chicago.norm = list(data= data.norm, domain= LN )

DATA = chicago.norm # chicago

{
x11()
plot(mesh, pch=".")
points(DATA$data$x, DATA$data$y, col = "red", pch=16, cex=1.5)
}
### FAMILY - poisson ###
delta = 0.25

mesh = refine.mesh.1.5D(mesh, delta=delta)
FEMbasis = create.FEM.basis(mesh)
#new_to_old  = refine1D(mesh$nodes, mesh$edges, delta)$new_to_old

#setting regions 
nnodes = nrow(mesh$nodes)
nregion = 12 #8 #6 #10 
nedges = nrow(mesh$edges)
ndata = length(DATA$data$x)

PP = runiflpp(n=nregion/2, L = chicago.norm$domain) 
centroids_ = cbind(PP$data$x, PP$data$y)
centroids_ = projection.points.1.5D(mesh=mesh, location=centroids_)

# lines == true edges of the network / edges == edges of the discretized network
data_ = cbind(DATA$data$x, DATA$data$y)
result_ <- kmeans(x=data_, centers=nregion, iter.max = 100)
centroids_ = rbind( centroids_, projection.points.1.5D(mesh, locations= result_$centers) )

{
x11()
plot(mesh, pch=".")
points(data_, col="red", pch=16, cex=2)
points(result_$centers, col="blue", pch = 16, cex=2)
points(centroids_, col="green4", pch=16, cex=2)
legend("topright", legend=c("data", "2D", "1.5D"), 
        col=c("red", "blue", "green4"), pch = 16, cex=1.25 )
}

lines_to_region <- set_region(centroids_, mesh=mesh, LN = DATA)

nregion = nrow(centroids_)

incidence_matrix = matrix(0, nrow=nregion, ncol=nedges)

# for(i in 1:nrow(new_to_old)){
#   incidence_matrix[ lines_to_region[new_to_old[i]],i] = 1
# }

for(i in 1:nedges){
   incidence_matrix[ lines_to_region[i],i] = 1
}

mask_= matrix(0, nrow=nregion, ncol=1)
for( i in 1:nregion){
  if( !sum(incidence_matrix[i,])){
    mask_[i] = 1
    print(mask_[i])
  }
}

if( sum(mask_)){
  incidence_matrix = incidence_matrix[-mask_,]
  nregion = nregion - sum(mask_)
  centroids_ = centroids_[-mask_,]
  lines_to_region <- set_region(centroids_, mesh=mesh, LN = PP)
}


response = rep(0, times= nregion)

for( i in DATA$data$seg){
  response[lines_to_region[i]] = response[lines_to_region[i]] + 1
}
range(response)

{
x11()
plot_region(lines_to_region, response, LN=DATA,mesh=mesh)
}
lambda_GSR = 10^seq(from=-2,to=1,length.out=250)

GSR_PDE <- smooth.FEM(observations = response,
                         covariates = NULL,
                         FEMbasis = FEMbasis,
                         incidence_matrix = incidence_matrix,
                         lambda = lambda_GSR,
                         lambda.selection.criterion = "grid",
                         lambda.selection.lossfunction = "GCV",
                         DOF.evaluation = "exact",
                         family="poisson")

lambda_opt <- GSR_PDE$optimization$lambda_position
# sqrt(GCV_ * (sum(!is.na(observations)) - dof)/sum(!is.na(observations)))
{
x11()
plot(log10(lambda), GSR_PDE$optimization$GCV_vector, xlab="log10(lambda)", ylab="GCV")
}

Mass = fdaPDE::CPP_get.FEM.Mass.Matrix(FEMbasis)
coeff_ <- exp(GSR_PDE$solution$f[,lambda_opt])
coeff_ <- coeff_ / sum( Mass%*% coeff_ )
plot(FEM(coeff_, FEMbasis))
range(coeff_)

lambda_DE = 10^seq(from=-3.5, to=-2.5,length.out = 20)
DE_PDE = fdaPDE::DE.FEM(data = cbind(DATA$data$x, DATA$data$y), FEMbasis = FEMbasis,
                        lambda = lambda_DE,
                        preprocess_method ="RightCV",
                        nfolds = 10)
x11()
plot(log10(lambda_DE), DE_PDE$CV_err, xlab="log10(lambda)", ylab="CV")

sum(Mass%*% exp(DE_PDE$g))
plot(FEM(exp(DE_PDE$g), FEMbasis))
range(exp(DE_PDE$g))

sum( Mass%*% (exp(DE_PDE$g) - coeff_)^2 )

################################################################################
K = 10

dataList = set_Kfold_data(DATA$data, seed = 0)
n = nrow(chicago$data)

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


CV_errors = matrix(0, nrow = K, ncol = 2)

for(i in 1:K){
  
  tmp = get_Kfold_data(dataList, iter = i)
  train_data = tmp$train_data
  test_data = tmp$test_data
  
  # DE-PDE 
  DE_PDE = fdaPDE::DE.FEM(data = cbind(train_data$x, train_data$y), FEMbasis = FEMbasis,
                        lambda = lambda_DE,
                        preprocess_method ="RightCV",
                        nfolds = 10)
  DE_PDE.FEM = FEM(coeff= exp(DE_PDE$g), FEMbasis= FEMbasis )
  
  CV_errors[i,1] = cv_error(FEM = DE_PDE.FEM, 
                            R0 = Mass, data.k = test_data)

  # GSR-PDE 
  response = rep(0, times= nregion)
  for(k in train_data$seg){
    response[lines_to_region[k]] = response[lines_to_region[k]] + 1
  }

  GSR_PDE <- smooth.FEM(observations = response,
                         covariates = NULL,
                         FEMbasis = FEMbasis,
                         incidence_matrix = incidence_matrix,
                         lambda = lambda_GSR,
                         lambda.selection.criterion = "grid",
                         lambda.selection.lossfunction = "GCV",
                         DOF.evaluation = "exact",
                         family="poisson")

  lambda_opt <- GSR_PDE$optimization$lambda_position
  
  coeff_ <- exp(GSR_PDE$solution$f[,lambda_opt])
  coeff_ <- coeff_ / sum( Mass%*% coeff_ )

  GSR_PDE.FEM = FEM(coeff_, FEMbasis) 
  CV_errors[i,2] = cv_error(FEM = GSR_PDE.FEM, 
                            R0 = Mass, data.k = test_data)


}

save(CV_errors, date_, folder.name,
     file = paste(folder.name, "CV_error.RData", sep=""))

boxplot(CV_errors, main = "CV error")
