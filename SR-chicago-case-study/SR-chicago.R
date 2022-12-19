library(fdaPDE)
library(spatstat)
setwd("SR-chicago-case-study/")

source("refine1D.R")
data("chicago")

vertices = cbind(chicago$domain$vertices$x, chicago$domain$vertices$y)
edges = cbind( chicago$domain$from, chicago$domain$to)

mesh = create.mesh.1.5D(nodes = vertices, edges = edges)

plot(mesh, pch=".")
points(chicago$data$x, chicago$data$y, col = "red", pch=16)

### FAMILY - poisson ###
delta = 65
new_to_old  = refine1D(mesh$nodes, mesh$edges, delta)$new_to_old
mesh = refine.mesh.1.5D(mesh, delta=delta)
FEMbasis = create.FEM.basis(mesh)
incidence_matrix = matrix(0, nrow=chicago$domain$lines$n, ncol=nrow(mesh$edges))

for(i in 1:nrow(new_to_old)){
  incidence_matrix[ new_to_old[i],i] = 1
}

response = rep(0, times= chicago$domain$lines$n)

for( i in chicago$data$seg){
  response[i] = response[i] + 1
}
range(response)
lambda = 10^seq(from=-3,to=5,length.out=500)

#lambda = runif(10,min=650,max=1000)

output_CPP <- smooth.FEM(observations = response,
                         covariates = NULL,
                         FEMbasis = FEMbasis,
                         incidence_matrix = incidence_matrix,
                         lambda = lambda,
                         lambda.selection.criterion = "grid",
                         lambda.selection.lossfunction = "GCV",
                         DOF.evaluation = "exact",
                         family="poisson")

plot(log10(lambda), output_CPP$optimization$GCV_vector, xlab="log10(lambda)", ylab="GCV")

lambda_opt <- output_CPP$optimization$lambda_position
plot(FEM(output_CPP$fit.FEM$coeff[,lambda_opt],FEMbasis = FEMbasis))

Mass = fdaPDE::CPP_get.FEM.Mass.Matrix(FEMbasis)

I_vec = Mass%*%output_CPP$f$coeff[,lambda_opt]

I_region = matrix(0, nrow= chicago$domain$lines$n, ncol=1)

for(i in 1:chicago$domain$lines$n){
   I_region[i] = sum( I_vec[incidence_matrix[i,]])
}

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

fitted = inv.link(I_region)
