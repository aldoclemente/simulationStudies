library(fdaPDE)
library(spatstat)
setwd("SR-chicago-case-study/")

source("refine1D.R")
data("chicago")

vertices = cbind(chicago$domain$vertices$x, chicago$domain$vertices$y)
edges = cbind( chicago$domain$from, chicago$domain$to)

mesh = create.mesh.1.5D(nodes = vertices, edges = edges)

x11()
plot(mesh, pch=".")
points(chicago$data$x, chicago$data$y, col = "red", pch=16, cex=2.5)

### FAMILY - poisson ###
delta = 65
new_to_old  = refine1D(mesh$nodes, mesh$edges, delta)$new_to_old
mesh = refine.mesh.1.5D(mesh, delta=delta)
FEMbasis = create.FEM.basis(mesh)

#setting regions 
nregion = 10
# lines == true edges of the network / edges == edges of the discretized network
# lines_to_region <- set_region(sources_ = c(47, 91, 67, 124, 151, 169, 195, 251, 275, 302))

incidence_matrix = matrix(0, nrow=nregion, ncol=nrow(mesh$edges))

for(i in 1:nrow(new_to_old)){
  incidence_matrix[ lines_to_region[new_to_old[i]],i] = 1
}

response = rep(0, times= nregion)

for( i in chicago$data$seg){
  response[lines_to_region[i]] = response[lines_to_region[i]] + 1
}
range(response)

x11()
plot_region(lines_to_region, response)

lambda = 10^seq(from=-5,to=-4,length.out=10)

output_CPP <- smooth.FEM(observations = response,
                         covariates = NULL,
                         FEMbasis = FEMbasis,
                         incidence_matrix = incidence_matrix,
                         lambda = lambda,
                         lambda.selection.criterion = "grid",
                         lambda.selection.lossfunction = "GCV",
                         DOF.evaluation = "exact",
                         family="poisson")


lambda_opt <- output_CPP$optimization$lambda_position

x11()
plot(log10(lambda), output_CPP$optimization$GCV_vector, xlab="log10(lambda)", ylab="GCV")

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
