library(fdaPDE)
library(spatstat)
setwd("SR-chicago-case-study/")

source("refine1D.R")
source("../utils.R")
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

x11()
plot(mesh, pch=".")
points(DATA$data$x, DATA$data$y, col = "red", pch=16, cex=1.5)

### FAMILY - poisson ###
delta = 150

mesh = refine.mesh.1.5D(mesh, delta=delta)
FEMbasis = create.FEM.basis(mesh)
new_to_old  = refine1D(mesh$nodes, mesh$edges, delta)$new_to_old

#setting regions 
nregion = 9 #8 #6 #10
ndata = length(DATA$data$x)
# lines == true edges of the network / edges == edges of the discretized network
data_ = cbind(DATA$data$x, DATA$data$y)
result_ <- kmeans(x=data_, centers=nregion, iter.max = 100)
centroids_ = projection.points.1.5D(mesh, locations= result_$centers)

x11()
plot(mesh, pch=".")
points(data_, col="red", pch=16, cex=2)
points(result_$centers, col="blue", pch = 16, cex=2)
points(centroids_, col="green4", pch=16, cex=2)
legend("topright", legend=c("data", "2D", "1.5D"), 
        col=c("red", "blue", "green4"), pch = 16, cex=1.25 )

lines_to_region <- set_region(centroids_, mesh=mesh, LN = DATA)

incidence_matrix = matrix(0, nrow=nregion, ncol=nrow(mesh$edges))

for(i in 1:nrow(new_to_old)){
  incidence_matrix[ lines_to_region[new_to_old[i]],i] = 1
}

response = rep(0, times= nregion)

for( i in DATA$data$seg){
  response[lines_to_region[i]] = response[lines_to_region[i]] + 1
}
range(response)

x11()
plot_region(lines_to_region, response, LN=DATA,mesh=mesh)

lambda = 10^seq(from=0,to=5,length.out=150)

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
# sqrt(GCV_ * (sum(!is.na(observations)) - dof)/sum(!is.na(observations)))
x11()
plot(log10(lambda), output_CPP$optimization$GCV_vector, xlab="log10(lambda)", ylab="GCV")

density.FEM = FEM(coeff= output_CPP$fit.FEM$coeff[,lambda_opt] / ndata, FEMbasis= FEMbasis )

plot(FEM(output_CPP$fit.FEM$coeff[,lambda_opt],FEMbasis = FEMbasis))
plot(density.FEM) 

Mass = fdaPDE::CPP_get.FEM.Mass.Matrix(FEMbasis)

I_vec = Mass%*%density.FEM$coeff
sum(I_vec)

I_region = matrix(0, nrow= nregion, ncol=1)

for(i in 1:nregion){
   I_region[i] = sum( I_vec[incidence_matrix[i,]])
}

l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

fitted = inv.link(I_region)
