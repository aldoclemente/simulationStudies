### TEST AREALE ###
source("refine1D.R")

eps = 1/2
x = c(0., 1)
y = c(0.,eps)
vertices = expand.grid(x,y)
vertices = cbind(vertices[,1], vertices[,2])
edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)

delta = 0.0125
mesh = create.mesh.1.5D(vertices, edges)

# areale
nregion=3
I_region = matrix(0, nrow=nregion, ncol=1)

new_to_old  = refine1D(mesh$nodes, mesh$edges, delta)$new_to_old

incidence_matrix = matrix(0, nrow=nrow(I_region), ncol=nrow(new_to_old))
for(i in 1:nrow(new_to_old)){
  incidence_matrix[ new_to_old[i],i] = 1
}


mesh = refine.mesh.1.5D(mesh,delta=delta)

nnodes=dim(mesh$nodes)[1]
FEMbasis=create.FEM.basis(mesh)
plot(mesh, asp=1)

# Exact solution (pointwise at nodes)
aux.4 = function(x,y){
  h = 1
  source = 4 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta < h ){
      coef[i] = 3*exp(1/(delta^2 / h^2 -1)+1 ) - 1
      
    }
  }
  
  return(coef)
}
aux.3 = function(x,y){
  
  h = eps
  source = 3 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,1]-mesh$nodes[source,1]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,2] - mesh$nodes[source,2])
    if(delta < h ){
      coef[i] = -1 - 1/h*delta
    }
    
    
  }
  return(coef)
}
aux.1 = function(x,y){
  
  h = 1
  source = 1 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta <= h ){
      coef[i] = -2 - 1/h*delta
    }
  }
  return(coef)
  
}  
AUX = function(x,y){
  
  res = aux.1(x,y) + aux.3(x,y) + aux.4(x,y)
  return(res)
}

sol_exact=AUX(mesh$nodes[,1],mesh$nodes[,2])
plot(FEM(sol_exact, FEMbasis))

# Add error to simulate data
set.seed(32)
ran=range(sol_exact)
data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*abs(ran[2]-ran[1]))

# Set smoothing parameter
lambda = 10^seq(-4,-2,length.out=10)
Mass = fdaPDE::CPP_get.FEM.Mass.Matrix(FEMbasis)

output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       lambda.selection.criterion='grid', DOF.evaluation='exact', lambda.selection.lossfunction='GCV')

plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(FEM(output_CPP$fit.FEM$coeff,FEMbasis))
L2.error = sqrt( sum((Mass%*%(sol_exact - output_CPP$fit.FEM$coeff)^2 ) ) )

# test areale #

I_vec_exact = Mass %*% sol_exact

# ok
I1 = -2.5       # down 
I2 = -0.75  
I3 = 0.81035    # up
I_ex = I1 + I2 + I3
abs( I_ex - sum(I_vec_exact))
# 
I_region[1] = I1
I_region[2] = I2
I_region[3] = I3

lens = matrix(0, nrow=nrow(mesh$edges), ncol=1) 
for(i in 1:nrow(mesh$edges)){
  lens[i] = norm( mesh$nodes[mesh$edges[i,1],]-mesh$nodes[mesh$edges[i,2],], "2")
}

lens_region = matrix(0, nrow=nregion, ncol=1)
for(i in 1:nregion){
  lens_region[i] = sum( lens[which(incidence_matrix[i,]==1)] )
}

data = I_region/lens_region + rnorm(nregion, mean=0, sd = 0.05*abs(range(I_region)))
lambda = 10^seq(-4,-2,length.out=10)

output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                       incidence_matrix = incidence_matrix,
                       lambda.selection.criterion='grid', 
                       DOF.evaluation='exact',
                       lambda.selection.lossfunction='GCV')

plot(log10(lambda), output_CPP$optimization$GCV_vector)
plot(output_CPP$fit.FEM)
L2.error = sqrt( sum((Mass%*%(sol_exact - output_CPP$fit.FEM$coeff)^2 ) ) )

M = 30
rmse = matrix(0, nrow=M, ncol=1)
L2.error = matrix(0, nrow=M, ncol=1)
nnodes = nrow(mesh$nodes)
points_ = refine.mesh.1.5D(mesh, delta= delta/8)$nodes
npoints = nrow(points_)
for(i in 1:M){
  data = I_region/lens_region + rnorm(nregion, mean=0, sd = 0.05*abs(range(I_region)))
  output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                         incidence_matrix = incidence_matrix,
                         lambda.selection.criterion='grid', 
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV')
  rmse[i] = sqrt( 1/npoints * sum((eval.FEM(FEM(sol_exact,FEMbasis), points_) - eval.FEM(output_CPP$fit.FEM, points_))^2 ) )
  L2.error[i] = sqrt( sum((Mass%*%(sol_exact - output_CPP$fit.FEM$coeff)^2 ) ) )
  
}
boxplot(rmse, main="RMSE")
boxplot(L2.error, main="L2 error")
################################################################################
source("refine1D.R")
eps = 1/2
x = c(0., 1)
y = c(0.,eps)
vertices = expand.grid(x,y)
vertices = cbind(vertices[,1], vertices[,2])
edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)

delta = 0.0125/2
mesh = create.mesh.1.5D(vertices, edges)
mesh = refine.mesh.1.5D(mesh,delta=40*delta)

# areale
new_to_old  = refine1D(mesh$nodes, mesh$edges, delta)$new_to_old

nregion = nrow(mesh$edges)
idxs = list()
for(i in 1:nregion){
  idxs[[i]] = (i-1)*nrow(new_to_old)/nregion + seq(from=1,to=nrow(new_to_old)/nregion)
}

mesh = refine.mesh.1.5D(mesh,delta=delta)
nnodes=dim(mesh$nodes)[1]
FEMbasis=create.FEM.basis(mesh)
plot(mesh, asp=1)

I_region = matrix(0, nrow=nregion, ncol=1)

incidence_matrix = matrix(0, nrow=nrow(I_region), ncol=nrow(mesh$edges))
for(i in 1:nrow(new_to_old)){
  incidence_matrix[ new_to_old[i],i] = 1
}

# Exact solution (pointwise at nodes)
aux.4 = function(x,y){
  h = 1
  source = 4 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta < h ){
      coef[i] = 3*exp(1/(delta^2 / h^2 -1)+1 ) - 1
      
    }
  }
  
  return(coef)
}
aux.3 = function(x,y){
  
  h = eps
  source = 3 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,1]-mesh$nodes[source,1]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,2] - mesh$nodes[source,2])
    if(delta < h ){
      coef[i] = -1 - 1/h*delta
    }
    
    
  }
  return(coef)
}
aux.1 = function(x,y){
  
  h = 1
  source = 1 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta <= h ){
      coef[i] = -2 - 1/h*delta
    }
  }
  return(coef)
  
}  
AUX = function(x,y){
  
  res = aux.1(x,y) + aux.3(x,y) + aux.4(x,y)
  return(res)
}

sol_exact=AUX(mesh$nodes[,1],mesh$nodes[,2])
for(i in 1:nregion){
  for(e in idxs[[i]]){
    I_region[i] = I_region[i] + trapezoidal(c(sol_exact[mesh$edges[e,1]], sol_exact[mesh$edges[e,2]]),
                                              mesh$nodes[mesh$edges[e,1],],
                                              mesh$nodes[mesh$edges[e,2],]) 
  }
}
sum(I_region)

lens = matrix(0, nrow=nrow(mesh$edges), ncol=1) 
for(i in 1:nrow(mesh$edges)){
  lens[i] = norm( mesh$nodes[mesh$edges[i,1],]-mesh$nodes[mesh$edges[i,2],], "2")
}

lens_region = matrix(0, nrow=nregion, ncol=1)
for(i in 1:nregion){
  lens_region[i] = sum( lens[which(incidence_matrix[i,]==1)] )
}


Mass = fdaPDE::CPP_get.FEM.Mass.Matrix(FEMbasis)
I_ex = sum(Mass %*% sol_exact)
abs(I_ex - sum(I_region))

M = 30
rmse = matrix(0, nrow=M, ncol=1)
L2.error = matrix(0, nrow=M, ncol=1)
nnodes = nrow(mesh$nodes)
points_ = refine.mesh.1.5D(mesh, delta= delta/8)$nodes
npoints = nrow(points_)
for(i in 1:M){
  data = I_region/lens_region + rnorm(nregion, mean=0, sd = 0.05*abs(range(I_region)))
  
  param 
  output_CPP<-smooth.FEM(observations=data, FEMbasis=FEMbasis, lambda=lambda,
                         incidence_matrix = incidence_matrix,
                         lambda.selection.criterion='grid', 
                         DOF.evaluation='exact',
                         lambda.selection.lossfunction='GCV')
  rmse[i] = sqrt( 1/npoints * sum((eval.FEM(FEM(sol_exact,FEMbasis), points_) - 
                                   eval.FEM(output_CPP$fit.FEM$coeff[,lambda_opt], points_))^2 ) )
  L2.error[i] = sqrt( sum((Mass%*%(sol_exact - output_CPP$fit.FEM$coeff$coeff[,lambda_opt])^2 ) ) )
  
}
boxplot(rmse, main="RMSE")
boxplot(L2.error, main="L2 error")




################################################################################
### TEST AREALE GLM ###
# family
setwd("SR-chicago-case-study/")
library(fdaPDE)
library(purrr)
rm(list=ls())
graphics.off()
source("refine1D.R")

#family
FAMILY = "poisson"
l<-make.link("log")
link<-l$linkfun
inv.link<-l$linkinv

# FAMILY = "binomial"
# logit <- function(x){qlogis(x)}
# inv.logit <- function(x){plogis(x)}
# link = logit
# inv.link = inv.logit

# scale param
scale.param = 1

# mesh
eps = 1 / 2
x = c(0., 1)
y = c(0.,eps)
vertices = expand.grid(x,y)
vertices = cbind(vertices[,1], vertices[,2])
edges = matrix(c(1,2,1,3,3,4), nrow=3,ncol=2, byrow=T)

delta = 0.0125/2
mesh = create.mesh.1.5D(vertices, edges)
mesh = refine.mesh.1.5D(mesh,delta=40*delta)

# areale
new_to_old  = refine1D(mesh$nodes, mesh$edges, delta)$new_to_old

nregion = nrow(mesh$edges)
idxs = list()
for(i in 1:nregion){
  idxs[[i]] = (i-1)*nrow(new_to_old)/nregion + seq(from=1,to=nrow(new_to_old)/nregion)
}

mesh = refine.mesh.1.5D(mesh,delta=delta)
nnodes=dim(mesh$nodes)[1]
FEMbasis=create.FEM.basis(mesh)
plot(mesh, asp=1)

incidence_matrix = matrix(0, nrow=nregion, ncol=nrow(mesh$edges))
for(i in 1:nrow(new_to_old)){
  incidence_matrix[ new_to_old[i],i] = 1
}

# 1.5D spatial field (function f)
aux.4 = function(x,y){
  h = 1
  source = 4 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta < h ){
      coef[i] = (3*exp(1/(delta^2 / h^2 -1)+1 ) - 1)
      
    }
  }
  
  return(coef)
}
aux.3 = function(x,y){
  
  h = eps
  source = 3 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,1]-mesh$nodes[source,1]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,2] - mesh$nodes[source,2])
    if(delta < h ){
      coef[i] = -1 - 1/h*delta
    }
    
    
  }
  return(coef)
}
aux.1 = function(x,y){
  
  h = 1
  source = 1 
  points_ = cbind(x,y)
  idx.ok = which( abs(points_[,2]-mesh$nodes[source,2]) < 10 * .Machine$double.eps )
  coef = vector(mode="numeric", length=length(x))
  
  for(i in idx.ok){
    delta = abs(points_[i,1] - mesh$nodes[source,1])
    if(delta <= h ){
      coef[i] = -2 - 1/h*delta
    }
  }
  return(coef)
  
}  
f = function(x,y){
  
  res = (aux.1(x,y) + aux.3(x,y) + aux.4(x,y)) + 1.
  return(res)
}

sol_exact=f(mesh$nodes[,1],mesh$nodes[,2])
I_region = matrix(0, nrow= nregion, ncol= 1) 
for(i in 1:nregion){
  for(e in idxs[[i]]){
    I_region[i] = I_region[i] + trapezoidal(c(sol_exact[mesh$edges[e,1]], sol_exact[mesh$edges[e,2]]),
                                            mesh$nodes[mesh$edges[e,1],],
                                            mesh$nodes[mesh$edges[e,2],]) 
  }
}
sum(I_region)

lens = matrix(0, nrow=nrow(mesh$edges), ncol=1) 
for(i in 1:nrow(mesh$edges)){
  lens[i] = norm( mesh$nodes[mesh$edges[i,1],]-mesh$nodes[mesh$edges[i,2],], "2")
}

lens_region = matrix(0, nrow=nregion, ncol=1)
for(i in 1:nregion){
  lens_region[i] = sum( lens[which(incidence_matrix[i,]==1)] )
}

# samoling response
#X = matrix(rbeta(nregion,1,2), nrow=nregion, ncol=1)
#beta_ = 1.
param= I_region  #/lens_region # + desmat%*%betas_truth 
ran=range(param) 

mu<-inv.link(param)
range(mu)
response <- rpois(nregion, lambda = mu)
#response <- rbernoulli(nregion, mu)
range(response)

#lambda = 10^seq(-4.,0,length.out = 250)

lambda = 10^seq(-2.,2,length.out = 500)
output_CPP<- smooth.FEM(observations = as.numeric(response), FEMbasis =FEMbasis, 
                        covariates = NULL,
                        incidence_matrix =incidence_matrix,
                        max.steps=30, family=FAMILY, mu0=NULL, scale.param=NULL,
                        lambda = lambda, lambda.selection.criterion = 'grid', DOF.evaluation = 'exact', lambda.selection.lossfunction = 'GCV')
plot(log10(lambda), output_CPP$optimization$GCV_vector)
lambda_opt <- output_CPP$optimization$lambda_position
plot(FEM(output_CPP$fit.FEM$coeff[,lambda_opt], FEMbasis))
plot(FEM(sol_exact, FEMbasis))
range(output_CPP$fit.FEM$coeff[,lambda_opt])

locations = refine.by.splitting.mesh.1.5D(mesh)$nodes
prediction = inv.link( eval.FEM( FEM(output_CPP$fit.FEM$coeff[,lambda_opt], FEMbasis), locations = locations))
cbind(prediction, inv.link(f(locations[,1],locations[,2])) )

Mass = CPP_get.FEM.Mass.Matrix(FEMbasis = FEMbasis)
M = 30
rmse = matrix(0, nrow=M, ncol=1)
rmse.2 = matrix(0,nrow=M, ncol=1)
L2.error = matrix(0, nrow=M, ncol=1)
rmse.beta = matrix(0, nrow=M, ncol=1)
nnodes = nrow(mesh$nodes)
points_ = refine.mesh.1.5D(mesh, delta= delta/8)$nodes
npoints = nrow(points_)
set.seed(NULL)
set.seed(314156)
for(i in 1:M){
  
  response <- rpois(nregion, lambda = mu)
  
  output_CPP<- smooth.FEM(observations = as.numeric(response), 
                          FEMbasis =FEMbasis,
                          covariates = X,
                          incidence_matrix =incidence_matrix,
                          max.steps=30, family=FAMILY,
                          lambda = lambda, lambda.selection.criterion = 'grid', 
                          DOF.evaluation = 'exact', 
                          .selection.lossfunction = 'GCV')
  lambda_opt = output_CPP$optimization$lambda_position
  rmse[i] = sqrt( mean((eval.FEM(FEM(sol_exact,FEMbasis), points_) - 
                        eval.FEM(FEM(output_CPP$fit.FEM$coeff[,lambda_opt], FEMbasis), points_))^2))
  
  rmse.2[i] = sqrt(mean( ))
  rmse.beta[i] = 1/M * (output_CPP$solution$beta-beta_)^2
  L2.error[i] = sqrt( sum((Mass%*%(sol_exact - output_CPP$fit.FEM$coeff)^2 ) ) )
  
}
boxplot(rmse, main="f")
boxplot(rmse.beta, main="Beta")
boxplot(L2.error, main="L2 error")

