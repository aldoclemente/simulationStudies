rm(list=ls())

library(femR)
library(Matrix)
source("utils.R")

# model grid
dx = 0.0125
riverlen = 1
Nb = riverlen/dx # number of interval
                 # Nb + 1 mesh nodes
xgrid = seq(from=0, to=riverlen, length.out=(Nb+1)) # mesh nodes
edges = cbind( seq(from=1,to=Nb,by=1), seq(from=2,to=Nb+1,by=1))
edges

neigh = neigh1D(edges)
boundary = matrix(0,nrow=(Nb+1), ncol=1)
boundary[1] = boundary[Nb+1] = 1

River <- list(nodes= as.matrix(xgrid), edges= edges, elements= edges, neigh= neigh, boundary= boundary)

# parameters 
mu = 1.
beta     = 5.
Q = 0.25

# model solution and plotting
Yini = sin(pi*xgrid)

times = seq(from=0, to= 1., by=0.025)

PDE_parameters <- list("diffusion" = mu, 
                       "transport" = as.matrix(beta),
                        "reaction" = Q)
dt_ = 0.125

C_ = matrix(0, nrow=2*(Nb+1),ncol=length(times))
C_[1:(Nb+1),1] = Yini

PDE <- new(PDE_1D_isotropic_ORDER_1, River)
PDE$set_PDEparameters(PDE_parameters)
quadrature_nodes = PDE$get_quadrature_nodes()

dirichletBC       = matrix(0, nrow=Nb+1,ncol=1) 
dirichletBC[1]    = 0.
dirichletBC[Nb+1] = 0.

for(k in 2:length(times)){

PDE$set_forcingTerm(  matrix(0., nrow= nrow(quadrature_nodes), ncol=1))
PDE$set_dirichletBC(dirichletBC)

res = PDE$solve() # matrices
Stiff = res$Stiff
Mass  = res$Mass 

# Implicit Euler
rhs_ = Mass%*%C_[1:(Nb+1),k-1]
rhs_ = 1/dt_ * rhs_

tmp1 = Mass %*% C_[(Nb+2):(2*(Nb+1)),(k-1)]
tmp1 = Q * tmp1

rhs_ = rhs_ + tmp1

A_ = (1/dt_ * Mass + Stiff)
 
C_[1:(Nb+1),k] = as.vector(solve(A_, rhs_))

C_[(Nb+2):(2*(Nb+1)),k] =       (1.+ Q*dt_)*C_[(Nb+2):(2*(Nb+1)),(k-1)] + 
                                dt_*Q   * C_[1:(Nb+1),k]    

}


x11()
filled.contour(x=xgrid, y=times, C_[1:(Nb+1),], main="PDE")

x11()
filled.contour(x=xgrid, y=times, C_[(Nb+2):(2*(Nb+1)),], main="ODE")

x11()
filled.contour(x=xgrid, y=times, (C_[1:(Nb+1),] + C_[(Nb+2):(2*(Nb+1)),]), main="Somma")