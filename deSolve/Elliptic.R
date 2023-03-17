rm(list=ls())
library(femR)
source('utils.R')

dx = 0.0125
riverlen = 1
Nb = riverlen/dx # number of interval
                 # Nb + 1 mesh nodes

xgrid = seq(from=0, to=riverlen, length.out=(Nb+1)) # mesh nodes
edges = cbind( seq(from=1,to=Nb,by=1), seq(from=2,to=Nb+1,by=1))

neigh = neigh1D(edges)
boundary = matrix(0,nrow=(Nb+1), ncol=1)
boundary[1] = boundary[Nb+1] = 1

River <- list(nodes= as.matrix(xgrid), edges= edges, elements= edges, neigh= neigh, boundary= boundary)
PDE_parameters = list("diffusion" = 1., 
                       "transport" = as.matrix(0.),
                        "reaction" = 1.)

U_t_0 = 0. 
U_t_1 = 0.

dt = 0.025
T = 1.5
times = seq(from=0, to=T, by=dt)

U_ = matrix(0, nrow=Nb+1, ncol=length(times))
U_[,1] = sin(pi*xgrid)

forcing_ = function(x,t){
    return((pi^2-1.)*exp(-t)*sin(pi*x))}

dirichletBC       = matrix(0, nrow=Nb+1,ncol=1) 
dirichletBC[1]    = U_t_0
dirichletBC[Nb+1] = U_t_1

PDE <- new(PDE_1D_isotropic_ORDER_1, River)
PDE$set_PDEparameters(PDE_parameters)
quadrature_nodes = PDE$get_quadrature_nodes()

for(k in 2:length(times)){

force_ = forcing_(quadrature_nodes,rep(times[k], times=nrow(quadrature_nodes)))
PDE$set_forcingTerm( as.matrix(force_) ) 

PDE$set_dirichletBC(dirichletBC)

res = PDE$solve() # matrices
Stiff = res$Stiff
Mass  = res$Mass 

# Implicit Euler
rhs_ = Mass%*%U_[,k-1]
rhs_ = 1/dt * rhs_
rhs_ = rhs_ + PDE$force()
rhs_ = as.vector(rhs_)
A_ = (1/dt * Mass + Stiff)
 
U_[,k] = solve(A_, rhs_)
}

U_ex = function(x,t){
    res = matrix(0, nrow=length(x), ncol=length(t))
    for(i in 1:length(x))
        for(j in 1:length(t))
            res[i,j] = exp(-t[j]) * sin(pi*x[i])

    return(res)        

}

Solu_ = U_ex(x=xgrid, t=times)

x11()
filled.contour(x=xgrid, y=times, U_)
x11()
filled.contour(x=xgrid, y=times, Solu_)


# Non-homogeneous Dirichlet BC
rm(list=ls())
source('utils.R')
dx = 0.0125
riverlen = 1
Nb = riverlen/dx # number of interval
                 # Nb + 1 mesh nodes

xgrid = seq(from=0, to=riverlen, length.out=(Nb+1)) # mesh nodes
edges = cbind( seq(from=1,to=Nb,by=1), seq(from=2,to=Nb+1,by=1))

neigh = neigh1D(edges)
boundary = matrix(0,nrow=(Nb+1), ncol=1)
boundary[1] = boundary[Nb+1] = 1

River <- list(nodes= as.matrix(xgrid), edges= edges, elements= edges, neigh= neigh, boundary= boundary)
PDE_parameters = list("diffusion" = 1., 
                       "transport" = as.matrix(0.),
                        "reaction" = 0.)

U_t_0 = 0. 
U_t_1 = function(t) {
    return(sin(1)*exp(-t))
}

dt = 0.025
T = 1.5
times = seq(from=0, to=T, by=dt)

U_ = matrix(0, nrow=Nb+1, ncol=length(times))
U_[,1] = sin(xgrid)

forcing_ = function(x,t){
    return(0.0 * x^0)}

for(k in 2:length(times)){

PDE <- new(PDE_1D_isotropic_ORDER_1, River)
PDE$set_PDEparameters(PDE_parameters)
quadrature_nodes = PDE$get_quadrature_nodes()

force_ = forcing_(quadrature_nodes,rep(times[k], times=nrow(quadrature_nodes)))
PDE$set_forcingTerm( as.matrix(force_) ) 

dirichletBC       = matrix(0, nrow=Nb+1,ncol=1) 
dirichletBC[1]    = U_t_0
dirichletBC[Nb+1] = U_t_1(times[k])
PDE$set_dirichletBC(dirichletBC)

res = PDE$solve() # matrices
Stiff = res$Stiff
Mass  = res$Mass 

# Implicit Euler
rhs_ = Mass%*%U_[,k-1]
rhs_ = 1/dt * rhs_
rhs_ = rhs_ + PDE$force()
rhs_ = as.vector(rhs_)
A_ = (1/dt * Mass + Stiff)
 
U_[,k] = solve(A_, rhs_)
}

x11()
filled.contour(x=xgrid, y=times, U_)

U_ex = function(x,t){
    res = matrix(0, nrow=length(x), ncol=length(t))
    for(i in 1:length(x))
        for(j in 1:length(t))
            res[i,j] = exp(-t[j]) * sin(x[i])

    return(res)        

}

Solu_ = U_ex(x=xgrid, t=times)

x11()
filled.contour(x=xgrid, y=times, U_)
x11()
filled.contour(x=xgrid, y=times, Solu_)
