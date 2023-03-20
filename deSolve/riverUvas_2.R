rm(list=ls())

library(femR)
library(Matrix)
source("utils.R")

# model grid
dx = 1
riverlen = 669
Nb = riverlen/dx # number of interval
                 # Nb + 1 mesh nodes
xgrid = seq(from=0, to=riverlen, length.out=(Nb+1)) # mesh nodes
edges = cbind( seq(from=1,to=Nb,by=1), seq(from=2,to=Nb+1,by=1))

neigh = neigh1D(edges)
boundary = matrix(0,nrow=(Nb+1), ncol=1)
boundary[1] = boundary[Nb+1] = 1

River <- list(nodes= as.matrix(xgrid), edges= edges, elements= edges, neigh= neigh, boundary= boundary)

# parameters 
area  = 0.52; area2 = 1.56; ased = 1.5
Vol   = area*dx
disp  = 0.40
D     = disp*area/dx
Q     = 0.0125
alpha = 4.5e-5
qlat  = 2.151e-6; Clat = 0.13
sorp  = 5.6e-5; sorpS = 1
kd    = 7e-5
CsB   = 0.13

K1 = (qlat/area + alpha + ased/area*sorp)
K2 = -(alpha*area/area2 + sorpS)
# Upstream concentration
Time = c(8.25, 8.399, 8.400, 11.399, 11.4, 50)*3600
C_up = c(0.13, 0.13, 1.73, 1.73, 0.13, 0.13)
Cfun = approxfun(Time, C_up, rule=2)

# model solution and plotting
Yini = c(rep(0.13, times=2*(Nb+1)), rep(0., times=Nb+1))
times = seq(from=8.25*3600, to= 24*3600, by=120)

PDE_parameters <- list("diffusion" = D, 
                       "transport" = as.matrix(Q/area),
                        "reaction" = qlat/area)
dt_ = 120

C_ = matrix(0, nrow=3*(Nb+1),ncol=length(times))
C_[,1] = Yini

PDE <- new(PDE_1D_isotropic_ORDER_1, River)
PDE$set_PDEparameters(PDE_parameters)
quadrature_nodes = PDE$get_quadrature_nodes()

dirichletBC = matrix(0, nrow=Nb+1,ncol=1) 

forcing_ = function(x){
    return(qlat*Clat/area* x^0)
}

#### Setting forcing term ####
force_ = forcing_(quadrature_nodes)
PDE$set_forcingTerm( as.matrix(force_) ) 

# Assembling matrices
PDE$init()
Stiff = PDE$Stiffness() 
Mass = PDE$Mass()
neumannBC = PDE$Stiffness()[Nb+1,]
#################################################

for(k in 2:length(times)){

dirichletBC[1]    = Cfun(times[k])
dirichletBC[Nb+1] = C_[(Nb+1),k-1]    #C_[(Nb+1),k-1]    # default deSolve (?)

PDE$set_dirichletBC(dirichletBC)      #C_[2*(Nb+1)+1, 3*(Nb+1), (k-1)]

#### setting Dirichlet BC on the first node ####
#Stiff[1,] = 0. * Stiff[1,]
#Stiff[1,1] = 1.0
#force_[1] = Cfun(times[k])

#Stiff[Nb+1,] = 0. * Stiff[Nb+1,]
#Stiff[Nb+1,Nb+1] = 1.0
#force_[Nb+1] = Cfun(times[k])
################################################

PDE$set_forcingTerm( as.matrix(force_)) 
res = PDE$solve() # matrices
Stiff = res$Stiff
Mass  = res$Mass 

# Implicit Euler
#rhs_ = Mass%*%C_[1:(Nb+1),k-1]
#rhs_ = 1/dt_ * rhs_
#rhs_ = rhs_ + PDE$force() 

#tmp1 = Mass %*% C_[(Nb+2):(2*(Nb+1)),(k-1)]
#tmp1 = alpha* tmp1

#tmp2 = Mass %*% C_[(2*(Nb+1)+1): (3*(Nb+1)), (k-1)]
#tmp2 = -ased/area * sorp * tmp2

#rhs_ = rhs_ + tmp1 + tmp2 
#rhs_ = as.vector(rhs_)

#A_ = (1/dt_ * Mass + Stiff)
 
#C_[1:(Nb+1),k] = as.vector(solve(A_, rhs_))

#C_[(Nb+2):(2*(Nb+1)),k] =       (1.+ dt_*K2)*C_[(Nb+2):(2*(Nb+1)),(k-1)] + 
#                                dt_*area/area2*alpha * C_[1:(Nb+1),k]    +
#                                dt_*sorpS*CsB 

#C_[(2*(Nb+1)+1):(3*(Nb+1)), k] = (1. - dt_*sorp)*C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] + 
#                               dt_*sorp*C_[1:(Nb+1),k]


# PDE - Implicit Euler
exchange = alpha*(C_[(Nb+2):(2*(Nb+1)),(k-1)] - C_[1:(Nb+1),(k-1)] )
sorption = sorp*(C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] - C_[1:(Nb+1),(k-1)])
forcingVec = PDE$force()

rhs_ = forcingVec + 1/dt_*(Mass%*%C_[1:(Nb+1),(k-1)] ) + Mass%*%exchange + ased/area* (Mass%*%sorption) 
A_ = (1/dt_ * Mass + Stiff)

C_[1:(Nb+1),k] = as.vector(solve(A_, rhs_))
C_[(Nb+2):(2*(Nb+1)),k] = C_[(Nb+2):(2*(Nb+1)),(k-1)] - dt_* area/area2*exchange 
                                                      + dt_* sorpS* (CsB - C_[(Nb+2):(2*(Nb+1)),(k-1)])  

C_[(2*(Nb+1)+1):(3*(Nb+1)), k] = C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] - dt_*sorption

}

x11()
filled.contour(x=times/3600, y=xgrid, t(C_[1:(Nb+1),]))

x11()
filled.contour(x=times/3600, y=xgrid, t(C_[(Nb+2):(2*(Nb+1)),]))

x11()
filled.contour(x=times/3600, y=xgrid, t(C_[(2*(Nb+1)+1):(3*(Nb+1)),]))

pdf("RiverUvas_femR.pdf")
filled.contour(x=times/3600, y=xgrid, t(C_[1:(Nb+1),]))
filled.contour(x=times/3600, y=xgrid, t(C_[(Nb+2):(2*(Nb+1)),]))
filled.contour(x=times/3600, y=xgrid, t(C_[(2*(Nb+1)+1):(3*(Nb+1)),]))
dev.off()




for(k in 2:length(times)){

dirichletBC[1]    = Cfun(times[k])
dirichletBC[Nb+1] = C_[(Nb+1),k-1]    #C_[(Nb+1),k-1]    # default deSolve (?)

PDE$set_dirichletBC(dirichletBC)      #C_[2*(Nb+1)+1, 3*(Nb+1), (k-1)]

PDE$set_forcingTerm( as.matrix(force_)) 
res = PDE$solve() # matrices
Stiff = res$Stiff
Mass  = res$Mass 

# PDE - Implicit Euler
exchange = alpha*(C_[(Nb+2):(2*(Nb+1)),(k-1)] - C_[1:(Nb+1),(k-1)] )
sorption = sorp*(C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] - C_[1:(Nb+1),(k-1)])
forcingVec = PDE$force()

rhs_ = forcingVec + 1/dt_*(Mass%*%C_[1:(Nb+1),(k-1)] ) + Mass%*%exchange + ased/area* (Mass%*%sorption) 
A_ = (1/dt_ * Mass + Stiff)

C_[1:(Nb+1),k] = as.vector(solve(A_, rhs_))

exchange = alpha*(C_[(Nb+2):(2*(Nb+1)),(k-1)] - C_[1:(Nb+1),k] )
sorption = sorp*(C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] - C_[1:(Nb+1),k])

C_[(Nb+2):(2*(Nb+1)),k] = C_[(Nb+2):(2*(Nb+1)),(k-1)] - dt_* area/area2*exchange 
                                                      + dt_ * sorpS * (CsB - C_[(Nb+2):(2*(Nb+1)),k])  

C_[(2*(Nb+1)+1):(3*(Nb+1)), k] = C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] - dt_*sorption

}

pdf("RiverUvas_backwardEuler_ODE.pdf")
filled.contour(x=times/3600, y=xgrid, t(C_[1:(Nb+1),]))
filled.contour(x=times/3600, y=xgrid, t(C_[(Nb+2):(2*(Nb+1)),]))
filled.contour(x=times/3600, y=xgrid, t(C_[(2*(Nb+1)+1):(3*(Nb+1)),]))
dev.off()


{

library(fdaPDE)
mesh_aux = create.mesh.1.5D(cbind(River$nodes, rep(0,nrow=(Nb+1))), River$elements)

FEMbasis = create.FEM.basis(mesh_aux)

PDE_parameters <- list("diffusion" = D, 
                       "transport" = as.matrix(Q/area),
                        "reaction" = (qlat/area+alpha+ ased/area*sorp))
dt_ = 120

C_ = matrix(0, nrow=3*(Nb+1),ncol=length(times))
C_[,1] = Yini

PDE <- new(PDE_1D_isotropic_ORDER_1, River)
PDE$set_PDEparameters(PDE_parameters)
quadrature_nodes = PDE$get_quadrature_nodes()

dirichletBC = matrix(0, nrow=Nb+1,ncol=1) 

forcing_ = function(x){
    return(qlat*Clat/area* x^0)
}

#### Setting forcing term ####
force_ = forcing_(quadrature_nodes)
PDE$set_forcingTerm( as.matrix(force_) ) 

for(k in 2:length(times)){

dirichletBC[1]    = Cfun(times[k])
dirichletBC[Nb+1] = C_[(Nb+1),k-1]    #C_[(Nb+1),k-1]    # default deSolve (?)

PDE$set_dirichletBC(dirichletBC)      #C_[2*(Nb+1)+1, 3*(Nb+1), (k-1)]

forcing_TMP = as.matrix(force_) 
              + eval.FEM( FEM(coeff= alpha          *C_[(Nb+2):(2*(Nb+1)),(k-1)],       FEMbasis=FEMbasis), locations=quadrature_nodes )
              + eval.FEM( FEM(coeff= ased/area*sorp *C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)], FEMbasis=FEMbasis), locations=quadrature_nodes )
PDE$set_forcingTerm(forcing_TMP) 
res = PDE$solve() # matrices
Stiff = res$Stiff
Mass  = res$Mass 

# PDE - Implicit Euler
forcingVec = PDE$force() 

rhs_ = forcingVec + 1/dt_*(Mass%*%C_[1:(Nb+1),(k-1)] ) #+ alpha*(Mass%*%C_[(Nb+2):(2*(Nb+1)),(k-1)]) 
                                                       #+ ased/area*sorp * (Mass%*%C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)]) 
A_ = (1/dt_ * Mass + Stiff)

C_[1:(Nb+1),k] = as.vector(solve(A_, rhs_))

#exchange = alpha*(C_[(Nb+2):(2*(Nb+1)),(k-1)] - C_[1:(Nb+1),(k-1)] )
#sorption = sorp*(C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] - C_[1:(Nb+1),(k-1)])
C_[(Nb+2):(2*(Nb+1)),k] = (         1/dt_ * C_[(Nb+2):(2*(Nb+1)),(k-1)] 
                         + alpha*area/area2*C_[(Nb+2):(2*(Nb+1)),k] 
                         +           sorpS* CsB) / (1/dt_ + alpha*area/area2 + sorpS)  

C_[(2*(Nb+1)+1):(3*(Nb+1)), k] = (1/dt_*C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] + sorp*C_[1:(Nb+1),(k-1)] ) / (1/dt_ + sorp)

}

x11()
filled.contour(x=times/3600, y=xgrid, t(C_[1:(Nb+1),]))

x11()
filled.contour(x=times/3600, y=xgrid, t(C_[(Nb+2):(2*(Nb+1)),]))

x11()
filled.contour(x=times/3600, y=xgrid, t(C_[(2*(Nb+1)+1):(3*(Nb+1)),]))

pdf("RiverUvas_femR_mody.pdf")
image(x=times/3600, y=xgrid, t(C_[1:(Nb+1),]))
image(x=times/3600, y=xgrid, t(C_[(Nb+2):(2*(Nb+1)),]))
image(x=times/3600, y=xgrid, t(C_[(2*(Nb+1)+1):(3*(Nb+1)),]))

x11()
plot(times/3600, C_[38,], type="l", xlab="time, hour", ylab="Conc", main="Strontium",ylim=range(C_))
points(times/3600, C_[38,])

points(times/3600, C_[281,], type="l", xlab="time, hour", ylab="Conc", main="Strontium")
points(times/3600, C_[281,], pch=16)
legend("topright", legend=c("38m", "281m"), pch=c(1,16))




}
