rm(list=ls())

library(femR)
library(Matrix)
source("utils.R")

# model grid
dx = 1
riverlen = 669
Nb = riverlen/dx # number of interval
                 # Nb + 1 mesh nodes
nodes = seq(from=0, to=riverlen, length.out=(Nb+1)) # mesh nodes
edges = cbind( seq(from=1,to=Nb,by=1), seq(from=2,to=Nb+1,by=1))

neigh = neigh1D(edges)
boundary = matrix(0,nrow=(Nb+1), ncol=1)
boundary[1] = boundary[Nb+1] = 1

River <- list(nodes= as.matrix(nodes), edges= edges, elements= edges, neigh= neigh, boundary= boundary)

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

# Upstream concentration
Time = c(8.25, 8.399, 8.400, 11.399, 11.4, 50)*3600
C_up = c(0.13, 0.13, 1.73, 1.73, 0.13, 0.13)
Cfun = approxfun(Time, C_up, rule=2)

# model solution and plotting
Yini = c(rep(0.13, times=2*(Nb+1)), rep(0., times=Nb+1))
times = seq(from=8.25*3600, to= 24*3600, by=120) 

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
dirichletBC[Nb+1] = C_[(Nb+1),k-1]    

PDE$set_dirichletBC(dirichletBC)      

PDE$set_forcingTerm(force_) 
res = PDE$solve() # matrices
Stiff = res$Stiff
Mass  = res$Mass 

# PDE - Implicit Euler
forcingVec = PDE$force() 

rhs_ = forcingVec + 1/dt_*(Mass%*%C_[1:(Nb+1),(k-1)] ) + alpha*(Mass%*%C_[(Nb+2):(2*(Nb+1)),(k-1)]) + ased/area*sorp * (Mass%*%C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)]) 
A_ = (1/dt_ * Mass + Stiff)

C_[1:(Nb+1),k] = as.vector(solve(A_, rhs_))

exchange = alpha*(C_[(Nb+2):(2*(Nb+1)),(k-1)] - C_[1:(Nb+1),k] )
sorption = sorp*(C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] - C_[1:(Nb+1),k])

C_[(Nb+2):(2*(Nb+1)),k] = (1+alpha*area/ased*dt_+sorpS*dt_)^(-1) * (C_[(Nb+2):(2*(Nb+1)),(k-1)] + alpha*area/ased*dt_* C_[1:(Nb+1),k] + sorpS*dt_*CsB)
C_[(2*(Nb+1)+1):(3*(Nb+1)), k] = (1+sorp*dt_)^(-1) * (C_[(2*(Nb+1)+1):(3*(Nb+1)), (k-1)] + dt_*sorp*C_[1:(Nb+1), k])
}

x11()
filled.contour(x=times/3600, y=nodes, t(C_[1:(Nb+1),]))

x11()
filled.contour(x=times/3600, y=nodes, t(C_[(Nb+2):(2*(Nb+1)),]))

x11()
filled.contour(x=times/3600, y=nodes, t(C_[(2*(Nb+1)+1):(3*(Nb+1)),]))

pdf("RiverUvas_femR_Implicit_FullForward_Euler.pdf")
image(x=times/3600, y=nodes, t(C_[1:(Nb+1),]))
image(x=times/3600, y=nodes, t(C_[(Nb+2):(2*(Nb+1)),]))
image(x=times/3600, y=nodes, t(C_[(2*(Nb+1)+1):(3*(Nb+1)),]))
legend("topright", legend=c("38m", "281m"), pch=1)
dev.off()

library(ReacTran)

RiverModel = function(time, state, pars){
    C    = state[1:Nb]
    Cs   = state[(Nb+1):(2*Nb)]
    Csed = state[(2*Nb+1):(3*Nb)]

    Cup = Cfun(time)

    tranC = tran.volume.1D(C=C, C.up = Cup, flow.lat= qlat, C.lat = Clat,
                            flow=Q, Disp = D, V = Vol)

    exchange = alpha * (Cs- C)
    sorption = sorp * (Csed-C)

    dC    = tranC$dC + exchange + sorption * ased/area
    dCs   = -exchange * area/area2 + sorp*(CsB-Cs)
    dCsed = -sorption
    list(c(dC, dCs, dCsed)) 
}

dx = 1
riverlen = 669
Nb = riverlen/dx 

xgrid = setup.grid.1D(N=Nb, L= riverlen)

# Upstream concentration
Time = c(8.25, 8.399, 8.400, 11.399, 11.4, 50)*3600
C_up = c(0.13, 0.13, 1.73, 1.73, 0.13, 0.13)
Cfun = approxfun(Time, C_up, rule=2)

# model solution and plotting
Yini = c(rep(0.13, times=2*Nb), rep(0., times=Nb))
times = seq(from=8.25*3600, to= 24*3600, by=120)

Dyna = ode.1D(y = Yini, func = RiverModel, nspec=3, times= times)

C_deSolve = Dyna[,2:(3*Nb+1)]


library(latex2exp)
x11()
par(mfrow=c(1,2))
image(x=times/3600, y=xgrid$x.mid, C_deSolve[,1:Nb], xlab="x", ylab="time")
mtext("deSolve", side = 3, line = 1)
image(x=times/3600, y=nodes, t(C_[1:(Nb+1),]), xlab="x", ylab="time")
mtext("femR", side = 3, line = 1)
mtext(TeX("$\\textbf{C}$"), side = 3, line = -2, outer = TRUE)

x11()
par(mfrow=c(1,2))
image(x=times/3600, y=xgrid$x.mid, C_deSolve[,(Nb+1):(2*Nb)], xlab="x", ylab="time")
mtext("deSolve", side = 3, line = 1)
image(x=times/3600, y=nodes, t(C_[(Nb+2):(2*(Nb+1)),]), xlab="x", ylab="time")
mtext("femR", side = 3, line = 1)
mtext(TeX("$\\textbf{C_{LS}}$"), side = 3, line = -2, outer = TRUE)

x11()
par(mfrow=c(1,2))
image(x=times/3600, y=xgrid$x.mid, C_deSolve[,(2*Nb+1):(3*Nb)],xlab="x", ylab="time")
mtext("deSolve", side = 3, line = 1)
image(x=times/3600, y=nodes, t(C_[(2*(Nb+1)+1):(3*(Nb+1)),]), xlab="x", ylab="time")
mtext("femR", side = 3, line = 1)
mtext(TeX("$\\textbf{C_{sed}}$"), side = 3, line = -2, outer = TRUE)


par(mfrow=c(1,1))
plot(times/3600, C_deSolve[,38], type="l", lwd=1.5, col="red3", xlab="time, hour", ylab="Conc",ylim=range(C_deSolve),
    main="Strontium (x = 38m)")
points(times/3600, C_[38,], type="l", lwd=1.5, col="blue2")
legend("topright", legend=c("deSolve", "femR"), 
        lty=c(1,1), lwd=3, col=c("red3", "blue2"))

par(mfrow=c(1,1))
plot(times/3600, C_deSolve[,281], type="l", lwd=1.5, col="red3", xlab="time, hour", ylab="Conc",ylim=range(C_deSolve),
    main="Strontium (x = 281m)")
points(times/3600, C_[281,], type="l", lwd=1.5, col="blue2")
legend("topright", legend=c("deSolve", "femR"), 
        lty=c(1,1), lwd=3, col=c("red3", "blue2"))