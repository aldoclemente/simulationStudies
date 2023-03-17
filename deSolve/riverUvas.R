rm(list=ls())
library(Matrix)

neigh1D <- function(edges){
    
    neigh = matrix(0, nrow=nrow(edges), ncol=2)
    for(i in 1:nrow(edges))
        for(j in 1:2)
            if( (i == 1 && j ==2) || ( i == nrow(edges) && j == 1))
                neigh[i,j] = -1
            else if(  (i != nrow(edges) && j==1)){   
                neigh[i,j] = i+1
            }                                             
            else if(  (i != 1 && j==2)){
                neigh[i,j] = i-1
            }

    return(neigh)
}

init_t = function(mesh_list, PDE_parameters, BC_t, f1_, f2_){

    PDE <- new(PDE_1D_isotropic_ORDER_1, mesh_list)
    PDE$set_PDEparameters(PDE_parameters)

    dirichletBC = matrix(0, nrow=nrow(mesh_list$nodes))
    dirichletBC[1] = BC_t
    dirichletBC[nrow(dirichletBC)] = BC_t # (?)    
    PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes = PDE$get_quadrature_nodes()
    F1_ = as.matrix( f1_ * rep(1, times=nrow(quadrature_nodes)) )
    PDE$set_forcingTerm(F1_)
    
    res = PDE$solve() # dummy
    RHS1 = PDE$force()
    Mass = res$Mass
    Stiff = res$Stiff

    F2_ = as.matrix( f2_ * rep(1, times=nrow(quadrature_nodes)) )
    PDE$set_forcingTerm(F2_)
    invisible(PDE$solve())
    RHS2 = PDE$force()
    
    return(list(Mass = Mass, Stiff = Stiff, RHS1 = RHS1, RHS2 = RHS2))

}

Assemble = function(mesh_list, coeffs_, matrices_, U, k){
    nnodes = nrow(mesh_list$nodes)
    K_ = Matrix(0.0,nrow=3*nnodes, ncol=3*nnodes, sparse=T)

    # "K_11"
    K_[1:nnodes,1:nnodes]                            = (1/coeffs_$dt * matrices_$Mass +  matrices_$Stiff)
    
    # "K_22"
    K_[(nnodes+1):(2*nnodes), (nnodes+1):(2*nnodes)]     = coeffs_$k_22 * matrices_$Mass
    
    # # "K_33"
    K_[(2*nnodes+1):(3*nnodes), (2*nnodes+1):(3*nnodes)] = coeffs_$k_33 * matrices_$Mass
    
    # "K_12"
    K_[1:nnodes,(nnodes+1):(2*nnodes)]                 =  coeffs_$k_12 * matrices_$Mass
    # "K_13"nnodes = nrow(mesh_list$nodes)
    K_ = Matrix(0.0,nrow=3*nnodes, ncol=3*nnodes, sparse=T)

    # "K_11"
    K_[1:nnodes,1:nnodes]                            = (1/coeffs_$dt * matrices_$Mass +  matrices_$Stiff)
    
    # "K_22"
    K_[(nnodes+1):(2*nnodes), (nnodes+1):(2*nnodes)]     = coeffs_$k_22 * matrices_$Mass
    
    # # "K_33"
    K_[(2*nnodes+1):(3*nnodes), (2*nnodes+1):(3*nnodes)]       = coeffs_$k_13 * matrices_$Mass
    
    # "K_21"
    K_[1:nnodes,(nnodes+1):(2*nnodes)]                 = coeffs_$k_21 * matrices_$Mass 

    # "K_23"
    K_[(nnodes+1):(2*nnodes),(2*nnodes+1):(3*nnodes)]    = coeffs_$k_23 * matrices_$Mass

    # "K_31"
    K_[(2*nnodes+1):(3*nnodes),(nnodes+1):(2*nnodes)]    = coeffs_$k_31 * matrices_$Mass

    # "K_32"
    K_[(2*nnodes+1):(3*nnodes),(nnodes+1):(2*nnodes)]    = coeffs_$k_32 * matrices_$Mass

    RHS = Matrix(0, nrow=3*nnodes, ncol=1, sparse =T)
    RHS[1:nnodes]              = matrices_$RHS1 + 1/coeffs_$dt * matrices_$Mass %*% U[1:nnodes,k]

    RHS[(nnodes+1):(2*nnodes)]   = matrices_$RHS2 + 1/coeffs_$dt * matrices_$Mass %*% U[(nnodes+1):(2*nnodes),k]
    
    RHS[(2*nnodes+1):(3*nnodes)] =                 1/coeffs_$dt * matrices_$Mass %*% U[(nnodes+1):(2*nnodes),k]
    return(list(K_ = K_, RHS_ = RHS))
}

Solve = function(mesh_list, PDE_parameters,
                 U, U_t, times, coeffs_,
                f1_=qlat*Clat/area, f2_=sorpS*CsB){

        
    for(i in 2:length(times)){
    Cup = U_t(times[i])

    matrices_ = init_t( mesh_list = mesh_list, PDE_parameters = PDE_parameters,
                        BC_t = Cup, f1_ = f1_, f2_ = f2_)
        
    
    
    assembler_ = Assemble(mesh_list =mesh_list, coeffs_ = coeffs_, matrices_=matrices_, U=U, k=(i-1))
    #tmp = solve(assembler_$K_, assembler_$RHS_)
    U[,i] = as.vector(solve(assembler_$K_, assembler_$RHS_))
    }

    return(list(U=U))
}

library(femR)

# model grid
dx = 1
riverlen = 669
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
                        "reaction" = K1)
dt_ = 120
coeffs_ = list(dt = dt_, 
                                              k_12 = (-alpha),      k_13 = (ased/area * sorp), 
                k_21 = (-alpha * area/area2), k_22 = (1/dt_ - K2), k_23 = 0., 
                k_31 = (-sorp),               k_32 = 0. ,         k_33 = (1/dt_ + sorp))
C_ = matrix(nrow=3*(Nb+1), ncol= length(times))
C_[,1] = Yini

result = Solve(mesh_list = River, PDE_parameters = PDE_parameters ,
               U = C_, U_t = Cfun, times=times, coeffs_=coeffs_)


x11()
filled.contour(x=times, y=xgrid, t(result$U[1:(Nb+1),]))
