rm(list=ls())
library(femR)
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
dx = 0.01
riverlen = 1
Nb = riverlen/dx # number of interval
                 # Nb + 1 mesh nodes

xgrid = seq(from=0, to=riverlen, length.out=(Nb+1)) # mesh nodes
edges = cbind( seq(from=1,to=Nb-1,by=1), seq(from=2,to=Nb,by=1))

neigh = neigh1D(edges)
boundary = matrix(0,nrow=(Nb+1), ncol=1)
boundary[1] = boundary[Nb+1] = 1

River <- list(nodes= as.matrix(xgrid), edges= edges, elements= edges, neigh= neigh, boundary= boundary)
PDE_parameters = list("diffusion" = 1., 
                       "transport" = as.matrix(0.),
                        "reaction" = 0.)

PDE <- new(PDE_1D_isotropic_ORDER_2, River)
PDE$set_PDEparameters(PDE_parameters)