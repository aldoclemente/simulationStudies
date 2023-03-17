library(fdaPDE)
rm(list=ls())
graphics.off()

N = 100
x_ = seq(0,pi,length.out=N)
y_ = rep(0, times = N)
nodes_ = cbind(x_,y_)

#edges_ = cbind(1:(N-1), 2:N)
edges_ = cbind(2:N, 1:(N-1))

mesh_ = create.mesh.1.5D(nodes=nodes_, edges=edges_, order=2)
FEMbasis = create.FEM.basis(mesh_)

x11()
plot(mesh_, pch=".")

f <- function(x) {
    return(cos(2*x))
}
n_ = 250
locs = cbind(runif(n_,0,pi), rep(0,times=n_))

observations_ = as.matrix( f(locs[,1]))  #+ rnorm(n_,mean=0, sd=0.05*diff(range(f(locs[,1])))) )

SR_PDE = smooth.FEM(observations= observations_ , locations = locs, FEMbasis=FEMbasis, lambda=1e-1)
plot(SR_PDE$fit.FEM)

ordering_ = order(mesh_$nodes[,1])

x11()
plot(mesh_$nodes[ordering_,1], SR_PDE$fit.FEM$coeff[ordering_], lwd=4, type="l",
    ylim=c( min(SR_PDE$fit.FEM$coeff, f(mesh_$nodes[ordering_])),
            max(SR_PDE$fit.FEM$coeff, f(mesh_$nodes[ordering_]))  ))
points(mesh_$nodes[ordering_,1], f(mesh_$nodes[ordering_]), lwd=2, type="l", col="red")

# # # Bif # # #
library(fdaPDE)
rm(list=ls())
graphics.off()

source("../utils.R")

N = 25
x_1 = seq(0,pi/2, length.out=N)

h = pi*cos(pi/4) / (N-1)
x_2 = x_1[N] + seq(0, pi*cos(pi/4), by=h)
x_2 = x_2[-1]
x_3 = x_2

y_1 = rep(0, times=N)
h_y = pi*sin(pi/4) / (N-1)
y_2 = seq(0, pi*sin(pi/4), by=h_y)
y_2 = y_2[-1]
y_3 = -y_2

E1_ = cbind(x_1, y_1)
E2_ = cbind(x_2, y_2)
E3_ = cbind(x_3, y_3)

nodes_ = rbind(E1_, E2_, E3_)

EDGE1 = cbind(1:(N-1), 2:N)
#EDGE2 = cbind(N:(2*N-2), (N+1): (2*N-1))
EDGE2 = cbind((N+1): (2*N-1), N:(2*N-2))
EDGE3 = cbind( c(N, seq(2*N,(3*N-3),by=1)), seq(2*N, (3*N-2), by=1))

edges_ = rbind(EDGE1, EDGE2, EDGE3)

mesh_ = create.mesh.1.5D(nodes= nodes_, edges= edges_, order=2)
FEMbasis = create.FEM.basis(mesh_)

x11()
plot(mesh_, pch=".")

f <- function(x,y){
    ret_ = matrix(0, nrow = length(x), ncol=1)
    
    for(i in 1:length(x)){
        if(abs(y[i]) < .Machine$double.eps)
            ret_[i]=cos(x[i])
        else 
            ret_[i] = sin( (x[i]-x_1[N])/(2*cos(pi/4)))

        }
    return(ret_)
}

FEM_ = FEM(coeff= f(mesh_$nodes[,1], mesh_$nodes[,2]), FEMbasis)
plot(FEM_)

ordering_1 = 1:N
ordering_2 = (N+1):(2*N-1)
ordering_3 = (2*N):(3*N-2)
x11()
plot(mesh_$nodes[ordering_1,1], FEM_$coeff[ordering_1,], type="l" ,lwd=3, main="cos(x)", xlab="", ylab="")
points(mesh_$nodes[ordering_1,1], cos(seq(0, pi/2, length.out=length(ordering_1))), type="l", lwd=2, col="red")
(FEM_$coeff[ordering_1[N],] - FEM_$coeff[ordering_1[N-1],]) / h
(f(pi/2 - 10*.Machine$double.eps , 0.) - f(pi/2, 0)) / (-10*.Machine$double.eps)


x11()
plot(c(mesh_$nodes[N,1], mesh_$nodes[ordering_2,1]), c(FEM_$coeff[N,], FEM_$coeff[ordering_2,]), type="l" ,lwd=4, main="sin(x/2)", xlab="", ylab="")
points(c(mesh_$nodes[N,1], mesh_$nodes[ordering_2,1]), sin(seq(0, pi, length.out=(length(ordering_2)+1))/2), type="l", lwd=2, col="red")

x11()
plot(c(mesh_$nodes[N,1], mesh_$nodes[ordering_3,1]), c(FEM_$coeff[N,], FEM_$coeff[ordering_3,]), type="l" ,lwd=4, main="sin(x/2)", xlab="", ylab="")
points(c(mesh_$nodes[N,1], mesh_$nodes[ordering_3,1]), sin(seq(0, pi, length.out=(length(ordering_3)+1))/2), type="l", lwd=2, col="red")

(f(pi/2 + 10*.Machine$double.eps, 10*.Machine$double.eps) - f(pi/2, 10*.Machine$double.eps)) / (10*.Machine$double.eps)


x11()
plot(mesh_$nodes[ordering_3,1], SR_PDE$fit.FEM$coeff[ordering_3,], type="l" ,lwd=4, 
     ylim=c( min(SR_PDE$fit.FEM$coeff, sin(seq(0, pi, length.out=(length(ordering_3)))/2 )  ),
             max(SR_PDE$fit.FEM$coeff, sin(seq(0, pi, length.out=(length(ordering_3)))/2 )  )) )
points(mesh_$nodes[ordering_3,1], sin(seq(0, pi, length.out=(length(ordering_3)))/2), type="l", lwd=2, col="red")

spat.stat <- as.spatstat.linnet.fdaPDE(mesh_)

n_ = 1000
locs = runiflpp(n_, spat.stat)

locs_ = cbind(locs$data$x, locs$data$y)

observations_ = as.matrix( f(locs_[,1], locs_[,2]) )

SR_PDE = smooth.FEM(observations= observations_ , locations = locs_, FEMbasis=FEMbasis, lambda=1e-2)
plot(SR_PDE$fit.FEM)

ordering_1 = 1:N
ordering_2 = (N+1):(2*N-1)
ordering_3 = (2*N):(3*N-2)
x11()
plot(mesh_$nodes[ordering_1,1], SR_PDE$fit.FEM$coeff[ordering_1,], type="l" ,lwd=2)

dF_h = eval.FEM(FEM= SR_PDE$fit.FEM, 
                locations= cbind( c( pi/2 - 10^5*.Machine$double.eps , pi/2 + 10^5*.Machine$double.eps, pi/2 + 10^5*.Machine$double.eps),
                                  c(0., 10^5*.Machine$double.eps, 10^5*.Machine$double.eps) ) )
dF_  = eval.FEM(FEM= SR_PDE$fit.FEM, 
                locations= cbind( c( pi/2, pi/2, pi/2),
                                  c(0., 0, 0) ) )

H_ =  c(-10^5*.Machine$double.eps, 10^5*.Machine$double.eps, 10^5*.Machine$double.eps)

sum( (dF_h - dF_)/ H_ ) 

SR_PDE$fit.FEM$coeff[N,]

x11()
plot(mesh_$nodes[ordering_1,1], SR_PDE$fit.FEM$coeff[ordering_1,], type="l" ,lwd=4, 
     ylim=c( min(SR_PDE$fit.FEM$coeff, cos(seq(0, pi/2, length.out=(length(ordering_1))))),
             max(SR_PDE$fit.FEM$coeff, cos(seq(0, pi/2, length.out=(length(ordering_1))))) ))
points(mesh_$nodes[ordering_1,1], cos(seq(0, pi/2, length.out=(length(ordering_1)))), type="l", lwd=2, col="red")

x11()
plot(mesh_$nodes[ordering_2,1], SR_PDE$fit.FEM$coeff[ordering_2,], type="l" ,lwd=4, 
     ylim=c( min(SR_PDE$fit.FEM$coeff, sin(seq(0, pi, length.out=(length(ordering_2)))/2 )  ),
             max(SR_PDE$fit.FEM$coeff, sin(seq(0, pi, length.out=(length(ordering_2)))/2 )  )) )
points(mesh_$nodes[ordering_2,1], sin(seq(0, pi, length.out=(length(ordering_2)))/2), type="l", lwd=2, col="red")

x11()
plot(mesh_$nodes[ordering_3,1], SR_PDE$fit.FEM$coeff[ordering_3,], type="l" ,lwd=4, 
     ylim=c( min(SR_PDE$fit.FEM$coeff, sin(seq(0, pi, length.out=(length(ordering_3)))/2 )  ),
             max(SR_PDE$fit.FEM$coeff, sin(seq(0, pi, length.out=(length(ordering_3)))/2 )  )) )
points(mesh_$nodes[ordering_3,1], sin(seq(0, pi, length.out=(length(ordering_3)))/2), type="l", lwd=2, col="red")

# # #  # # # 
library(femR)
library(fdaPDE)
rm(list=ls())
source("/home/aldo/Desktop/fdaPDE/wrappers/femR/tests/utils.R")
{
N = 25
x_1 = seq(0,pi/2, length.out=N)

h = pi*cos(pi/4) / (N-1)
x_2 = x_1[N] + seq(0, pi*cos(pi/4), by=h)
x_2 = x_2[-1]
x_3 = x_2

y_1 = rep(0, times=N)
h_y = pi*sin(pi/4) / (N-1)
y_2 = seq(0, pi*sin(pi/4), by=h_y)
y_2 = y_2[-1]
y_3 = -y_2

E1_ = cbind(x_1, y_1)
E2_ = cbind(x_2, y_2)
E3_ = cbind(x_3, y_3)

nodes_ = rbind(E1_, E2_, E3_)

EDGE1 = cbind(1:(N-1), 2:N)
#EDGE2 = cbind(N:(2*N-2), (N+1): (2*N-1))
EDGE2 = cbind((N+1): (2*N-1), N:(2*N-2))
EDGE3 = cbind( c(N, seq(2*N,(3*N-3),by=1)), seq(2*N, (3*N-2), by=1))

edges_ = rbind(EDGE1, EDGE2, EDGE3)
}

mesh = create.mesh.1.5D(nodes= nodes_, edges= edges_)
x11()
plot(mesh, pch=".")
graph_ <- list(nodes= as.matrix(mesh$nodes), edges= mesh$edges, elements= mesh$edges, neigh= sparseNeigh(mesh), boundary= mesh$nodesmarkers)

PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)
PDE <- new(PDE_1_5D_isotropic_ORDER_1, graph_)

PDE$set_PDEparameters(PDE_parameters)

f <- function(x,y){
    ret_ = matrix(0, nrow = length(x), ncol=1)
    
    for(i in 1:length(x)){
        if(abs(y[i]) < .Machine$double.eps)
            ret_[i]=cos(x[i])
        else 
            ret_[i] = sin( (x[i]-x_1[N])/(2*cos(pi/4)))

        }
    return(ret_)
}
