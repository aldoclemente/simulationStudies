setwd("deSolve/")
library(fdaPDE)
library(ReacTran)
library(latex2exp)
rm(list=ls())
graphics.off()

rmse <- function(x,y){return(sqrt(mean( (x-y)^2)))}

mesh <- fdaPDE::create.mesh.2D(nodes = matrix(c(0,0,1,0,1,1,0,1), nrow=4,ncol=2, byrow=T))
mesh <- fdaPDE::refine.mesh.2D(mesh=mesh, minimum_angle=30, maximum_area=0.1)
plot(mesh, main=paste("int. nodes = ", nrow(mesh$nodes) - nrow(mesh$nodes[mesh$nodesmarkers,]),sep=""))

for(i in 1:5){
  mesh <- fdaPDE::refine.by.splitting.mesh.2D(mesh = mesh)
  plot(mesh, main = paste("int. nodes = ", nrow(mesh$nodes) - nrow(mesh$nodes[mesh$nodesmarkers,]),sep="")) 
}

#2D
exact <- function(x,y) sin(2*pi*x)*sin(2*pi*y)
f <- function(x,y) 8*pi^2*sin(2*pi*x)*sin(2*pi*y)
Dy    <- Dx <- -1.   # diffusion coeff, X- and Y-direction

N = c(5,10,20,40,80) # int. nodes N[i]^2

errors.l2 <- rep(0, times = length(N))
h <- rep(0, times = length(N))

for(i in 1:length(N)){
  x.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
  y.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
  grid2D    <- setup.grid.2D(x.grid, y.grid)

  h[i] = 1/N[i]^2

  D.grid    <- setup.prop.2D(value = Dx, y.value = Dy, grid = grid2D)
  v.grid    <- setup.prop.2D(value = 0, grid = grid2D)
  A.grid    <- setup.prop.2D(value = 1, grid = grid2D)
  VF.grid   <- setup.prop.2D(value = 1, grid = grid2D)

  C.x.up    <- rep(0, times = N[i])
  C.y.up    <- rep(0, times = N[i])

  C.x.down  <- rep(0, times = N[i])
  C.y.down  <- rep(0, times = N[i])

  forcing = matrix(nrow=N[i], ncol=N[i], data=0)
  y.ex = matrix(nrow=N[i], ncol=N[i], data=0)
  for(k in 1:N[i]){
    for(l in 1:N[i]){
      forcing[k,l] = f(grid2D$x.mid[k], grid2D$y.mid[l])
      y.ex[k,l] = exact(grid2D$x.mid[k], grid2D$y.mid[l])  
    }
  }

  Diff2Db <- function (t, y, parms)  {
  
    Y  <- matrix(nrow = N[i], ncol = N[i], data = y)
  
    dY <- tran.2D(C = Y, 
                     C.x.up = C.x.up, C.x.down = C.x.down,
                     C.y.up = C.y.up, C.y.down = C.y.down,
                     grid = grid2D, 
                     D.grid = D.grid, 
                     A.grid = A.grid, 
                     VF.grid = VF.grid, 
                     v.grid = v.grid)$dC
  
    dY <- dY - forcing 
  
    return (list(dY))
  }

  y = matrix(data = rep(1,times=N[i]*N[i]))
  y = matrix(data = rnorm(N[i]*N[i], mean = 1, sd =0.5))
  # lrw (stodes) 
  print(system.time(
    Y <- steady.2D(y=y,
                   dimens = c(N[i],N[i]), 
                   time = 0, 
                   func = Diff2Db, parms=NULL, lrw=1e8) ) )

  y.hat <- matrix(Y$y, nrow=N[i], ncol=N[i])
  par(mfrow=c(1,2))
  image(y.ex, ask = FALSE, main = "exact")
  image(y.hat, ask = FALSE, main = "estimate")
  errors.l2[i] = rmse(Y$y, as.vector(y.ex))
}

q = log2(errors.l2[1:(length(N)-1)]/errors.l2[2:length(N)])
cat("order = ", q, "\n")

par(mfrow=c(1,1))
plot(log2(h), log2(errors.l2), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
     ylim = c(min(log2(h^2), log2(errors.l2)), max(log2(h), log2(errors.l2))+2),
     xlab = TeX("$h$"), ylab="", cex.lab=1.25)
lines(log2(h), log2(h^2), col = "black", type = "b", pch = 16, lwd = 3, lty =2, cex = 2 )
legend("topleft", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$h^2$")), 
       col=c("red", "black"), 
       lty = 1.5, 
       cex=0.75)
