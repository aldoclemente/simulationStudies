setwd("deSolve/")
library(fdaPDE)
library(ReacTran)
library(latex2exp)

rm(list=ls())
graphics.off()

rmse <- function(x,y){return(sqrt(mean( (x-y)^2)))}

exact_solution <- function(points){
  return( exp(points[,1] + points[,2]) )
}

forcing <- function(points){
  return(0.) 
}

N = c(5,10,20,40) #,50) # int. nodes N[i]^2

Dy    <- Dx <- -1.   # diffusion coeff, X- and Y-direction
errors.l2 <- rep(0, times = length(N))
h <- rep(0, times = length(N))

for(i in 1:length(N)){
  x.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
  y.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
  grid2D    <- setup.grid.2D(x.grid, y.grid)
  
  h[i] = max(grid2D$dx, grid2D$dy)
  
  D.grid    <- setup.prop.2D(value = Dx, y.value = Dy, grid = grid2D)
  v.grid    <- setup.prop.2D(value = 0., y.value=0., grid = grid2D) # on the x direction "minus alpha". the "minus"is embedded into the model...
  A.grid    <- setup.prop.2D(value = 1, grid = grid2D)
  VF.grid   <- setup.prop.2D(value = 1, grid = grid2D)
  
  C.x.up    <- exact_solution(cbind(rep(x.grid$x.up, times=N[i]), y.grid$x.mid))
  C.y.up    <- exact_solution(cbind(x.grid$x.mid, rep(y.grid$x.up, times=N[i])))
  
  C.x.down  <- exact_solution(cbind(rep(x.grid$x.down, times=N[i]), y.grid$x.mid))
  C.y.down  <- exact_solution(cbind(x.grid$x.mid, rep(y.grid$x.down, times=N[i])))
  
  f = matrix(nrow=N[i], ncol=N[i], data=0)
  y.ex = matrix(nrow=N[i], ncol=N[i], data=0)
  for(k in 1:N[i]){
    for(l in 1:N[i]){
      f[k,l] = forcing( cbind(grid2D$x.mid[k], grid2D$y.mid[l]) )
      y.ex[k,l] = exact_solution( cbind(grid2D$x.mid[k], grid2D$y.mid[l]))  
    }
  }
  
  Reac.Diff <- function (t, y, parms)  {
    
    Y  <- matrix(nrow = N[i], ncol = N[i], data = y)
    
    dY <- tran.2D(C = Y, 
                  C.x.up = C.x.up, C.x.down = C.x.down,
                  C.y.up = C.y.up, C.y.down = C.y.down,
                  grid = grid2D, 
                  D.grid = D.grid, 
                  A.grid = A.grid, 
                  VF.grid = VF.grid, 
                  v.grid = v.grid)$dC
    
    dY <- dY + 2*Y 
    
    return (list(dY))
  }
  
  y = matrix(data = rep(1,times=N[i]*N[i]))
  y = matrix(data = rnorm(N[i]*N[i], mean = 1, sd =0.5))
  # lrw (stodes) 
  print(system.time(
    Y <- steady.2D(y=y,
                   dimens = c(N[i],N[i]), 
                   time = 0, 
                   func = Reac.Diff, parms=NULL, lrw=1e8) ) )
  
  y.hat <- matrix(Y$y, nrow=N[i], ncol=N[i])
  par(mfrow=c(1,2))
  image(y.ex, ask = FALSE, main = "exact")
  image(y.hat, ask = FALSE, main = "estimate")
  errors.l2[i] = rmse(Y$y, as.vector(y.ex))
}

q = log2(errors.l2[1:(length(N)-1)]/errors.l2[2:length(N)])
cat("order = ", q, "\n")

imgdir_ = "imgs/"
if(!dir.exists(imgdir_))
  dir.create(imgdir_)

title_ = "reac_diff_rates_order_1_deSolve.pdf"

pdf(paste(imgdir_,title_,sep=""))
par(mfrow=c(1,1))
plot(log2(h), log2(errors.l2), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
     ylim = c(min(log2(h), log2(errors.l2)), max(log2(h), log2(errors.l2))+2),
     xlab = TeX("$h$"), ylab="", cex.lab=1.25)
lines(log2(h), log2(h^2), col = "black", type = "b", pch = 16, lwd = 3, lty =2, cex = 2 )
legend("topleft", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$h^2$")), 
       col=c("red", "black"), 
       lty = 2, 
       cex=1.25)
dev.off()
