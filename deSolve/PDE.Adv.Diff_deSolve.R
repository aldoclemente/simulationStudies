#setwd("deSolve/")
library(fdaPDE)
library(ReacTran)
library(latex2exp)
rm(list=ls())
graphics.off()

rmse <- function(x,y){return(sqrt(mean( (x-y)^2)))}

W_ <- 1.
R_ <- 1.
H_ <- 1.
beta_ <- 1

alpha_ <- H_ * beta_ / R_
gamma_ <- pi * W_ / R_  

lambda1 <- -alpha_/2 - sqrt((alpha_/2)^2 + pi^2)
lambda2 <- -alpha_/2 + sqrt((alpha_/2)^2 + pi^2)

p_ <- (1-exp(lambda2))/(exp(lambda1)-exp(lambda2))

exact_solution <- function(points){
  return( -gamma_/pi^2 * ( p_ * exp(lambda1 * points[,1]) + (1 - p_) * exp(lambda2 * points[,1]) - 1. ) * sin(pi * points[,2] ) )
}

# Diffusion > 0 if check = T
forcing <- function(points){
  return(gamma_ * sin(pi * points[,2])) 
}

N = 4*c(5,10,20,40) #,50) # int. nodes N[i]^2

Dy    <- Dx <- -1.   # diffusion coeff, X- and Y-direction
errors.l2 <- rep(0, times = length(N))
h <- rep(0, times = length(N))

for(i in 1:length(N)){
x.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
y.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
grid2D    <- setup.grid.2D(x.grid, y.grid)

h[i] = max(grid2D$dx, grid2D$dy)

D.grid    <- setup.prop.2D(value = Dx, y.value = Dy, grid = grid2D)
v.grid    <- setup.prop.2D(value = +alpha_, y.value=0., grid = grid2D) # on the x direction "minus alpha". the "minus"is embedded into the model...
A.grid    <- setup.prop.2D(value = 1, grid = grid2D)
VF.grid   <- setup.prop.2D(value = 1, grid = grid2D)

C.x.up    <- rep(0, times = N[i])
C.y.up    <- rep(0, times = N[i])

C.x.down  <- rep(0, times = N[i])
C.y.down  <- rep(0, times = N[i])

f = matrix(nrow=N[i], ncol=N[i], data=0)
y.ex = matrix(nrow=N[i], ncol=N[i], data=0)
for(k in 1:N[i]){
  for(l in 1:N[i]){
    f[k,l] = forcing( cbind(grid2D$x.mid[k], grid2D$y.mid[l]) )
    y.ex[k,l] = exact_solution( cbind(grid2D$x.mid[k], grid2D$y.mid[l]))  
  }
}

Adv.Diff <- function (t, y, parms)  {
  
  Y  <- matrix(nrow = N[i], ncol = N[i], data = y)
  
  dY <- tran.2D(C = Y, 
                   C.x.up = C.x.up, C.x.down = C.x.down,
                   C.y.up = C.y.up, C.y.down = C.y.down,
                   grid = grid2D, 
                   D.grid = D.grid, 
                   A.grid = A.grid, 
                   VF.grid = VF.grid, 
                   v.grid = v.grid)$dC #full.output = TRUE, full.check = TRUE
  
  dY <- dY - f 
  
  return (list(dY))
}

y = matrix(data = rep(1,times=N[i]*N[i]))
y = matrix(data = rnorm(N[i]*N[i], mean = 1, sd =0.5))
# lrw (stodes) 
print(system.time(
  Y <- steady.2D(y=y,
                 dimens = c(N[i],N[i]), 
                 time = 0, 
                 func = Adv.Diff, parms=NULL, lrw=1e8) ) )

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

if(beta_ == 2e2){ 
    title_ = "adv_dominated_diff_rates_order_1_deSolve.pdf"
}else{ 
    title_ = "adv_diff_rates_order_1_deSolve.pdf"    
}

pdf(paste(imgdir_,title_,sep=""))
plot(log2(h), log2(errors.l2), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
        ylim = c(min(log2(h), log2(errors.l2)), max(log2(h), log2(errors.l2))+2),
        xlab = TeX("$h$"), ylab="", cex.lab=1.25)
lines(log2(h), log2(h), col = "black", type = "b", pch = 16, lwd = 3, lty =2, cex = 2 )
legend("topleft", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$h$")), 
        col=c("red", "black"), 
        lty = 2, 
        cex=1.25)
dev.off()

pdf(paste(imgdir_, "advection_dominated.pdf",sep=""))
par(mfrow=c(1,2))
image(y.ex, ask = FALSE, main = "exact")
image(y.hat, ask = FALSE, main = "estimate")
dev.off()
#microbenchmark("solution" = {Y <- steady.2D(y=y$y,
#                                            dimens = c(N[i],N[i]), 
#                                            time = 0, 
#                                            func = Adv.Diff, parms=NULL, lrw=1e8) })
