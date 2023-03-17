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
Nb = riverlen/dx # number of interval
                 # Nb + 1 mesh nodes
#xgrid = seq(from=0, to=riverlen, length.out=(Nb+1)) # mesh nodes
xgrid = setup.grid.1D(N=Nb, L= riverlen)

# parameters 
area  = 0.52; area2 = 1.56; ased = 1.5
Vol   = area*dx
disp  = 0.40
D     = disp*area/dxK1 = (qlat/area + alpha + ased/area*sorp)
K2 = -(alpha*area/area2 + sorpS)

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
Yini = c(rep(0.13, times=2*Nb), rep(0., times=Nb))
times = seq(from=8.25*3600, to= 24*3600, by=120)

Dyna = ode.1D(y = Yini, func = RiverModel, nspec=3, times= times)

image(Dyna, grid= xgrid$x.mid)

C_deSolve = Dyna[,2:(3*Nb+1)]
x11()
filled.contour(x=times/3600, y=xgrid$x.mid, C_deSolve[,1:Nb], main="C")

x11()
filled.contour(x=times/3600, y=xgrid$x.mid, C_deSolve[,(Nb+1):(2*Nb)], main="C_LS")

x11()
filled.contour(x=times/3600, y=xgrid$x.mid, C_deSolve[,(2*Nb+1):(3*Nb)], main="C_sed")


pdf("RiverUvas_deSolve.pdf")
filled.contour(x=times/3600, y=xgrid$x.mid, C_deSolve[,1:Nb], main="C")
filled.contour(x=times/3600, y=xgrid$x.mid, C_deSolve[,(Nb+1):(2*Nb)], main="C_LS")
filled.contour(x=times/3600, y=xgrid$x.mid, C_deSolve[,(2*Nb+1):(3*Nb)], main="C_sed")
dev.off()