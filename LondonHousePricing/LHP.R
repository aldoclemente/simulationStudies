# creo grafo direttamente da fdaPDE
data(LNNT)
data(LNHP)
library(fdaPDE)
source("./utils.R")
setwd("LondonHousePricing/")

mesh = as.fdaPDE.SpatialLinesDataFrame(LN.nt) # molto veloce
FEMbasis  = create.FEM.basis(mesh)
#mesh = as.fdaPDE.SpatialLinesDataFrame.shp2graph(LN.nt) # molto lento
#save(mesh, file="data_fdaPDE_only/LHP_mesh.RData")

locs = cbind(LN.prop$X, LN.prop$Y)

mean.purchase = mean(LN.prop$PURCHASE)
dev.sd.purchase = sd(LN.prop$PURCHASE)
PURCHASE.norm = (LN.prop$PURCHASE - mean.purchase)/dev.sd.purchase

#shapiro.test(PURCHASE.norm)
#boxcox(lm(LN.prop$PURCHASE ~ 1))
#PURCHASElog = log(LN.prop$PURCHASE)


locs.proj = projection.points.1.5D(mesh,locs)
plot(mesh, pch=".")
observations = log(LN.prop$PURCHASE)

plot(mesh, type="n")
points(locs, pch=16, cex = 1,col="red")

plot(mesh,type="n")
points(locs.proj, pch=16, cex=1,col="blue")

W = cbind(LN.prop$FLOORSZ, LN.prop$PROF, LN.prop$BATH2)
lambda = 10^seq(from=0.5,to=2.5,by=0.0725)[16:28] #28
output_CPP = smooth.FEM(observations = observations, 
                        locations = locs.proj,
                        FEMbasis = FEMbasis,
                        covariates = W,
                        lambda = lambda,
                        lambda.selection.criterion = "grid",
                        lambda.selection.lossfunction = "GCV",
                        DOF.evaluation = "stochastic")
plot(log10(lambda), output_CPP$optimization$GCV_vector)
sqrt(mean((output_CPP$solution$z_hat - observations)^2))

dataFrame = LN.prop
dataFrame$PURCHASE = log(dataFrame$PURCHASE)
dataFrame@coords = locs.proj

bw.ND = bw.gwr(PURCHASE ~ FLOORSZ + PROF + BATH2, 
               data = dataFrame, 
               approach="AIC", 
               kernel="gaussian")# ,
#               dMat = train_ND)

GWR.ND = gwr.basic(PURCHASE ~ FLOORSZ + PROF + BATH2, 
                     data = dataFrame, 
                     kernel = "gaussian",
                     bw = bw.ND)
sqrt(mean( (GWR.ND$SDF$yhat- observations)^2 ))

# vedo se i punti sono nei nodi e mi salvo gli indici 

# normalizzo la mesh (tanto la mappa 1 a 1 fra nodi e osservazioni rimane uguale)

# 