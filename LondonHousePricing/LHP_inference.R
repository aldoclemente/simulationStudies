rm(list=ls())
graphics.off()

setwd("LondonHousePricing/")
filename = "data/2023-01-11-14_14_56/LHP.RData"

load(filename)
library(viridis)
library(colorspace)
source("../utils.R")
source("LNH_utils.R")
library(GWmodel)
data(LNNT)
data(LNHP)

spat.stat.linnet = maptools::as.linnet.SpatialLines(LN.nt)

x.m = mean(spat.stat.linnet$vertices$x)
y.m = mean(spat.stat.linnet$vertices$y)

x.sd = sd(spat.stat.linnet$vertices$x)
y.sd = sd(spat.stat.linnet$vertices$y)

x.norm = (spat.stat.linnet$vertices$x - x.m)/x.sd
y.norm = (spat.stat.linnet$vertices$y - y.m)/y.sd

coords_ = cbind(x.norm, y.norm)

Windows_ = owin(xrange=c(min(x.norm), max(x.norm)), 
                yrange=c(min(y.norm), max(y.norm)) )

spat.stat.linnet = linnet(vertices=as.ppp(coords_, W = Windows_), 
                          edges = cbind(spat.stat.linnet$from, spat.stat.linnet$to),
                          sparse = T)

locs = LN.prop@coords

locs = cbind( (locs[,1]-x.m)/x.sd, (locs[,2]-y.m)/y.sd)
LPP = lpp(locs, spat.stat.linnet)

dataFrame = LN.prop
dataFrame@coords = cbind(LPP$data$x, LPP$data$y)
dataFrame$DATA.IDX = 1:nrow(dataFrame)
dataFrame$PURCHASE = dataFrame$PURCHASE/10^3 # k pounds

x11()
plot(spat.stat.linnet)
points(locs_[,1], locs_[,2], pch=16, col="red3")
points(locs[,1], locs[,2], pch=16, col="blue")

n_cov <- ncol(dataFrame) - 4 # PURCHASE----X-Y-DATA.IDX

observations = dataFrame$PURCHASE
locations = cbind(LPP$data$x, LPP$data$y)
X = as.matrix(dataFrame@data[1:nrow(dataFrame), 2:(n_cov+1)])

inference.data.object <- inferenceDataObjectBuilder(test = "oat", 
                                                    interval = "oat", 
                                                    component = "parametric",
                                                    type = "w",  
                                                    dim = 2,
                                                    n_cov = n_cov)

lambda = 10^seq(from=-3,to=-1.5,length.out=20) 
# NOPE
output_CPP = smooth.FEM(observations = observations, 
                        locations = locations,
                        FEMbasis = FEMbasis,
                        covariates = X,
                        lambda = 1.,
                        lambda.selection.criterion = "newton_fd",
                        lambda.selection.lossfunction = "GCV",
                        DOF.evaluation = "exact",
                        inference.data.object = inference.data.object)

plot(log10(lambda), output_CPP$optimization$GCV_vector)



