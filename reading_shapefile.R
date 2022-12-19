library(sf)

path = "C:/Users/Aldo/Downloads/inc_strad_milano_2016/inc_strad_milano_2016.shp"

# leggi shapefile 
incidenti <- st_read(path)

# printing
incidenti

# 
SHPs <- list.files(path="C:/Users/Aldo/Downloads/DBT2012_STRATO01_E0", pattern=".shp", full.names=T)

# A010101 Area di circolazione veicolare, 
# A010102 Area di circolazione pedonale, 
# A010103 Area di circolazione ciclabile, 
# A010104 Area stradale, 
# A010105 Viabilità mista secondaria, 
# L010105 Viabilità mista secondaria, 
# L010107 Elemento stradale, 
# L010202 Elemento ferroviario, 
# L010204 Elemento tranviario, 
# L010206 Elemento di metropolitana,
# L010210 Binario industriale, 
# P010108 Giunzione stradale, 
# P010203 Giunzione ferroviaria,
# P010205 Giunzione tranviaria,
# P010207 Giunzione di metropolitana, 
# LIM010102 Contorno area di circolazione pedonale, 
# LIM010103 Contorno area di circolazione ciclabile, 
# LIM010104 Contorno area stradale, 
# LIM010105 Contorno viabilità mista secondaria. 
SHPs.objects <- sapply(SHPs, function(x) st_read(x))

network.sf <- as.data.frame(SHPs.objects[[6]])

source("utils.R")

# converting reference system
network.sf <- st_transform(SHPs.objects[[6]], crs = st_crs(incidenti))

# drop the z coordinates
network.sf <- st_zm(network.sf, drop=TRUE, what ="ZM")

# Sp Object
SpLinesDataFrame_ <- as_Spatial(network.sf)
class(SpLinesDataFrame_)

# fdaPDE mesh "it takes a while"
mesh = as.fdaPDE.SpatialLinesDataFrame.shp2graph( SpLineDataFrame_)
milanoLN = mesh

save(milanoLN, file="data/MilanoLN.RData")

###

load("data/MilanoLN.RData")
nrow(milanoLN$nodes)
nrow(milanoLN$edges)
nrow(incidenti)

source("Auxiliary/R_plot_graph.ggplot2.R")

pdf("milanoIncidenti.pdf")
plot(milanoLN, pch=".")
points(incidenti$longitudin, incidenti$latitudine, 
       pch=16, col="red3",cex=0.25)

R_plot_point_pattern.ggplot(mesh=milanoLN)
R_plot_point_pattern.ggplot(mesh=milanoLN, 
                            points_ = cbind(incidenti$longitudin, incidenti$latitudine),
                            points.size = 1)

dev.off() 
