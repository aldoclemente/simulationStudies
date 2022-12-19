library(fdaPDE)
library(ggplot2)

nodes = matrix(c(0,0,1,0,(1+sqrt(2)/2), sqrt(2)/2, (1+sqrt(2)/2), -sqrt(2)/2), nrow=4,ncol=2, byrow = T)
edges = matrix(c(1,2,2,3,2,4), nrow=3,ncol=2, byrow=T)

mesh = create.mesh.1.5D(nodes, edges)

R_plot_mesh.ggplot = function(mesh, alpha = 1, line.size=0.75,
                              points_ = NULL,
                              points.size = 7,
                              points.color="black",
                              title = "",
                              ratio = diff(range(mesh$nodes[,1]))/diff(range(mesh$nodes[,2])),
                              title.size = 26,
                              legend.pos = "right",
                              line.color="black",
                              points2.color= NULL, points2_ = NULL)
{
  x=vector(mode="double")
  y=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  MyTheme <- theme(
    axis.text = element_text(size=(title.size-2)),
    axis.title = element_text(size=title.size),
    title = element_text(size=title.size),
    legend.text = element_text(size=(title.size-6)),
    legend.key.size = unit(1,"cm"),
    legend.position = legend.pos
  )
  
    data=data.frame(x,y,grp.nodes)
    
    num_points = nrow(points_)
    
    x.points = points_[,1]
    y.points = points_[,2]
    data.points = data.frame(x.points, y.points)
    
    PLOT <- ggplot(data=NULL) + 
        geom_point(data=data, aes(x=x,y=y,group=grp.nodes),
                   alpha=0.0) + 
        geom_line(data=data, aes(x=x,y=y,group=grp.nodes), 
                  size=line.size, alpha = 0.5) +
        geom_point(data=data.points,aes(x=x.points,y=y.points ), color = points.color,
                   size=points.size) +
        labs(x="",y="",color="", title=title) + 
        coord_fixed(ratio=ratio) + 
        theme_void() +
        MyTheme +   
        theme(plot.title = element_text(hjust=0.5),
              legend.title = element_blank(),
              axis.title = element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              legend.key.size = unit(1,"cm"),
              legend.key.width = unit(0.5,"cm"))
    
    if(!is.null(points2_)){
      num_points.2 = nrow(points2_)
      
      x.points2 = points2_[,1]
      y.points2 = points2_[,2]
      data.points = data.frame(x.points2, y.points2)
      
      PLOT <- PLOT + 
               geom_point(data=data.points,aes(x=x.points2,y=y.points2 ), color = points2.color,
                          size=points.size)
  }
    return(PLOT)
}

mesh.ref = refine.by.splitting.mesh.1.5D(mesh)

pdf("network_notation.pdf")
R_plot_mesh.ggplot(mesh=mesh, points_ = mesh$nodes, 
                   title = bquote(G[]))
R_plot_mesh.ggplot(mesh, points_ = mesh$nodes[mesh$nodesmarkers, ], 
                   points.color = "red2", 
                   title = bquote(W[B]) )
R_plot_mesh.ggplot(mesh, points_ = matrix( mesh$nodes[!mesh$nodesmarkers, ], nrow=1,ncol=2), 
                   points.color = "green4", 
                   title = bquote(W[I]) )

R_plot_mesh.ggplot(mesh, points_ = mesh$nodes, 
                   title = "")
R_plot_mesh.ggplot(mesh, points_ = mesh$nodes[mesh$nodesmarkers, ], 
                   points.color = "red2", 
                   title = "" )
R_plot_mesh.ggplot(mesh, points_ = matrix( mesh$nodes[!mesh$nodesmarkers, ], nrow=1,ncol=2), 
                   points.color = "green4", 
                   title = "" )

R_plot_mesh.ggplot(mesh.ref, points_ = mesh.ref$nodes,
                   points.color = "red2",
                   title = "")

R_plot_mesh.ggplot(mesh.ref, points_ = mesh.ref$nodes, 
                   title = "")
R_plot_mesh.ggplot(mesh.ref, points_ = mesh.ref$nodes[mesh.ref$nodesmarkers, ], 
                   points.color = "red2", 
                   title = "" )

R_plot_mesh.ggplot(mesh.ref, points_ = matrix( mesh.ref$nodes[!mesh.ref$nodesmarkers, ], nrow=nrow(mesh.ref$nodes[!mesh.ref$nodesmarkers, ]),ncol=2), 
                   points.color = "green4", 
                   title = "" )

R_plot_mesh.ggplot(mesh.ref, 
                   points_ = matrix( mesh.ref$nodes[!mesh.ref$nodesmarkers, ], nrow=nrow(mesh.ref$nodes[!mesh.ref$nodesmarkers, ]),ncol=2), 
                   points.color = "green4", 
                   points2_ = mesh.ref$nodes[mesh.ref$nodesmarkers, ],
                   points2.color = "red2",
                   title = "" )


dev.off()


x.points2 = mesh.ref$nodes[mesh.ref$nodesmarkers,1]
y.points2 = mesh.ref$nodes[mesh.ref$nodesmarkers,2]
data.points = data.frame(x.points2, y.points2)

esempio + geom_point(data=data.points,aes(x=x.points2,y=y.points2 ), color = "red2",
                  size=7)


### London Observations ###
library(viridis)
library(ggplot2)
R_plot_observations.ggplot = function(mesh, points_, points.value,
                                      alpha = 1, line.size=0.25,
                              points.size = 1,
                              points.palette= viridis,
                              title = "",
                              ratio = diff(range(mesh$nodes[,1]))/diff(range(mesh$nodes[,2])),
                              title.size = 26,
                              legend.pos = "right",
                              line.color="black")
{
  x=vector(mode="double")
  y=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  MyTheme <- theme(
    axis.text = element_text(size=(title.size-2)),
    axis.title = element_text(size=title.size),
    title = element_text(size=title.size),
    legend.text = element_text(size=(title.size-6)),
    legend.key.size = unit(1,"cm"),
    legend.position = legend.pos
  )
  
  data=data.frame(x,y,grp.nodes)
  
  num_points = nrow(points_)
  
  x.points = points_[,1]
  y.points = points_[,2]
  data.points = data.frame(x.points, y.points, points.value)
  
  p <- points.palette(n=56,alpha=1)
  
  PLOT <- ggplot(data=NULL) + 
    geom_point(data=data, aes(x=x,y=y,group=grp.nodes),
               alpha=0.0) + 
    geom_line(data=data, aes(x=x,y=y,group=grp.nodes), 
              size=line.size, alpha = 0.75) +
    geom_point(data=data.points,aes(x=x.points,y=y.points, color=points.value ),
               size=points.size) +
    labs(x="",y="",color="", title=title) + 
    coord_fixed(ratio=ratio) + 
    scale_color_gradientn(colours=p)+
    theme_void() +
    MyTheme +   
    theme(plot.title = element_text(hjust=0.5),
          legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.key.size = unit(1,"cm"),
          legend.key.width = unit(0.5,"cm"))
  
  return(PLOT)
}


library(fdaPDE)
library(shp2graph)

source("utils.R")
source("LondonHousePricing/LNH_utils.R")
data(LNNT)
data(LNHP)

spat.stat.linnet = maptools::as.linnet.SpatialLines(LN.nt)

mu_x = mean(spat.stat.linnet$vertices$x)
mu_y = mean(spat.stat.linnet$vertices$y)
coords_ = cbind(spat.stat.linnet$vertices$x-mu_x, 
                spat.stat.linnet$vertices$y-mu_y)/10^3
Windows_ = owin(xrange=c((min(spat.stat.linnet$vertices$x)-mu_x)/10^3,
                         (max(spat.stat.linnet$vertices$x)-mu_x)/10^3), 
                yrange=c((min(spat.stat.linnet$vertices$y)-mu_y)/10^3,
                         (max(spat.stat.linnet$vertices$y)-mu_y)/10^3)
)

spat.stat.linnet = linnet(vertices=as.ppp(coords_, W = Windows_), 
                          edges = cbind(spat.stat.linnet$from, spat.stat.linnet$to),
                          sparse = T)

locs = LN.prop@coords
#which.duplicated = which(duplicated(locs))
which.duplicated = which(duplicated(LN.prop@data))

locs = locs[-which.duplicated, ]
locs = cbind(locs[,1]-mu_x, locs[,2]-mu_y)/10^3

LPP = lpp(locs, spat.stat.linnet)

# fdaPDE mesh 
nodes = cbind(spat.stat.linnet$vertices$x, spat.stat.linnet$vertices$y)
edges = cbind(spat.stat.linnet$from, spat.stat.linnet$to)
mesh = create.mesh.1.5D(nodes=nodes,
                        edges=edges)
FEMbasis = create.FEM.basis(mesh)

dataFrame = LN.prop[-which.duplicated,]
dataFrame@coords = cbind(LPP$data$x, LPP$data$y)
dataFrame$DATA.IDX = 1:nrow(dataFrame)
dataFrame$PURCHASE = dataFrame$PURCHASE # k pounds

pdf("prova2.pdf")
R_plot_observations.ggplot(mesh, points_ = dataFrame@coords, points.value = dataFrame$PURCHASE,
                           alpha = 1, line.size=0.125,
                           points.size = 0.5)
R_plot_observations.ggplot(mesh, points_ = dataFrame@coords, points.value = dataFrame$PURCHASE,
                           alpha = 1, line.size=0.125,
                           points.size = 0.75)

R_plot_observations.ggplot(mesh, points_ = dataFrame@coords, points.value = dataFrame$PURCHASE,
                           alpha = 1, line.size=0.125,
                           points.size = 0.5)
R_plot_observations.ggplot(mesh, points_ = dataFrame@coords, points.value = dataFrame$PURCHASE,
                           alpha = 1, line.size=0.0725,
                           points.size = 1)


dev.off()