library(fdaPDE)
library(ggplot2)
library(viridis)
source("Auxiliary/R_plot_graph.ggplot2.R")

R_plot_graph.ggplot2.2<-function(FEM,
                                 title="", 
                                 line.size=0.75,
                                 legend.pos="right",
                                 color.min=min(FEM$coeff), 
                                 color.max=max(FEM$coeff),
                                 ratio = diff(range(mesh$nodes[,1]))/diff(range(mesh$nodes[,2]) ),
                                 palette = jet.col,
                                 num.color = 28,
                                 title.size = 26,
                                 color.gradient = T, 
                                 points_=NULL,
                                 points.col=NULL,
                                 points.size=7,
                                 is.hat.function=FALSE){
  
  mesh=FEM$FEMbasis$mesh
  x=vector(mode="double")
  y=vector(mode="double")
  coef=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    coef=append(coef, rep((FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]])/2,times=2) )  
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  p <- palette(n=num.color,alpha=1)
  
  MyTheme <- theme(
    axis.text = element_text(size=(title.size-2)),
    axis.title = element_text(size=title.size),
    title = element_text(size=title.size),
    legend.text = element_text(size=(title.size-6)),
    legend.key.size = unit(1,"cm"),
    legend.position = legend.pos,
  )
  
  data=data.frame(x,y,grp.nodes,coef)
  MyTheme <- theme(
    axis.text = element_text(size=(title.size-2)),
    axis.title = element_text(size=title.size),
    title = element_text(size=title.size),
    legend.text = element_text(size=(title.size-6)),
    legend.key.size = unit(1,"cm"),
    legend.position = legend.pos,
    legend.key.height=unit(2.5,"cm"),
    legend.key.width = unit(0.5,"cm")
  )
  
  x.points = points_[,1]
  y.points = points_[,2]
  data.points = data.frame(x.points, y.points)
  
  if( !is.null(points_) ){
  gplot <- ggplot() + 
    geom_point(data=data, aes(x=x,y=y,group=grp.nodes),
               alpha=0.0) + 
    geom_line(data=data, aes(x=x,y=y,group=grp.nodes,color=coef),
              size=line.size)+
    geom_point(data=data.points,aes(x=x.points,y=y.points),
               color = points.col,
               size=points.size) +
    scale_color_gradientn(colours=p, limits = c(color.min, color.max))+
    labs(x="",y="",color="", title=title) +  
    coord_fixed(ratio=ratio) + 
    theme_void() +
    MyTheme +
    theme(plot.title = element_text(hjust=0.5),
          legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.key.width = unit(0.3,"cm"),
          legend.position = legend.pos,
          legend.key.height=unit(2.5,"cm"))
  }else{
    
    gplot <- ggplot() + 
      geom_point(data=data, aes(x=x,y=y,group=grp.nodes),
                 alpha=0.0) + 
      geom_line(data=data, aes(x=x,y=y,group=grp.nodes,color=coef),
                size=line.size)+
      scale_color_gradientn(colours=p, limits = c(color.min, color.max))+
      labs(x="",y="",color="", title=title) +  
      coord_fixed(ratio=ratio) + 
      theme_void() +
      MyTheme +
      theme(plot.title = element_text(hjust=0.5),
            legend.title = element_blank(),
            axis.title = element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            legend.key.width = unit(0.3,"cm"),
            legend.position = legend.pos,
            legend.key.height=unit(2.5,"cm"))
    
  }
  if(!color.gradient){
    gplot <- gplot + theme(legend.position = "none") 
  }
  if(is.hat.function)
    gplot <- gplot + scale_color_gradientn(colours=p, limits = c(color.min, color.max),
                    breaks=c(0,1),labels = c("0", "1"))
  
  return (gplot)
}


nodes = matrix(c(0,0,1,0,(1+sqrt(2)/2), sqrt(2)/2, (1+sqrt(2)/2), -sqrt(2)/2), nrow=4,ncol=2, byrow = T)
edges = matrix(c(1,2,2,3,2,4), nrow=3,ncol=2, byrow=T)

mesh = create.mesh.1.5D(nodes, edges)
mesh.ref = refine.by.splitting.mesh.1.5D(mesh)

mesh.plot = refine.mesh.1.5D(mesh.ref, delta=0.01)

mesh.ref$nodes
mesh.plot$nodes[1:7,]

Hat = function(x,node=2){
  result = matrix(0, nrow=nrow(x$nodes), ncol=1)
  h = 0.5
  for(i in 1: nrow(x$nodes)){
    dist = norm(x$nodes[i,] - x$nodes[node,], type="2")
    if(dist<=h){
      result[i] = (1 - 2*dist) 
    }
  }
  return(result)
}

coeff = Hat(mesh.plot)

FEMbasis = create.FEM.basis(mesh.plot)
FEM_ = FEM(coeff, FEMbasis)
plot(FEM_)

points.col = viridis(2,begin=1, end=0)

LINE.SIZE = 2.5
FEM.ggplot.1 <- R_plot_graph.ggplot2.2(FEM_, line.size = LINE.SIZE, palette = viridis, 
                       points_ = mesh.ref$nodes[c(2,c(5:7)), ], 
                       points.col = c(points.col[1], rep(points.col[2], times=3)),
                       is.hat.function = T)
x11()
FEM.ggplot.1


FEM.ggplot.2 <- R_plot_graph.ggplot2.2(FEM_, line.size = LINE.SIZE, palette = viridis, 
                                     points_ = mesh.ref$nodes[c(2,c(5:7)), ], 
                                     points.col = c(points.col[1], rep(points.col[2], times=3)),
                                     color.gradient = F)

FEM.ggplot.2

pdf("hat_functions")
  R_plot_graph.ggplot2.2(FEM_, line.size = LINE.SIZE, palette = viridis, 
                                       points_ = mesh.ref$nodes[c(2,c(5:7)), ], 
                                       points.col = c(points.col[1], rep(points.col[2], times=3)),
                                       is.hat.function = T)

  R_plot_graph.ggplot2.2(FEM_, line.size = LINE.SIZE, palette = viridis, 
                                       points_ = mesh.ref$nodes[c(2,c(5:7)), ], 
                                       points.col = c(points.col[1], rep(points.col[2], times=3)),
                                       color.gradient = F, is.hat.function = T)

  R_plot_graph.ggplot2.2(FEM_, line.size = LINE.SIZE, palette = viridis, is.hat.function = T)
  R_plot_graph.ggplot2.2(FEM_, line.size = LINE.SIZE, palette = viridis,
                       color.gradient = F, is.hat.function = T)
  dev.off()

coeff = Hat(mesh.plot)

FEMbasis = create.FEM.basis(mesh.plot)
FEM_ = FEM(coeff, FEMbasis)
plot(FEM_)

line.size = c(1,1.5,1.75,2,2.5,3,3.5)
for( i in 1:length(line.size)){
  pdf(paste("hat_degree_3_",line.size[i],".pdf",sep=""))
  print(R_plot_graph.ggplot2.2(FEM_, line.size = line.size[i], palette = viridis, is.hat.function = T))
  print(R_plot_graph.ggplot2.2(FEM_, line.size = line.size[i], palette = viridis,
                       color.gradient = F, is.hat.function = T))
  dev.off()
}

library(cowplot)
library(grid)
library(gridExtra)

HAT_PLOT <- R_plot_graph.ggplot2.2(FEM_, line.size = line.size[i], palette = viridis, is.hat.function = T)
HAT_legend <- cowplot::get_legend(HAT_PLOT + theme(legend.key.height = unit(2,"cm")))

grid.newpage()
grid.draw(HAT_legend)
pdf("hat_legend.pdf", height=5, width = 1)
grid.draw(HAT_legend)
dev.off()

coeff = Hat(mesh.plot,node = 5)
FEMbasis = create.FEM.basis(mesh.plot)
FEM_ = FEM(coeff, FEMbasis)

for( i in 1:length(line.size)){
  pdf(paste("hat_degree_2_",line.size[i],".pdf",sep=""))
  print(R_plot_graph.ggplot2.2(FEM_, line.size = line.size[i], palette = viridis, is.hat.function = T))
  print(R_plot_graph.ggplot2.2(FEM_, line.size = line.size[i], palette = viridis,
                       color.gradient = F, is.hat.function = T))
  dev.off()
}
