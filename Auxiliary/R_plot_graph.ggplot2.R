library(ggplot2)
R_plot_graph.ggplot2<-function(FEM, title="", 
                               line.size=0.5,
                               legend.pos="right"){
  
  mesh=FEM$FEMbasis$mesh
  x=vector(mode="double")
  y=vector(mode="double")
  coef=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    coef=append(coef, rep( (FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]]) /2,times=2) )  
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  p <- jet.col(n=128,alpha=0.8)
  
  MyTheme <- theme(
    axis.text = element_text(size=24),
    axis.title = element_text(size=26),
    title = element_text(size=26),
    legend.text = element_text(size=20),
    legend.key.size = unit(1.,"cm"),
    legend.position = legend.pos
  )
  
  data=data.frame(x,y,grp.nodes,coef)
  ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
    geom_point(alpha=0.0) + 
    geom_line(aes(color=coef), size=line.size) +
    scale_color_gradientn(colours=p) + 
    labs(x="",y="",color="", title=title) + 
    theme(panel.background = element_rect(fill="white",color="gray50"),
          plot.title = element_text(hjust=0.5))+
    MyTheme
  
  
}

R_plot_graph.a.sym.ggplot2<-function(FEM,title="", 
                                     line.size=0.5,
                                     legend.pos="right"){
  
  mesh=FEM$FEMbasis$mesh
  x=vector(mode="double")
  y=vector(mode="double")
  coef=vector(mode="double")
  grp.nodes=vector(mode="integer")
  
  num_edges= dim(mesh$edges)[1]
  for(e in 1:num_edges){
    x = append(x, c(mesh$nodes[ mesh$edges[e,1], 1], mesh$nodes[ mesh$edges[e,2], 1]))
    y = append(y, c(mesh$nodes[ mesh$edges[e,1], 2], mesh$nodes[ mesh$edges[e,2], 2]))
    if( abs(FEM$coeff[mesh$edges[e,1]])<= 10*.Machine$double.eps 
        ||
        abs(FEM$coeff[mesh$edges[e,2]]) <= 10*.Machine$double.eps){
      coef = append(coef, c(0.0, 0.0) )
    }else{
      coef= append(coef, rep( (FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]]) /2,times=2) )  
    }
    grp.nodes = append(grp.nodes, rep(e,times=2))
  }
  
  p <- jet.col(n=128,alpha=0.8)
  
  MyTheme <- theme(
    axis.text = element_text(size=24),
    axis.title = element_text(size=26),
    title = element_text(size=26),
    legend.text = element_text(size=20),
    legend.key.size = unit(1,"cm"),
    legend.position = legend.pos
  )
  
  data=data.frame(x,y,grp.nodes,coef)
  ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
    geom_point(alpha=0.0) + 
    geom_line(aes(color=coef), size=line.size) +
    scale_color_gradientn(colours=p) + 
    labs(x="",y="",color="", title=title) + 
    theme(panel.background = element_rect(fill="white",color="gray50"),
          plot.title = element_text(hjust=0.5))+
    MyTheme
  
  
}

# FEM, oggetto FEM da plottare (fem = FEM(coeff, FEMbasis))
# title, titolo plot
# legend.pos, posizione legenda
# color.min, valore del colore più piccolo da inserire nella scala colore
# color.max, valore del colore più grande da inserire nella scala colore
# ratio, ratio fra asse x e asse y (consiglio lasciare valore di default)
# palette, color palette (es jet.col, viridis, magma, inferno, plasma)
# num.col, numero di colori della palette da generare
# title.size, dimensione titolo plot
R_plot_graph.ggplot2.2<-function(FEM,
                                 title="", 
                                 line.size=0.5,
                                 legend.pos="right",
                                 color.min=min(FEM$coeff), 
                                 color.max=max(FEM$coeff),
                                 ratio = diff(range(mesh$nodes[,1]))/diff(range(mesh$nodes[,2]) ),
                                 palette = jet.col,
                                 num.color = 28,
                                 title.size = 26,
                                 color.gradient = T){
  
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
      legend.position = legend.pos
    )
    
    gplot <- ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
    geom_point(alpha=0.0) + 
    geom_line(aes(color=coef), size=line.size)+
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
          legend.key.width = unit(0.5,"cm"),
          legend.position = legend.pos)
    
    if(!color.gradient){
       gplot <- gplot + theme(legend.position = "none") 
     }
    
  return (gplot)
}


# alternativa 

R_plot_graph.R <- function(FEM, title="", line.size="1"){
  
  mesh = FEM$FEMbasis$mesh
  coef = matrix(0, nrow=nrow(mesh$edges), ncol=1)
  for(e in 1: nrow(mesh$edges)){
    if( abs(FEM$coeff[mesh$edges[e,1]])<= 10*.Machine$double.eps 
        ||
        abs(FEM$coeff[mesh$edges[e,2]]) <= 10*.Machine$double.eps){
      coef[e,1] = 0.0
      
    }else{
      coef[e,1] = (FEM$coeff[mesh$edges[e,1]] + FEM$coeff[mesh$edges[e,2]]) /2   
    } 
  }
  

  heat <- jet.col(100)
  plot(mesh$nodes, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", type="n", main=title, cex.main=3)
  segments(mesh$nodes[mesh$edges[,1],1], mesh$nodes[mesh$edges[,1],2],
           mesh$nodes[mesh$edges[,2],1], mesh$nodes[mesh$edges[,2],2], 
           col=heat[round(99*(coef[,1]-min(coef[,1]))/diff(range(coef)))+1], 
           lwd=line.size)

}

R_plot_mesh.ggplot = function(mesh, alpha = 1, line.size=0.5,
                              mesh.2D = NULL, alpha.2D=0.9,
                              points_ = NULL,
                              points.size = 3,
                              mu = NULL,
                              color.min=min(mu), 
                              color.max=max(mu),
                              title = "",
                              ratio = diff(range(mesh$nodes[,1]))/diff(range(mesh$nodes[,2])),
                              palette = jet.col,
                              num.col = 128,
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
  
  if( is.null(mesh.2D) && is.null(points_) ){

  data=data.frame(x,y,grp.nodes)
  
  ret <- ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
    geom_point(alpha=0.0) + 
    geom_line(size=line.size, color=line.color)+
    labs(x="",y="",color="", title=title) +  
    coord_fixed(ratio=ratio) + 
    theme_void() +
    MyTheme +
    theme(plot.title = element_text(hjust=0.5),
          legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          legend.key.width = unit(0.5,"cm"))
  return(ret)
  }else if(is.null(mesh.2D) && !is.null(points_)){
    if(num.col>2){
     p = palette(n=num.col, alpha=alpha)
    }else{
     p = c("red3", "green4")
    }
    data=data.frame(x,y,grp.nodes)
    
    num_points = nrow(points_)
    coef.points = mu
    
    x.points = points_[,1]
    y.points = points_[,2]
    data.points = data.frame(x.points, y.points, coef.points)
    if(num.col>2){ 
    ggplot(data=NULL) + 
      geom_point(data=data, aes(x=x,y=y,group=grp.nodes), 
                 alpha=0.0) + 
      geom_line(data=data, aes(x=x,y=y,group=grp.nodes), 
                size=line.size) +
      geom_point(data=data.points,aes(x=x.points,y=y.points, color=coef.points),
                 size=points.size) +
      labs(x="",y="",color="", title=title) + 
      scale_color_gradientn(colours=p, limits = c(color.min, color.max))+ 
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
    }else{
      ggplot(data=NULL) + 
        geom_point(data=data, aes(x=x,y=y,group=grp.nodes), 
                   alpha=0.0) + 
        geom_line(data=data, aes(x=x,y=y,group=grp.nodes), 
                  size=line.size, alpha = 0.5) +
        geom_point(data=data.points,aes(x=x.points,y=y.points, color=coef.points),
                   size=points.size) +
        labs(x="",y="",color="", title=title) + 
        scale_color_gradientn(colours=p)+
        coord_fixed(ratio=ratio) + 
        theme_void() +
        MyTheme +
        theme(plot.title = element_text(hjust=0.5),
              legend.title = element_blank(),
              axis.title = element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              legend.key.width = unit(0.5,"cm"))
      
    }
  }else{
   p = palette(n=num.col, alpha=0.8)  
   
   data=data.frame(x,y,grp.nodes)
   
   x.2D=vector(mode="double")
   y.2D=vector(mode="double")
   grp.nodes.2D=vector(mode="integer")
   
   num_edges= dim(mesh.2D$edges)[1]
   for(e in 1:num_edges){
     x.2D = append(x.2D, c(mesh.2D$nodes[ mesh.2D$edges[e,1], 1], mesh.2D$nodes[ mesh.2D$edges[e,2], 1]))
     y.2D = append(y.2D, c(mesh.2D$nodes[ mesh.2D$edges[e,1], 2], mesh.2D$nodes[ mesh.2D$edges[e,2], 2]))
     grp.nodes.2D = append(grp.nodes.2D, rep(e,times=2))
   }
   data.2D = data.frame(x.2D, y.2D, grp.nodes.2D)
   
   x.2D=vector(mode="double")
   y.2D=vector(mode="double")
   grp.nodes.2D=vector(mode="integer")
   num_segments= dim(mesh.2D$segments)[1]
   for(e in 1:num_segments){
     x.2D = append(x.2D, c(mesh.2D$nodes[ mesh.2D$segments[e,1], 1], mesh.2D$nodes[ mesh.2D$segments[e,2], 1]))
     y.2D = append(y.2D, c(mesh.2D$nodes[ mesh.2D$segments[e,1], 2], mesh.2D$nodes[ mesh.2D$segments[e,2], 2]))
     grp.nodes.2D = append(grp.nodes.2D, rep(e,times=2))
   }
   data.2D.segments = data.frame(x.2D, y.2D, grp.nodes.2D)
   
  num_points = nrow(points_)
  coef.points = mu
  
   x.points = points_[,1]
   y.points = points_[,2]
   data.points = data.frame(x.points, y.points, coef.points)
   
   ggplot(data=NULL) + 
     geom_point(data=data, aes(x=x,y=y,group=grp.nodes), 
                alpha=0.0) + 
     geom_line(data=data, aes(x=x,y=y,group=grp.nodes), 
               size=line.size, alpha = 0.5) +
     geom_point(data=data.points,aes(x=x.points,y=y.points, color=coef.points),
                size=points.size) +
     labs(x="",y="",color="", title=title) + 
     scale_color_gradientn(colours=p)+ 
     coord_fixed(ratio=ratio) + 
     theme_void() +
     MyTheme + 
     theme(plot.title = element_text(hjust=0.5),
           legend.title = element_blank(),
           axis.title = element_blank(),
           axis.text.x=element_blank(),
           axis.text.y=element_blank(),
           legend.key.size = unit(1,"cm"))
  
   
   }
}

R_plot_point_pattern.ggplot <- function(mesh, 
                                       alpha = 1, 
                                       line.size=0.5,
                              points_ = NULL,
                              points.size = 3,
                              title = "",
                              ratio = diff(range(mesh$nodes[,1]))/diff(range(mesh$nodes[,2])),
                              col = "red3",
                              title.size = 26,
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
    legend.key.size = unit(1,"cm"))
  
  if(is.null(points_) ){
    
    data=data.frame(x,y,grp.nodes)
    
    ret <- ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
      geom_point(alpha=0.0) + 
      geom_line(size=line.size, color=line.color)+
      labs(x="",y="",color="", title=title) +  
      coord_fixed(ratio=ratio) + 
      theme_void() +
      MyTheme +
      theme(plot.title = element_text(hjust=0.5),
            legend.title = element_blank(),
            axis.title = element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            legend.key.width = unit(0.5,"cm"))
    return(ret)
  }else{
  
    data=data.frame(x,y,grp.nodes)
    
    num_points = nrow(points_)
    
    x.points = points_[,1]
    y.points = points_[,2]
    data.points = data.frame(x.points, y.points)
    
    ret <- ggplot(data=NULL) + 
        geom_point(data=data, aes(x=x,y=y,group=grp.nodes), 
                   alpha=0.0) + 
        geom_line(data=data, aes(x=x,y=y,group=grp.nodes), 
                  size=line.size) +
        geom_point(data=data.points,aes(x=x.points,y=y.points),
                   color=col,
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
            legend.key.width = unit(0.5,"cm"))
    
    return(ret)
    }
  }
  

