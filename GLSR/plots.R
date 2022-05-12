library(ggplot2)
R_plot_mesh.ggplot = function(mesh, 
                              alpha = 1, 
                              line.size=0.5,
                              mesh.2D = NULL, 
                              alpha.2D=0.9,
                              points_ = NULL,
                              mu = NULL,
                              title = "",
                              ratio = diff(range(mesh$nodes[,1]))/diff(range(mesh$nodes[,2])) )
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
    axis.text = element_text(size=24),
    axis.title = element_text(size=26),
    title = element_text(size=26),
    legend.text = element_text(size=20),
    legend.key.size = unit(1,"cm"),
    legend.position = "right"
  )
  
  if( is.null(mesh.2D) && is.null(points_) ){
    
    data=data.frame(x,y,grp.nodes)
    
    ggplot(data=data, aes(x=x,y=y,group=grp.nodes)) + 
      geom_point(alpha=0.0) + 
      geom_line(size=line.size)+
      labs(x="",y="",color="", title=title) +  
      coord_fixed(ratio=ratio) + 
      theme_void() +
      theme(plot.title = element_text(hjust=0.5),
            title = element_text(size=26),
            legend.title = element_blank(),
            axis.title = element_blank(),
            legend.text = element_text(size=20),
            legend.key.size = unit(1,"cm"),
            legend.key.width = unit(0.5,"cm"))
  
  }else if(is.null(mesh.2D) & !is.null(points_)){
    p = jet.col(n=128, alpha=0.8)  
    data=data.frame(x,y,grp.nodes)
    num_points = nrow(points_)
    coef.points = mu
    
    x.points = points_[,1]
    y.points = points_[,2]
    data.points = data.frame(x.points, y.points, coef.points)
    
    ggplot(data=data) + 
      geom_point( aes(x=x,y=y,group=grp.nodes), 
                 alpha=0.0) + 
      geom_line( aes(x=x,y=y,group=grp.nodes), 
                size=1, alpha = 0.5) +
      geom_point(data=data.points,aes(x=x.points,y=y.points, color=coef.points),
                 size=3) +
      labs(x="",y="",color="", title=title) + 
      scale_color_gradientn(colours=p)+ 
      coord_fixed(ratio=ratio) + 
      theme_void() +
      theme(plot.title = element_text(hjust=0.5),
            title = element_text(size=26),
            legend.title = element_blank(),
            axis.title = element_blank(),
            legend.text = element_text(size=20),
            legend.key.size = unit(1,"cm"),
            legend.key.width = unit(0.5,"cm"))
    
    
   
  }else if(!is.null(mesh.2D) & is.null(points_)){
  
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
    
    ggplot(data=data) + 
      geom_point( aes(x=x,y=y,group=grp.nodes), 
                 alpha=0.0) + 
      geom_line( aes(x=x,y=y,group=grp.nodes), 
                size=1, alpha = 0.5) +
      geom_line(data=data.2D, aes(x=x.2D, y=y.2D, group=grp.nodes.2D), 
                size=1, alpha=0.2) +
      geom_line(data=data.2D.segments, aes(x=x.2D, y=y.2D, group=grp.nodes.2D), 
                size=1, color="red", alpha=0.3) +
      labs(x="",y="",color="", title=title) + 
      scale_color_gradientn(colours=p)+ 
      coord_fixed(ratio=ratio) + 
      theme_void() +
      theme(plot.title = element_text(hjust=0.5),
            title = element_text(size=26),
            legend.title = element_blank(),
            axis.title = element_blank(),
            legend.text = element_text(size=20),
            legend.key.size = unit(1,"cm"),
            legend.key.width = unit(0.5,"cm"))
    
    
  
  }else{
    
    p = jet.col(n=128, alpha=0.8)  
    
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
    
    ggplot(data=data) + 
      geom_point( aes(x=x,y=y,group=grp.nodes), 
                 alpha=0.0) + 
      geom_line( aes(x=x,y=y,group=grp.nodes), 
                size=1, alpha = 0.5) +
      geom_line(data=data.2D, aes(x=x.2D, y=y.2D, group=grp.nodes.2D), 
                size=1, alpha=0.2) +
      geom_line(data=data.2D.segments, aes(x=x.2D, y=y.2D, group=grp.nodes.2D), 
                size=1, color="red", alpha=0.3) +
      geom_point(data=data.points,aes(x=x.points,y=y.points, color=coef.points),
                 size=3) +
      labs(x="",y="",color="", title=title) + 
      scale_color_gradientn(colours=p)+ 
      coord_fixed(ratio=5) + 
      theme_void() +
      theme(plot.title = element_text(hjust=0.5),
            title = element_text(size=26),
            legend.title = element_blank(),
            axis.title = element_blank(),
            legend.text = element_text(size=20),
            legend.key.size = unit(1,"cm"),
            legend.key.width = unit(0.5,"cm"))
    
    
  }
}
