library(fdaPDE)
library(spatstat)
library(maptools)
library(shp2graph)
library(igraph)
library(rgeos)
library(GWmodel)
library(diffusionMaps)


#' Utility to convert fdaPDE Linear Network mesh into SpatialLinesDataFrame. 
#' It is used to solve GWR, GGWR or Kriging over Linear Networks.
#' @param mesh, fdaPDE Linear Network mesh
#' 
#' @return SpatialLinesDf, It contains the length of each edges in the network.
#'   
as.SpatialLinesDataFrame.fdaPDE <- function(mesh){
  
  nodes = mesh$nodes
  edges = mesh$edges
  
  x0 = vector(mode="numeric",length=nrow(edges))
  y0 = x0
  x1 = x0
  y1 = x0
  
  for(e in 1:nrow(edges)){
    x0[e] = nodes[edges[e,1],1]
    y0[e] = nodes[edges[e,1],2]
    x1[e] = nodes[edges[e,2],1]
    y1[e] = nodes[edges[e,2],2]
  }
  
  psp_ =psp(x0,y0,x1,y1, window=owin(xrange=c(min(x0,x1),max(x0,x1)), 
                               yrange=c(min(y0,y1),max(y0,y1))) ) 
  # maptools 
  SpatialLines_ = as.SpatialLines.psp(psp_)
  
  # rgeos
  df <- data.frame(len = sapply(1:length(SpatialLines_), function(i) gLength(SpatialLines_[i, ])))
  rownames(df) <- sapply(1:length(SpatialLines_), function(i) SpatialLines_@lines[[i]]@ID)
  
  # object to return
  SpatialLinesDf <- SpatialLinesDataFrame(SpatialLines_, data = df)
  
  return(SpatialLinesDf)
}

#' Utility to convert SpatialLinesDataFrame into fdaPDE Linear Network mesh. 
#' It is used to solve GWR, GGWR or Kriging over Linear Networks.
#' @param SpLinesDf, SpatialLinesDataFrame
#' 
#' @return mesh, It contains the length of each edges in the network.
#'   
as.fdaPDE.SpatialLinesDataFrame <- function(SpLinesDf){
  #maptools
  
  spat.stat.linnet = maptools::as.linnet.SpatialLines(SpLinesDf)
  nodes = cbind(spat.stat.linnet$vertices$x, spat.stat.linnet$vertices$y)
  edges = cbind(spat.stat.linnet$from, spat.stat.linnet$to)
  
  mesh = create.mesh.1.5D(nodes, edges)
  
  return(mesh)
}

#' Utility to convert SpatialLinesDataFrame into fdaPDE Linear Network mesh. 
#' It is used to solve GWR, GGWR or Kriging over Linear Networks.
#' @param SpLinesDf, SpatialLinesDataFrame
#' 
#' @return mesh, It contains the length of each edges in the network.
#'   
as.fdaPDE.SpatialLinesDataFrame.shp2graph<-function(SpLinesDf){
  
  rt.NEL = readshpnw(SpLinesDf, 
                    ELComputed = TRUE,
                    Detailed = TRUE,
                    ea.prop = names(SpLinesDf))
  g = nel2igraph(nodelist = rt.NEL[[2]], 
                 edgelist = rt.NEL[[3]],
                 weight   = rt.NEL[[4]])
  
  DataFrame = as_data_frame(g, what="both")
  mesh = create.mesh.1.5D(nodes = as.matrix(DataFrame$vertices),
                          edges = as.matrix(DataFrame$edges[,1:2]))

  return(mesh)
}

fdaPDEnodes2Spnodes<-function(mesh.fdaPDE, sp.edges){
  
  fdaPDE.edges = mesh.fdaPDE$edges
  nedges = nrow(fdaPDE.edges)
  nnodes = nrow(mesh.fdaPDE$nodes)
  
  # 0-vector initilization
  fdaPDE2Sp = vector(mode="integer", length=nnodes)
  set_ = list()
  for(e in 1:nedges){
    
    if( !(0 %in% fdaPDE2Sp)){
      break
    }
    
    flag1 = sp.edges[e,1] %in% fdaPDE2Sp
    flag2 = sp.edges[e,2] %in% fdaPDE2Sp
    
    if( !flag1 ){
      set_ = append(set_, sp.edges[e,1])
      fdaPDE2Sp[fdaPDE.edges[e,1]] = sp.edges[e,1]
    }
    
    if( !flag2 ){
      set_ = append(set_, sp.edges[e,2])
      fdaPDE2Sp[fdaPDE.edges[e,2]] = sp.edges[e,2]
    }
  }
  
  return(fdaPDE2Sp)
}

#' Compute the Network distance matrix assuming that locations are a subset of mesh
#' nodes. It is used to solve GWR, GGWR or Kriging over Linear Networks.
#' @param mesh, fdaPDE Linear Network mesh
#' @param nodes.idx, integer vector containing locations indices
#' 
#' @return ND, network distance matrix
compute_NetworkDist_matrix <- function(mesh, nodes.idx){
  
  SpatialLinesDf = as.SpatialLinesDataFrame.fdaPDE(mesh)
  
  # shp2graph
  rt.NEL = readshpnw(SpatialLinesDf)
  igr = nel2igraph(rt.NEL[[2]], rt.NEL[[3]])
  
  # indices map from fdaPDE to Sp  
  fdaPDE2Sp = fdaPDEnodes2Spnodes(mesh, rt.NEL[[3]][,2:3])
  
  # igraph
  ND = distances(igr, v = fdaPDE2Sp[nodes.idx], 
                      to = fdaPDE2Sp[nodes.idx], 
                      weights = SpatialLinesDf$len)
  
  return(ND)
}

#' Compute the Euclidean distance matrix assuming that locations are a subset of mesh
#' nodes. It is used to solve GWR, GGWR or Kriging over Linear Networks.
#' @param mesh, fdaPDE Linear Network mesh
#' @param nodes.idx, integer vector containing locations indices
#' 
#' @return ED, euclidean distance matrix
computed_EuclideanDist_matrix <- function(mesh, nodes.idx){
  
  coords_ = mesh$nodes[nodes.idx,]
  
  # GWmodel
  ED = gw.dist(dp.locat = coords_)
  
  return(ED)
}

normalize_mesh <-function(mesh){
  
  x.m = mean(mesh$nodes[,1])
  y.m = mean(mesh$nodes[,2])
  
  x.sd = sd(mesh$nodes[,1])
  y.sd = sd(mesh$nodes[,1])
  
  x.norm = (mesh$nodes[,1] - x.m)/x.sd
  y.norm = (mesh$nodes[,2] - y.m)/y.sd

  mesh.norm = create.mesh.1.5D(nodes=cbind(x.norm,y.norm), edges = mesh$edges)  
  return(mesh.norm)
}

#' Compute the Network distance matrix assuming that locations are a subset of mesh
#' nodes. It is used to solve GWTR over Linear Networks.
#' @param mesh, fdaPDE Linear Network mesh
#' @param nodes.idx, integer vector containing locations indices
#' @param time.vec, a vector of time tags for each observations
#' 
#' @return ND, network distance matrix
compute_NetworkDist_matrix.time <- function(mesh, nodes.idx, time.vec){
  
  SpatialLinesDf = as.SpatialLinesDataFrame.fdaPDE(mesh)
  
  # shp2graph
  rt.NEL = readshpnw(SpatialLinesDf)
  igr = nel2igraph(rt.NEL[[2]], rt.NEL[[3]])
  
  # indices map from fdaPDE to Sp  
  fdaPDE2Sp = fdaPDEnodes2Spnodes(mesh, rt.NEL[[3]][,2:3])
  
  # igraph
  ND = distances(igr, v = fdaPDE2Sp[nodes.idx], 
                 to = fdaPDE2Sp[nodes.idx], 
                 weights = SpatialLinesDf$len)
  
  # GWmodel
  x.loc = rep(mesh$nodes[nodes.idx,1], times=length(time.vec))
  y.loc = rep(mesh$nodes[nodes.idx,2], times=length(time.vec))
  space.locations = cbind(x.loc, y.loc)
  time.locations = rep(time.vec, each=nrow(space.locations))
  
  ND.time = st.dist(dp.locat = space.locations, 
                    obs.tv= time.locations, 
                    s.dMat=ND)
  
  return(ND.time)
}

#' Compute the Euclidean distance matrix assuming that locations are a subset of mesh
#' nodes. It is used to solve GWR, GGWR or Kriging over Linear Networks.
#' @param mesh, fdaPDE Linear Network mesh
#' @param nodes.idx, integer vector containing locations indices
#' @param time.vec, a vector of time tags for each observations
#' 
#' @return ED, euclidean distance matrix
computed_EuclideanDist_matrix.time <- function(mesh, nodes.idx, time.vec){
  
  nnodes = length(nodes.idx)
  
  # GWmodel
  x.loc = rep(mesh$nodes[nodes.idx,1], times=length(time.vec))
  y.loc = rep(mesh$nodes[nodes.idx,2], times=length(time.vec))
  space.locations = cbind(x.loc, y.loc)
  time.locations = rep(time.vec, each=nrow(mesh$nodes[nodes.idx,]))
  
  ED = gw.dist(dp.locat = mesh$nodes[nodes.idx,])
  
  ED.time = st.dist(dp.locat = space.locations,
                    obs.tv = time.locations,
                    s.dMat = ED.space)
  
  
  return(ED.time)
}

#' Utility to convert fdaPDE Linear Network mesh into lattice objects. 
#' It is used to solve GWR, GGWR or Kriging over Linear Networks.
#' @param mesh, fdaPDE 1.5D mesh
#' 
#' @return a list containing a matrix of lattice nodes and a sparse adjacency matrix.
#'   
as.lattice.fdaPDE<-function(mesh){
  
  # diffusionMaps nb. z coords are set to 0. 
  nodes.lattice = cbind(mesh$nodes, rep(0,times=nnodes))
  
  # Sparse matrix (spam package)
  adj_matrix = spam(0,nrow=nnodes, ncol=nnodes)
  nedges = nrow(mesh$edges)
  for( e in 1:nedges){
    sx = mesh$edges[e,1]
    dx = mesh$edges[e,2]
    
    adj_matrix[sx, dx] = 1
    adj_matrix[dx, sx] = 1
  }
  adj_matrix = spam::cleanup(adj_matrix)
  
  out = list(nodes.lattice = nodes.lattice,
             adj_matrix = adj_matrix)
  
  return(out)
}

#' Utility to convert fdaPDE Linear Network mesh into spatstat linnet. 
#' @param mesh, fdaPDE Linear Network mesh
#' 
#' @return SpatialLinesDf, It contains the length of each edges in the network.
#'   
as.spatstat.linnet.fdaPDE <- function(mesh){
  
  vertices = ppp(x = mesh$nodes[,1], 
                 y = mesh$nodes[,2], 
                 window = owin(xrange = c(min(mesh$nodes[,1]),max(mesh$nodes[,1])),
                               yrange = c(min(mesh$nodes[,2]),max(mesh$nodes[,2]))))
  spat.stat.linnet = linnet(vertices = vertices, edges = mesh$edges)
 
 return(spat.stat.linnet)
}

#' Utility to convert fdaPDE Linear Network mesh into spatstat linnet. 
#' @param mesh, fdaPDE Linear Network mesh
#' 
#' @return SpatialLinesDf, It contains the length of each edges in the network.
#'   
as.fdaPDE.spatstat.linnet <- function(spatstat.linnet){
  
  nodes = cbind(spat.stat.linnet$vertices$x, spat.stat.linnet$vertices$y)
  edges = cbind(spat.stat.linnet$from, spat.stat.linnet$to)
  
  mesh = create.mesh.1.5D(nodes, edges)
  
  return(mesh)
}
