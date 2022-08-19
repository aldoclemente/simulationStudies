auxiliary = function(x, y, seg, tp, sigma= 0.125, 
               nodes.lpp = ppp(x = mesh$nodes[,1], y = mesh$nodes[,2], 
               window = owin(xrange = c(min(mesh$nodes[,1]),max(mesh$nodes[,1])),
                             yrange = c(min(mesh$nodes[,2]),max(mesh$nodes[,2])))),
               source=63)
{ 
  PP = ppp(x = x, y = y, window = nodes.lpp$window)
  ND = crossdist(nodes.lpp, PP)
  
  return( 1/sqrt(2*pi*sigma^2) * exp(-ND[source,]^2/(2*sigma^2)) )
  
}

setting <-function(network = "" ){
  
  if(network == "ontario"){
    data("ORN")
    
    mesh = as.fdaPDE.SpatialLinesDataFrame(ORN.nt)
    mesh = normalize_mesh(mesh)
    
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    
    spat.stat.linnet = as.spatstat.linnet.fdaPDE(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spat.stat.linnet = spat.stat.linnet)
    
    return(res)
  }else if( network == "estevan"){
    
    data("ERN_OSM_correct")
    
    mesh = as.fdaPDE.SpatialLinesDataFrame(ERN_OSM_cor.nt)
    mesh = normalize_mesh(mesh)
    
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    
    spat.stat.linnet = as.spatstat.linnet.fdaPDE(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spat.stat.linnet = spat.stat.linnet)
  }else if( network == "simplenet"){
    mesh = create.mesh.1.5D(nodes = cbind(simplenet$vertices$x, simplenet$vertices$y),
                            edges = cbind(simplenet$from, simplenet$to))
    
    mesh = refine.mesh.1.5D(mesh, delta = 0.025)
    FEMbasis = create.FEM.basis(mesh)
    
    nnodes = nrow(mesh$nodes)
    
    spat.stat.linnet = as.spatstat.linnet.fdaPDE(mesh)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               spat.stat.linnet = spat.stat.linnet)
    
    return(res)
  }
  
}
