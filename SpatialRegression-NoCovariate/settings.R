f.sin = function(ND,source=63,sigma=1e-5){
  
  #sigma = 0.125/24
  nnodes = nrow(ND)
  res = vector(mode="numeric", length=nnodes)
  
  source.1 = source
  distances.1 = ND[source.1,]
  other.1 = vector(mode="integer")
  for(i in 1:nnodes){
    
    res[i] = sin(2*pi*distances.1[i]/sigma)*cos(2*pi*distances.1[i]/sigma)
    
  }
  
  return(res)
}

f.exp = function(ND,source=63,sigma=0.125){
  
  nnodes = nrow(ND)
  res = vector(mode="numeric", length=nnodes)
  
  source.1 = source
  distances.1 = ND[source.1,]
  other.1 = vector(mode="integer")
  for(i in 1:nnodes ){
    res[i] = exp(-distances.1[i]/sigma)  
  }
  
  return(res)
}

fun <-function(f = "exp"){
  if(f == "exp")
    return (f.exp)
  if(f == "sin")
    return (f.sin)
  
}

setting <-function(network = "" ){
  
  if(network == "ontario"){
  data("ORN")
  
  mesh = as.fdaPDE.SpatialLinesDataFrame(ORN.nt)
  mesh = normalize_mesh(mesh)
  mesh = refine.mesh.1.5D(mesh, delta=0.05)
  
  mesh.2D = create.mesh.2D(nodes=mesh$nodes)
  mesh.2D = refine.mesh.2D(mesh.2D, minimum_angle = 20, maximum_area = 0.05/4)
  
  FEMbasis = create.FEM.basis(mesh)
  FEMbasis.2D = create.FEM.basis(mesh.2D)
  
  nnodes = nrow(mesh$nodes)
  
  res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
             mesh.2D = mesh.2D, FEMbasis.2D = FEMbasis.2D)
  
  
  return(res)
  }else if( network == "estevan"){
    
    data("ERN_OSM_correct")
    
    mesh = as.fdaPDE.SpatialLinesDataFrame(ERN_OSM_cor.nt)
    mesh = normalize_mesh(mesh)
    mesh = refine.mesh.1.5D(mesh, delta=0.05)
    
    mesh.2D = create.mesh.2D(mesh$nodes)
    mesh.2D = refine.mesh.2D(mesh.2D, minimum_angle = 20, maximum_area = 0.0125/4 )
    
    FEMbasis = create.FEM.basis(mesh)
    FEMbasis.2D = create.FEM.basis(mesh.2D)
    
    nnodes = nrow(mesh$nodes)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               mesh.2D = mesh.2D, FEMbasis.2D = FEMbasis.2D)
    
  }
  
  
}

setting.NO.norm <-function(network="estevan"){
  if(network=="estevan"){
  data("ERN_OSM_correct")
  
  mesh = as.fdaPDE.SpatialLinesDataFrame(ERN_OSM_cor.nt)
  mesh = refine.mesh.1.5D(mesh, delta=0.0125/8)
  
  FEMbasis = create.FEM.basis(mesh)
  FEMbasis.2D = create.FEM.basis(mesh.2D)
  
  nnodes = nrow(mesh$nodes)
  
  res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
             mesh.2D = mesh.2D, FEMbasis.2D = FEMbasis.2D, f_ = f.exp)
  }
  else if(network=="ontario"){
    
    data("ORN")
    
    mesh = as.fdaPDE.SpatialLinesDataFrame(ORN.nt)
    mesh = normalize_mesh(mesh)
    mesh = refine.mesh.1.5D(mesh, delta=0.05)
    
    mesh.2D = create.mesh.2D(nodes=mesh$nodes)
    mesh.2D = refine.mesh.2D(mesh.2D, minimum_angle = 20, maximum_area = 0.05/4)
    
    FEMbasis = create.FEM.basis(mesh)
    FEMbasis.2D = create.FEM.basis(mesh.2D)
    
    nnodes = nrow(mesh$nodes)
    
    res = list(mesh = mesh, FEMbasis=FEMbasis, nnodes = nnodes,
               mesh.2D = mesh.2D, FEMbasis.2D = FEMbasis.2D)
  }
  
}
