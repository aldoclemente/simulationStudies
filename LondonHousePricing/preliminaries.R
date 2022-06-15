data(LNNT)
data(LNHP)

# STEP 1
# coordinates data
ptsxy<-coordinates(LN.prop)[,1:2]

# integrate an individual point data set into a network under different rules. 
# In particular, we map each point to the nearest node in the network/graph.
# res[[1]] -> A “nodelist” object (shp2graph)
# res[[2]] -> An “edgelist” object
# res[[3]] -> A vector of the cooresponding node ID for each point in “pointsxy”

res<-points2network(ntdata=LN.nt,pointsxy=ptsxy, approach=1)

# storing STEP 1 data
data.list = list(res=res, ptsxy=ptsxy) 
save(data.list, file="LHP_1.RData")

# STEP 2 - igraph object
# igraph linear network

# rt.NEL[[1]] TRUE if the output is under a “Detailed” mode,and “edgelist” will have a different structure;
# rt.NEL[[2]] A “nodelist” object
# rt.NEL[[3]] An “edgelist” object
# rt.NEL[[4]] A vector (of length equal to the number of edges) of edge lengths if “ELComputed” is TRUE;
# rt.NEL[[5]] A data frame of edge attributes, [EdgeID,...(items extracted from the “SpatialLinesDataFrame” object )...]
# rt.NEL[[6]] A vector contains X-coordinates of all the nodes
# rt.NEL[[7]] A vector contains Y-coordinates of all the nodes
rt.NEL = readshpnw(LN.nt, 
                   ELComputed = TRUE, 
                   Detailed = TRUE,
                   ea.prop = names(LN.nt))

igr = nel2igraph(nodelist =  rt.NEL[[2]], edgelist = rt.NEL[[3]], weight    = rt.NEL[[4]])

# nodes and edges list
data.frame.igr = as_data_frame(igr, what="both")

save(rt.NEL,igr, file="LHP_2.RData")


# STEP 3 - NO DUPLICATES
nloc = nrow(ptsxy)
tmp = vector(mode="integer",length=nloc)
for(i in 1:nloc)
  tmp[i] = data.list$res[[3]][[i]]
  
duplicated_idx = which(duplicated(tmp))
ptsxy.no.duplicated = ptsxy[-duplicated_idx,]
nodes.idx = tmp[-duplicated_idx]

data_ = data.frame(PURCHASE   = LN.prop$PURCHASE[-duplicated_idx],
                     FLOORSZ    = LN.prop$FLOORSZ[-duplicated_idx],
                     TYPEDETCH  = LN.prop$TYPEDETCH[-duplicated_idx], 
                     TPSEMIDTCH = LN.prop$TPSEMIDTCH[-duplicated_idx], 
                     TYPETRRD   = LN.prop$TYPETRRD[-duplicated_idx],
                     TYPEBNGLW  = LN.prop$TYPEBNGLW[-duplicated_idx],
                     TYPEFLAT   = LN.prop$TYPEFLAT[-duplicated_idx],
                     BLDPWW1    = LN.prop$BLDPWW1[-duplicated_idx],
                     BLDPOSTW   = LN.prop$BLDPOSTW[-duplicated_idx],
                     BLD60S     = LN.prop$BLD60S[-duplicated_idx],
                     BLD70S     = LN.prop$BLD70S[-duplicated_idx],
                     BLD80S     = LN.prop$BLD80S[-duplicated_idx],
                     BLD90S     = LN.prop$BLD90S[-duplicated_idx],
                     BATH2      = LN.prop$BATH2[-duplicated_idx],
                     BED2       = LN.prop$BEDS2[-duplicated_idx],
                     GARAGE1    = LN.prop$GARAGE1[-duplicated_idx],
                     CENTHEAT   = LN.prop$CENTHEAT[-duplicated_idx],
                     UNEMPLOY   = LN.prop$UNEMPLOY[-duplicated_idx],
                     PROF       = LN.prop$PROF[-duplicated_idx],
                     BLDINTW    = LN.prop$BLDINTW[-duplicated_idx],
                     X          = data.frame.igr$vertices[nodes.idx,1],
                     Y          = data.frame.igr$vertices[nodes.idx,2])
  
  coords_ = cbind(data_$X, data_$Y)
  
  LN.data = SpatialPointsDataFrame(coords = coords_,
                                   data = data_)
  
  save(ptsxy.no.duplicated,
       nodes.idx,
       LN.data, file="LHP_3.RData")

# STEP 4 - Network Distance matrix (ND)
  
  ND = distances(igr, v = nodes.idx, to = nodes.idx)
  
  save(ND, file = "LHP_4.RData")

  
# CHECK
  for(i in duplicated_idx){
    print(paste(""))
    
  }
  
  res.check<-points2network(ntdata=LN.nt,pointsxy=ptsxy.no.duplicated, approach=1)
