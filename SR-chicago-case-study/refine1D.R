refine1D<-function(Rnodes, Redges, Rdelta){
  
  edges_old = Redges
  nodes = Rnodes 
  delta = Rdelta
  lengths = matrix(0, nrow = nrow(edges_old), ncol=1)
  num_subs = matrix(0,nrow = nrow(edges_old), ncol=1)
  
  # vector #num_tot_edges. Contiene il valore del vecchio lato a cui appartiene 
  new_to_old = matrix(nrow=0,ncol=1)
  
  num_new_nodes = 0
  num_tot_edges = 0
  for(i in 1: nrow(edges_old))
  {
    x0 = nodes[edges_old[i,1],1]
    y0 = nodes[edges_old[i,1],2]
    x1 = nodes[edges_old[i,2],1]
    y1 = nodes[edges_old[i,2],2]
    
    lengths[i] = sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );
    
    # number of sub_edges
    if( lengths[i] > delta )
      num_subs[i] = ceiling( lengths[i]/delta )
    else 
      num_subs[i] = 1
    
    # number of new nodes
    num_new_nodes = num_new_nodes + (num_subs[i] - 1)
    # number of edges
    num_tot_edges = num_tot_edges + num_subs[i];
  }
  
  edges = matrix(0, nrow= num_tot_edges, ncol=2)
  nodes_new = matrix(0, nrow= num_new_nodes, ncol=2)
  
  total_nodes = nrow(nodes)
  edges_count = 0
  nodes_new_count = 0
  
  eps = 10 * .Machine$double.eps
  
  for( i in 1:nrow(edges_old))
  {
    if( num_subs[i] == 1)
    {   
        edges_count = edges_count + 1
        #Indexes in R starts from 1, in C++ from 0, needed transformations!
        edges[edges_count,1] = edges_old[i,1]
        edges[edges_count,2] = edges_old[i,2]
        new_to_old = rbind(new_to_old, i)
    }
    else # there is at least one new internal node!
      {
        x0 = nodes[edges_old[i,1],1]
        y0 = nodes[edges_old[i,1],2]
        x1 = nodes[edges_old[i,2],1]
        y1 = nodes[edges_old[i,2],2]
        
        cos_ = 0
        sin_ = 0
        
        if( abs(x1-x0) < eps && (y1-y0) > 0  ){
          cos_ =  0.
          sin_ =  1.
        }else if( abs( x1-x0) < eps && (y1-y0) < 0 ){
          cos_ =  0.
          sin_ = -1.
        }else if( abs( y1-y0)<eps && (x1-x0) > 0 ) {
          cos_ =  1.
          sin_ =  0.
        }else if( abs( y1-y0)<eps && (x1-x0) < 0 ){
          cos_ = -1.
          sin_ =  0.
        }else{
          slope = (y1-y0)/(x1-x0)
          if((x1-x0)>0.)
            cos_ = sqrt(1./(1.+slope*slope))
          else 
            cos_ = -sqrt(1./(1.+slope*slope))
          if((y1-y0)>0.)
            sin_ =sqrt(slope*slope/(1.+slope*slope))
          else 
            sin_ = -sqrt(slope*slope/(1.+slope*slope))
          }
        
        delta_ = lengths[i]/num_subs[i]
        #// vector storing the global numbers of the new nodes which belong to old current edges
        #// If N = num_subs[i] then there are N-1 internal nodes (the newest ones)
        nodes_global_num = matrix(0,  nrow = num_subs[i] + 1, ncol=1)
        
        #//Fixing the first and the last nodes
        nodes_global_num[1] = edges_old[i,1]
        nodes_global_num[ num_subs[i] + 1] = edges_old[i,2]
        
        #// Computing new nodes coordinates
        for( n in 2:num_subs[i]) # && n<(num_subs[i] - 1))
        {
          nodes_global_num[n] = (n - 1) + total_nodes
          nodes_new[nodes_new_count + n - 1, 1] = x0 + delta_ * cos_ * (n-1)
          nodes_new[nodes_new_count + n - 1, 2] = y0 + delta_ * sin_ * (n-1)
        }
        
        nodes_new_count =  nodes_new_count + num_subs[i] - 1;
        total_nodes = total_nodes + num_subs[i] - 1
        
        for(e in 1:num_subs[i])
        {
            #Indexes in R starts from 1, in C++ from 0, needed transformations!
            edges_count = edges_count + 1
            edges[edges_count,1] =  nodes_global_num[e]   
            edges[edges_count,2] = nodes_global_num[e+1] 
            new_to_old = rbind(new_to_old, i)
        }
      }
    
    
    
  }
  
  rownames(new_to_old) <- NULL
  result = list(edges=edges, nodes = rbind(nodes,nodes_new), new_to_old = new_to_old)
  return(result)
}

trapezoidal <- function(f,point1, point2){
 delta = norm(point1-point2, "2")
 I = (f[1] + f[2])/2 * delta
 return(I)
}
