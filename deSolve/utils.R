neigh1D <- function(edges){
    
    neigh = matrix(0, nrow=nrow(edges), ncol=2)
    for(i in 1:nrow(edges))
        for(j in 1:2)
            if( (i == 1 && j ==2) || ( i == nrow(edges) && j == 1))
                neigh[i,j] = -1
            else if(  (i != nrow(edges) && j==1)){   
                neigh[i,j] = i+1
            }                                             
            else if(  (i != 1 && j==2)){
                neigh[i,j] = i-1
            }

    return(neigh)
}