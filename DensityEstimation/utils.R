###########################
### DE case study utils ###
###########################

set_Kfold_data <-function(data, seed = 27182, K = 10){
  set.seed(seed) 
  
  data = data[sample(1:nrow(data)), ]
  
  listData = list()
  
  num_data = round(nrow(data)/K)
  for(i in 1:(K-1)){
    listData[[i]] = data[(1 + num_data*(i-1)):(num_data*i),]
  }
  listData[[K]] = data[(num_data*(K-1) + 1):nrow(data), ]
  
  return(listData)
}

get_Kfold_data <- function(dataList, iter, K = 10){
  
  train_data = list()
  for(i in 1:K){
    if( i == iter){
      test_data = dataList[[i]]
    }else{
      train_data = rbind(train_data, dataList[[i]])
      
    }
  }
  
  ret_list = list(train_data = train_data, test_data = test_data)
  return(ret_list)
}


#' Compute the CV error. 
#' @param FEM, fdaPDE function
#' @param R0, Mass matrix, NB. R0 = CPP_get.FEM.Mass.Matrix(FEMbasis)
#' @param data.k, k-th data fold
#'
cv_error <- function(FEM, R0, data.k){
  f = FEM$coeff
  f.eval.k = eval.FEM(FEM, locations = cbind(data.k$x, data.k$y))
  
  return( as.numeric( t(f^2) %*% R0 %*% f^2  - 2*mean(f.eval.k)) )
}
