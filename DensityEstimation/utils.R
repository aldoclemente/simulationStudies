###########################
### DE case study utils ###
###########################

library(ggplot2)
library(viridis)
library(colorspace)
library(grid)
library(gridExtra)
source("../Auxiliary/R_plot_graph.ggplot2.R")

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

boxplot_CV_error <-function(CV_errors,
                            methods,
                            methods.names,
                        title.size=20,
                        begin=0.95, #color
                        end=0.25,   #color
                        width =0.75,
                        title="CV error")
{
  
  METHODS = rep(methods.names[methods], each=nrow(CV_errors))
  RMSE =  as.vector(CV_errors)
  dataFrame = data.frame(RMSE=RMSE, METHODS = METHODS)
  
  MyTheme <- theme(
    axis.text = element_text(size=title.size-5),
    axis.title = element_text(size=title.size),
    title = element_text(size=title.size),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size=title.size-5),
    legend.key.size = unit(1,"cm"),
    legend.key.height = unit(1,"cm"),
    legend.title = element_blank(),
    legend.background = element_rect(fill="white", color="black",
                                     size=c(1,0.5))
  )
  
  border_col = darken(viridis(ncol(CV_errors), begin=end,end=begin), amount=0.25)
  fill_col = viridis(ncol(CV_errors), begin=end, end=begin)
  
  BORDER = c()
  FILL = c()
  for(i in 1:length(methods)){
    if(methods[i]){ 
      FILL = append(FILL, fill_col[i])
      BORDER = append(BORDER, border_col[i])
    }
  }
  
  dataFrame$METHODS = factor(dataFrame$METHODS, 
                             levels=methods.names) 
  
  
  p<-ggplot(dataFrame)+
    geom_boxplot(aes(x=METHODS,
                     y=RMSE, group=METHODS,
                     fill=METHODS,
                     color=METHODS), width=width)+
    scale_x_discrete(limits=methods.names[methods])+
    labs(x="", y="",
         title=title)+
    scale_fill_viridis(begin = end,
                       end = begin,
                       option = "viridis", discrete=T) + #ok
    scale_color_manual(values=border_col) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    MyTheme + 
    theme(#plot.title=element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none")
  return(p)  
}

# estimates, lista di oggetti FEM ( fem = FEM(coeff, FEMbasis))
# palette, color palette (es. jet.col, viridis, magma, inferno, plasma)
# line.size, dimensione segmenti network
# titles, vettore di titoli dei plot 
plot_estimates_CASE_STUDY <-function(estimates, # list of estimates
                          palette = jet.col, #palette = NULL (ggplot default colors), "viridis", "magma"
                          line.size=1,
                          titles = c(rep("",times=length(estimates))) )
  {
  
  max.col = -1e8
  min.col = 1e8
  for(i in 1:length(estimates)){
    max.col = max(estimates[[i]]$coeff, max.col)
    min.col = min(estimates[[i]]$coeff, min.col)
  }
  
  estimates.plot = list()
  
  for(i in 1:length(estimates)){
    
    estimates.plot[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                                 line.size = line.size,
                                                 color.min = min.col,
                                                 color.max = max.col,
                                                 title = titles[[i]],
                                                 palette=palette,
                                                 legend.pos = "right")
    
  }
  
  return(estimates.plot)
  
}

boxplot_RMSE <- function(RMSE,
                         methods,
                         methods.names,
                         nsim,
                         title.size=20,
                         begin=0.95, #color
                         end=0.25,   #color
                         width =0.75,
                         n = c(50, 100, 150, 250),
                         title="RMSE")
{
  RMSE = RMSE[, methods]
  
  METHODS = rep(methods.names[methods], each = nrow(RMSE)) 
  N = as.character(rep(rep(n, each = nsim), times = sum(methods)))
  
  RMSE =  as.vector(RMSE)
  dataFrame = data.frame(RMSE=RMSE, METHODS = METHODS, N = N)
  
  MyTheme <- theme(
    axis.text = element_text(size=title.size-5),
    axis.title = element_text(size=title.size),
    title = element_text(size=title.size),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size=title.size-5),
    legend.key.size = unit(1,"cm"),
    legend.key.height = unit(1,"cm"),
    legend.title = element_blank(),
    legend.background = element_rect(fill="white", color="black",
                                     size=c(1,0.5))
  )
  
  border_col = darken(viridis(length(methods), begin=end,end=begin), amount=0.25)
  fill_col = viridis(length(methods), begin=end, end=begin)
  
  BORDER = c()
  FILL = c()
  for(i in 1:length(methods)){
    if(methods[i]){ 
      FILL = append(FILL, fill_col[i])
      BORDER = append(BORDER, border_col[i])
      }
  }
  
  dataFrame$METHODS = factor(dataFrame$METHODS, 
                              levels=methods.names) 
  
  p<-ggplot(dataFrame)+
    geom_boxplot(aes(x=N,
                     y=RMSE, group=interaction(METHODS,N),
                     fill=METHODS, color = METHODS))+
    scale_x_discrete(limits=as.character(n))+
    labs(x="", y="",
         title=title)+
    scale_fill_manual(values = FILL) +
    scale_color_manual(values= BORDER) + 
    MyTheme + 
    theme(#plot.title=element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = c(0.85,0.85) )
  return(p)  
}

plot_density<-function(true.density,
                       max.col,
                       min.col,
                       main =  "Density", 
                       palette = NULL, #palette = NULL (ggplot default colors), "viridis", "magma"
                       line.size=1)
{
  
  if(!is.null(palette)){
    if(palette=="viridis")
      p=viridis
    else if(palette=="magma")
      p=magma
    else if(palette=="inferno")
      p=inferno
  }else{
    p=jet.col
  }
  
  max.col = max(max.col, true.density$coeff)
  min.col = min(min.col, true.density$coeff)
  
  density.plot <- R_plot_graph.ggplot2.2( true.density, # FEM object
                                          line.size = line.size,
                                          color.min = min.col,
                                          color.max = max.col,
                                          title = main,
                                          palette=p,
                                          legend.pos = "right")
  
  ret = list(max.col = max.col, min.col = min.col, density.plot = density.plot)  
  return(ret)
}


plot_estimates <-function(estimates, # list of estimates
                          true.density,
                          methods,
                          methods.names =  c(rep("",times=length(estimates))), 
                          true.density.main = "Density",
                          palette = NULL, #palette = NULL (ggplot default colors), "viridis", "magma"
                          line.size=1)
{
  
  if(!is.null(palette)){
    if(palette=="viridis")
      p=viridis
    else if(palette=="magma")
      p=magma
    else if(palette=="inferno")
      p=inferno
  }else{
    p=jet.col
  }
  
  max.col = -1e8
  min.col =  1e8
  
  for(i in 1:length(methods.names)){
    max.col = max(max.col, estimates[[i]]$coeff)
    min.col = min(min.col, estimates[[i]]$coeff)
  }
  
  ret.list = plot_density(true.density, max.col, min.col, main =  true.density.main, 
                         palette = palette, 
                         line.size=line.size)
  
  max.col = ret.list$max.col
  min.col = ret.list$min.col
  estimates.plot = list()
  
  for(i in 1:length(estimates)){
    if( methods[i]){
    estimates.plot[[i]] = R_plot_graph.ggplot2.2( estimates[[i]], # FEM object
                                                  line.size = line.size,
                                                  color.min = min.col,
                                                  color.max = max.col,
                                                  title = methods.names[[i]],
                                                  palette=p,
                                                  legend.pos = "right")
    }
  }
  
  ret = list(estimates.plot = estimates.plot, density.plot = ret.list$density.plot) 
  return(ret)
  
}