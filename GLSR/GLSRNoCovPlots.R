library(gridExtra)
library(ggplot2)
library(latex2exp)
source("C:/Users/Aldo/Documents/SimulationStudies/Auxiliary/R_plot_graph.ggplot2.R")

GLRNoCovPlots<-function(imgfile,
                        field, line.size.field = 0.5,   
                        response,
                        RMSE,legend.pos.RMSE = "right"){
  
  
  true.spatial.field<- R_plot_graph.ggplot2.2(FEM(field, FEMbasis),
                                              line.size = line.size.field,
                                              title = "True Spatial Field",
                                              ratio=1, 
                                              return.ggplot.object = T,
                                              legend.pos = "right") 
  

  sample_ = sample(1:nnodes, 350)
  points_ = mesh$nodes[sample_,]
  mu_ =  response[sample_]
  response.example <- R_plot_mesh.ggplot(mesh = mesh,points_ = points_, mu = mu_, 
                                         title = "Sample", num.col = (diff(range(mu_))+1) )
  rmse <- boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F), 
                       names_ = c("fdaPDE","GWR","",""),
                       legend.pos = legend.pos.RMSE)
  
  pdf(imgfile)
  print(true.spatial.field)
  print(response.example)
  print(rmse)
  dev.off()
}
