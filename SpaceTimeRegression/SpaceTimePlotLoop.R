library(viridis)
library(latex2exp)
library(grid)
library(gridExtra)
source("../Auxiliary/R_plot_graph.ggplot2.R")
source("../SpatialRegression-NoCovariate/NoCovariatesCore.R")

SpaceTimeLoop<-function(imgfile,
                        FEMbasis,
                        time_locations,
                        field,
                        mean.field.fdaPDE,
                        palette,
                        line.size=1)
{

  if(!is.null(palette)){
    if(palette=="viridis")
      p=viridis
    else if(palette=="magma")
      p=magma
  }else{
    p=jet.col
  }
  
FEM_time = FEM.time(coeff=mean.field.fdaPDE[,4],
                    time_mesh=time_locations,
                    FEMbasis = FEMbasis,
                    FLAG_PARABOLIC = T) 
time_evals = seq(from=time_locations[1], 
                 to=time_locations[length(time_locations)], 
                 length.out=50)

tmp = eval.FEM.time(FEM.time=FEM_time, locations = FEMbasis$mesh$nodes, time.instants = time_evals)


coeff.evals = matrix(tmp, nrow=FEMbasis$nbasis, ncol=length(time_evals))

FEM_time = FEM.time(coeff=field,
                    time_mesh=time_locations,
                    FEMbasis = FEMbasis,
                    FLAG_PARABOLIC = T)
tmp = eval.FEM.time(FEM.time=FEM_time, locations = FEMbasis$mesh$nodes, time.instants = time_evals)
coeff.ex.evals = matrix(tmp, nrow=FEMbasis$nbasis, ncol=length(time_evals))

col.max = max(coeff.evals, coeff.ex.evals)
col.min = min(coeff.evals, coeff.ex.evals)

estimate.field = list()
true.field = list()
loop.list = list()

pdf(imgfile)
for(i in 1:length(time_evals)){
  title_ = TeX(sprintf("$t = %f", round(time_evals[i],digits = 2)))
  estimate.field[[i]] <- R_plot_graph.ggplot2.2(FEM(coeff.evals[,i], FEMbasis),
                                                line.size = line.size,
                                                color.max = col.max,
                                                color.min = col.min,
                                                title = bquote(hat(f)),
                                                return.ggplot.object = T,
                                                legend.pos = "right",
                                                palette=p,
                                                title.size = 16)
  
  true.field[[i]] <- R_plot_graph.ggplot2.2(FEM(coeff.ex.evals[,i], FEMbasis),
                                            line.size = line.size,
                                            color.max = col.max,
                                            color.min = col.min,
                                            title = bquote(f),
                                            return.ggplot.object = T,
                                            legend.pos = "right",
                                            palette=p,
                                            title.size = 16)
  grid.arrange(true.field[[i]], estimate.field[[i]], 
               ncol=2, nrow=1, top=textGrob(title_, gp=gpar(fontsize=26,font=8)) )
  
}
dev.off()


#begin=0.25
#end=0.75

#res = boxplot_RMSE(RMSE,n_data = n_data, model_=c(T,T,F,F), palette = "viridis", begin=0.35, end=0.65)

#tmp = gray.colors(2, start = begin,end=end)  
#res + scale_fill_grey(start=begin, end=end) + scale_color_manual(values=darken(tmp, amount = 0.25))
}


SpaceTime6Plots<-function(imgfile,
                        FEMbasis,
                        time_locations,
                        field,
                        mean.field.fdaPDE,
                        palette,
                        line.size=1){
  
  
  if(!is.null(palette)){
    if(palette=="viridis")
      p=viridis
    else if(palette=="magma")
      p=magma
  }else{
    p=jet.col
  }
  
  #FEM_time = FEM.time(coeff=mean.field.fdaPDE[,4],
  #                    time_mesh=time_locations,
  #                    FEMbasis = FEMbasis,
  #                    FLAG_PARABOLIC = T) 
  
  time_evals = time_locations
  
  #tmp = eval.FEM.time(FEM.time=FEM_time, locations = FEMbasis$mesh$nodes, time.instants = time_evals)
  
  
  coeff.evals = matrix(mean.field.fdaPDE[,4], 
                       nrow=FEMbasis$nbasis, 
                       ncol=length(time_evals))
  
  #FEM_time = FEM.time(coeff=field,
  #                    time_mesh=time_locations,
  #                    FEMbasis = FEMbasis,
  #                    FLAG_PARABOLIC = T)
  #tmp = eval.FEM.time(FEM.time=FEM_time, locations = FEMbasis$mesh$nodes, time.instants = time_evals)
  coeff.ex.evals = matrix(field, 
                          nrow=FEMbasis$nbasis, 
                          ncol=length(time_evals))
  
  col.max = max(coeff.evals, coeff.ex.evals)
  col.min = min(coeff.evals, coeff.ex.evals)
  
  estimate.field = list()
  true.field = list()
  loop.list = list()
  
  pdf(imgfile)
  for(i in 1:length(time_evals)){
    title_ = TeX(sprintf("$t = %f", round(time_evals[i],digits = 2)))
    estimate.field[[i]] <- R_plot_graph.ggplot2.2(FEM(coeff.evals[,i], FEMbasis),
                                                  line.size = line.size,
                                                  color.max = col.max,
                                                  color.min = col.min,
                                                  #title = bquote(hat(f)),
                                                  return.ggplot.object = T,
                                                  legend.pos = "right",
                                                  palette=p,
                                                  title = " ",
                                                  title.size = 30)
    
    true.field[[i]] <- R_plot_graph.ggplot2.2(FEM(coeff.ex.evals[,i], FEMbasis),
                                              line.size = line.size,
                                              color.max = col.max,
                                              color.min = col.min,
                                              #title = bquote(f),
                                              return.ggplot.object = T,
                                              legend.pos = "right",
                                              palette=p,
                                              title.size = 30,
                                              title = title_)
    
    print(true.field[[i]])
    print(estimate.field[[i]])
    
  }
  dev.off()
  
  
  
  
}

