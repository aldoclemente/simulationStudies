library(gridExtra)
library(ggplot2)
library(latex2exp)
source("../Auxiliary/R_plot_graph.ggplot2.R")

SpaceTimePlots <- function(imgfile, 
                           time_locations,
                           den, # time_locations: from=0 to= pi/den
                           field, 
                           observations, 
                           W,
                           betas,
                           RMSE,legend.pos.RMSE = "right"){

  min.col = min(field)
  max.col = max(field)
  coef.ex = matrix(field, 
                   nrow=length(field)/length(time_locations),
                   ncol=length(time_locations))
  
  ggplot_list = list()
  col.min = min(coef.ex)
  col.max = max(coef.ex)  
  ggplot_list = list()
  n_istants = length(time_locations)
  for( i in 1:n_time_locs){
    if(i == 1)
      title_ = TeX(sprintf("$t = 0"))
    else if(i==2)
      title_ = TeX(sprintf("$t = \\pi / %d$",den*(n_istants-1)))
    else if(i != n_time_locs)
      title_ = TeX(sprintf("$t =  %d \\pi / %d$",i-1,den*(n_istants-1)))
    else{
      if(den==1){
        title_ = TeX(sprintf("$t =  \\pi$"))
      }else{
        title_ = TeX(sprintf("$t =  \\pi / %d$", den))  
      }
    }
    
    #print(
    ggplot_list[[i]] <- R_plot_graph.ggplot2.2(FEM(coef.ex[,i], FEMbasis),
                                               line.size = 1,
                                               color.max = col.max,
                                               color.min = col.min,
                                               title = title_,
                                               return.ggplot.object = T,
                                               legend.pos = "right") 
  }  
  
  
pdf(imgfile, height=9,width = 9)
grid.arrange(grobs = ggplot_list, 
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)))




grid.arrange(grobs = ggplot_list)

tmp.sample_ = sample(x=(1:nnodes), size=n_data[1])
sample_ = vector(mode="integer")
for(k in 1:length(time_locations)){
  sample_ = rbind(sample_, tmp.sample_ + (k-1)*nnodes)
}
sample_ = as.vector(sample_)
points_ = mesh$nodes[tmp.sample_,]

W.1 = W[sample_,1]
W.1 = matrix(W.1, nrow = n_data[1], ncol=length(time_locations))
W.2 = W[sample_,2]
W.2 = matrix(W.2, nrow = n_data[1], ncol=length(time_locations))

firstCov <- R_plot_mesh.ggplot(mesh = mesh,points_ = points_, mu = W.1[,1], title = "First Covariate")
secondCov <- R_plot_mesh.ggplot(mesh = mesh,points_ = points_, mu = W.2[,2], title = "Second Covariate")

RMSE <- boxplot_RMSE(RMSE, n_data, model_ = c(T,T,F,F), 
             names_ = c("fdaPDE","GWR","",""),
             legend.pos = legend.pos.RMSE)

print(firstCov)
print(secondCov)
print(RMSE)
dev.off()

}
