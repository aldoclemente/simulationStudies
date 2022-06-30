#######################
##### LHP utils #######
#######################


set_Kfold_data <-function(SpPointsDataFrame, seed = 27182, K = 10){
set.seed(seed) # 31415 # 1234

SpPointsDataFrame = SpPointsDataFrame[sample(1:nrow(SpPointsDataFrame)), ]
  
listData = list()

num_data = round(nrow(SpPointsDataFrame)/K)
for(i in 1:(K-1)){
    listData[[i]] = SpPointsDataFrame[(1 + num_data*(i-1)):(num_data*i),]
}
listData[[K]] = SpPointsDataFrame[(num_data*(K-1) + 1):nrow(SpPointsDataFrame), ]

return(listData)
}
  
get_Kfold_data <- function(SpLinesDataSet_List, iter, K = 10){
  
  tmp_data = data.frame()
  tmp_coords = matrix(0,nrow=0,ncol=2)
  for(i in 1:K){
    if( i == iter){
      test_data = SpLinesDataSet_List[[i]]
    }else{
      tmp_data = rbind(tmp_data, SpLinesDataSet_List[[i]]@data)
      tmp_coords = rbind(tmp_coords, SpLinesDataSet_List[[i]]@coords)
    }
  }
  
  train_data = SpatialPointsDataFrame(coords = tmp_coords,
                                      data = tmp_data)
  
  ret_list = list(train_data = train_data, test_data = test_data)
  return(ret_list)
}

boxplot_RMSE <-function(rmse.SR.PDE, rmse.GWR,
                        title.size=20,
                        begin=0.95, #color
                        end=0.25,   #color
                        width =0.75,
                        title="RMSE")
{
  
  model = rep(c("GWR", "SR-PDE"), each=length(rmse.SR.PDE))
  RMSE =  c(as.vector(rmse.GWR), as.vector(rmse.SR.PDE))
  dataFrame = data.frame(RMSE=RMSE, model = model)
  
  MyTheme <- theme(
    axis.text = element_text(size=title.size-10),
    axis.title = element_text(size=title.size),
    title = element_text(size=title.size),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size=title.size-6),
    legend.key.size = unit(1,"cm"),
    legend.key.height = unit(1,"cm"),
    legend.title = element_blank(),
    legend.background = element_rect(fill="white", color="black",
                                     size=c(1,0.5))
  )
  
  border_col = darken(viridis(2, begin=begin,end=end), amount=0.25)
  
  p<-ggplot(dataFrame)+
    geom_boxplot(aes(x=model,
                     y=RMSE, group=model,
                     fill=model,
                     color=model), width=width)+
    scale_x_discrete(limits=c("GWR", "SR-PDE"))+
    labs(x="", y="",
         title=title,)+
    scale_fill_viridis(begin = begin,
                       end = end,
                       option = "viridis", discrete=T) +
    scale_color_manual(values=border_col) +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    MyTheme + 
    theme(#plot.title=element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none")
  return(p)  
}
