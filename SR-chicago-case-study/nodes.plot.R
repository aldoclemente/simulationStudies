
if( !dir.exists("imgs/") )
    dir.create("imgs")
pdf("nodes.pdf")
for(i in 1:nrow(mesh$nodes)){
   # png(paste("imgs/node-",i,sep=""))
    print(plot(mesh, main=paste("node ", i, sep="")))
    print(points(chicago$data$x,chicago$data$y, col="red", pch=16, cex =1.25))
    print(points(mesh$nodes[i,1], mesh$nodes[i,2], col = "green4", pch=16, cex=1.25))
    #dev.off()
}
dev.off()

sources_ = c(47, 91, 67, 124, 151, 169, 195, 251, 275, 302)
plot(mesh)
points(chicago$data$x, chicago$data$y, pch=16, col="red3", cex=2.5)
points(mesh$nodes[sources_,1], mesh$nodes[sources_,2], pch=16, col="blue", cex=2.5)

set_region <- function(sources_, LN = chicago){

    lpp_sources <- spatstat.linnet::lpp(X= ppp(LN$domain$vertices$x[sources_], y=LN$domain$vertices$y[sources_], window= LN$domain$window), 
                                        L= LN$domain)

    nedges = LN$domain$lines$n
    mid_points = matrix(0, nrow=nedges, ncol=2)
    for(e in 1:nedges){  
        mid_points[e,1] =  ( LN$domain$vertices$x[ LN$domain$from[e]] + LN$domain$vertices$x[ LN$domain$to[e]] )/2
        mid_points[e,2] =  ( LN$domain$vertices$y[ LN$domain$from[e]] + LN$domain$vertices$y[ LN$domain$to[e]] )/2
    }

    lpp_midpoints <- spatstat.linnet::lpp(X = ppp(x=mid_points[,1], y= mid_points[,2], window= LN$domain$window),
                                          L = LN$domain)

    network_dist = spatstat.linnet::crossdist.lpp(lpp_midpoints, lpp_sources)
    edges_to_region = matrix(0, nrow = nedges, ncol=1)

    for(e in 1:nedges) edges_to_region[e] = which( network_dist[e,] == min(network_dist[e,]) )

    return(edges_to_region)
}