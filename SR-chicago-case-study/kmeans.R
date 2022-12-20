source("refine1D.R")
source("../utils.R")

#
#centroids_ = cbind(LN$domain$vertices$x[sources_], LN$domain$vertices$y[sources_])

ncluster = 10
ndata = length(LN$data$x)
sample_ = sample(1:ndata, ncluster)
centroids_ = cbind(LN$data$x[sample_], LN$data$y[sample_])

kmeans_graphs<-function(centroids_, mesh, LN = chicago, tol = 1e-6, max_iter = 1e3){
    J_old <- 1e10
    J <- 0
    err <- 1    
    ncluster <- nrow(centroids_)
    ndata = length(LN$data$x)
    niter <- 0
    while( niter < max_iter && err > tol){
        lpp_centroid <- spatstat.linnet::lpp(X= ppp(centroids_[,1], y=centroids_[,2], 
                                        window= LN$domain$window), 
                                        L= LN$domain)
        lpp_data <- spatstat.linnet::lpp(X= ppp(x=LN$data$x, y= LN$data$y, 
                                     window= LN$domain$window),
                                     L= LN$domain)

        network_dist = spatstat.linnet::crossdist.lpp(lpp_data, lpp_centroid)

        data_to_centroid = matrix(0, nrow = ndata, ncol=1)
        for(e in 1:ndata) 
            data_to_centroid[e] = which(network_dist[e,] == min(network_dist[e,]) )

        for(i in 1:ncluster){
          mask_ <- which(data_to_centroid == i)
          J = J + sum(network_dist[mask_,i])
          
          centroids_[i,1] = sum(lpp_data$data$x[mask_])/ncluster
          centroids_[i,2] = sum(lpp_data$data$y[mask_])/ncluster
        }
        err = abs(J-J_old)/J_old
        
        centroids_ = fdaPDE::projection.points.1.5D(mesh= mesh, locations= centroids_)
        J_old = J
        niter = niter + 1

    }

    ret = list(niter = niter, centroids = centroids_, error= err, J= J)
    return(ret)
}

data_ = cbind(LN$data$x, LN$data$y)
result_ <- kmeans(x=data_, centers=ncluster, iter.max = 100)
x11()
plot(mesh, pch=".")
points(result_$centers, col="blue", pch = 16, cex=2)
points(data_, col="red", pch=16, cex=2)

centroids_ = fdaPDE::projection.points.1.5D(mesh, locations= result_$centers)
