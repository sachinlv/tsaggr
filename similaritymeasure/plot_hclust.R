library('dendroextras')

dist.mat.filename <<- 'distance_matrix_'
plot.file.name <<- 'hclust_'
error.file <<- '_shortterm_hclust_aggr.csv'

site.type.vec <- c('random', 'specific')
algo.vec <- c('neuralnet', 'gbm', 'brnn', 'mars')
simi.meas.vec <- c(
  'euclidean',
  'fourier',
  'pca',
  'lm',
  'weuclid')
cut.count <- 10

generate.dendro.plots <- function(meas, algo){
  error.file <- paste(error.file.path,
                      algo,"_shortterm_windspeed_hclust/",
                      meas,'/',
                      algo,'_shortterm_hclust_aggr.csv',
                      sep="")
  error <- data.frame(read.csv(error.file))

  plot.file <- paste(plot.file.path, plot.file.name, meas,'.pdf', sep="")
  dist.file <- paste(dist.mat.file.path, dist.mat.filename, meas,'.csv',sep="")
  dist.mat <- data.frame(read.csv(dist.file))
  dist.mat <- dist.mat[,2:11]
  colnames(dist.mat) <- paste('site',seq(1,10),sep="")

  par(mfrow=c(5,2),
      cex=0.4, cex.axis=0.5, cex.lab=1, cex.main=1.5,
      font.lab=5, lab=c(10,10,10),
      oma=c(2,2,2,2), mar=c(5.1,4.1,4.1,2.1))
  #par(mar=c(4.0, 4.0, 4.0, 4.0))
  plot.window(xlim=c(0,400),ylim=c(0,1),asp=1)

  for(k in error$cutsize){
    h = hclust(as.dist(dist.mat))
    p = color_clusters(h, k,groupLabels=TRUE)
    plot(p, main=paste("Forecast Error=", error$rmse[k]))
    dev.copy2pdf(file=plot.file, height=10, width=8)
  }
}

generate.plots <- function(){
  for(site.type in site.type.vec){
    dist.mat.file.path <<- paste('/home/freak/Dropbox/results2/resultsreview/',site.type,'_sites/distance/',sep="")
    plot.file.path <<- paste('/home/freak/Programming/Thesis/results/plots/hclust/',site.type,'_sites/',sep="")
    error.file.path <<- paste('/home/freak/Dropbox/results2/clustering/history10/',site.type,'_sites/',sep="")
    for(algo in algo.vec){
      for(meas in simi.meas.vec){
        generate.dendro.plots(meas, algo)
        break
      }
      break
    }
    break
  }
}

generate.plots()
