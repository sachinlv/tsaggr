require(RMySQL)
require(ff)
require(RSNNS)
require(TSdist)
require(forecast)
require(timeSeries)
require(TSclust)
require(corpcor)
require(MADAM)

sites.count <- 10
data.len <-7200 #52560
powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
threshold.rmse <<- 0.3
threshold.dist <<- 0.3
drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")
start.date <- '20061112'
end.date <- '20070101'
table.ip.type <- "random"#"specific"
preprocess <- 'raw'#c('raw', 'linear', 'expsmooth', 'wma')
err.type.vec <- c('rmse','mse','sd','cor')
folder.ip.type.vec <- c('aggr','mean','wmean')
#cfolder.ip.type.vec <-('exp_aggr', 'exp_mean','exp_wmean')


results.file.path <-'/home/freak/Programming/Thesis/results/history50/random_sites5/'
plot.file.path.all <- '/home/freak/Programming/Thesis/plots/history50/random_sites5/'

if(table.ip.type == "random"){
  t <- c('onshore_SITE_07726',
         'onshore_SITE_07791',
         'onshore_SITE_02148',
         'onshore_SITE_02797',
         'onshore_SITE_01351',
         'onshore_SITE_03986',
         'onshore_SITE_00324',
         'onshore_SITE_00342',
         'onshore_SITE_01562',
         'onshore_SITE_03069')
  lat=c(41.17, 41, 39.91, 40.13, 36.35, 36.29, 37.84, 37.55, 42.68, 42.71)
  long=c(-74.5, -74.55, -78.89, -78.71, -81.65, -81.72, -100.07, -99.35, -73.29, -73.04)
  tables <<- data.frame(cbind(numeric(0),t))
  coordinates <<- data.frame(long,lat)
}else{
  t <- c("onshore_SITE_00538",
         "onshore_SITE_00366",
         "onshore_SITE_00623",
         "onshore_SITE_00418",
         "onshore_SITE_00627",
         "onshore_SITE_00532",
         "onshore_SITE_00499",
         "onshore_SITE_00571",
         "onshore_SITE_03247",
         "onshore_SITE_00622")
  lat=c(34.6,34.96,35.01,35.09,35.2,35.05,35.24,35.32,33.86,34.44)
  long=c(-102.54,-102.72,-102.32,-102.62,-102.52,-103,-102.76,-102.52,-102.79,-101.7)
  tables <<- data.frame(cbind(numeric(0),t))
  coordinates <<- data.frame(long,lat)
}


algo.vec <<- c('neuralnet','brnn','gbm','mars')
simi.meas.vec <<- c(
  'euclidean',
  'minkowski',
  'manhattan',
  'fourier',
  'correlation',
  'pca',
  #'lcss',
  #'erp',
  #'dtw',
  #'edr',
  'coord',
  'lm',
  'weuclid',
  'maj',
  'shrinkage',
  'pdist'
  #'mahalanobis'
  #'fish'
)

setvals <- function(algorithm, measure, type, folder, err){
  algo <<- algorithm
  sim.meas <<- measure
  ip.type <<- type#'statistical'
  folder.ip.type <<- folder
  err.type <<- err
  combination.result <<- c()

  plot.file.path <- paste(plot.file.path.all,
                          folder.ip.type,'/',
                          err.type,'/',sep="")
  combination.file.path <- paste(plot.file.path.all,
                                 folder.ip.type,'/',
                                 err.type,'/',
                                 'combination/',sep="")

  file.path.all <<- paste(results.file.path,
                          algo,'_shortterm_',
                          folder.ip.type,'/',sep="")
  plot.file <<- paste(plot.file.path,
                      algo,'_',
                      ip.type,'_',
                      sim.meas,'.pdf',sep="")
  result.file <<- paste(combination.file.path,
                        algo,'_',
                        ip.type,'_',
                        sim.meas,'.csv',sep="")

  file.name.all <<- paste(algo,'_shortterm_aggr_combi',sep="")
}


load.data <- function(){
  for(indx in seq(1,sites.count)){
    tab <- as.character(tables[indx,1])
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= ",start.date," && mesdt < ",end.date,") LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- switch(preprocess,
                              "raw"={
                                as.double(data06[,1])
                              },
                              "linear"={
                                d <- as.double(data06[,1])
                                f <- filter(d,rep(1/mean(d),mean(d)))
                                f[is.na(f)] <- 0
                              },
                              "expsmooth"={
                                d <- as.double(data06[,1])
                                f <- tbats(d,
                                           use.trend=FALSE,
                                           use.damped.trend=FALSE,
                                           use.box.cox=FALSE,
                                           use.arma.errors=FALSE,
                                           use.parallel=TRUE,
                                           num.cores=4 )
                                f$fitted.values
                              },
                              "wma"={
                                c()
                              }
    )
    #powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]),type="0_1")
  }
}

load.err.data <- function(file.no){
  plotdata <<- c()

  file <- paste(file.name.all,file.no,'.csv', sep="")
  err.tbl <- read.csv(file,sep = ',')
  seq <- err.tbl$AggrNo.Seq
  err <- switch(err.type,
                'rmse'={
                  err.tbl$rmse
                },
                'mse'={
                  err.tbl$mse
                },
                'sd'={
                  err.tbl$sd
                },
                'cor'={
                  err.tbl$cor
                }
  )

  seq <- gsub("10","0",seq)#Replace 10 with 0 in vector
  plotdata <<- data.frame(seq,err)
}


gen.dist.mat <- function(measure){
  dist.mat <<- matrix(0, nrow=sites.count, ncol=sites.count)
  for(i in seq(1,sites.count)){
    for(j in seq(1, sites.count)){
      data.mat <- data.frame(c1=powdata[,i], c2=powdata[,j])

      dist.mat[i,j] <<- switch(measure,
                               "pca"={
                                 f <- formula("~ c1 + c2")
                                 pc <- princomp(f,data=data.mat, cor=FALSE)
                                 score <- data.frame(pc$scores)
                                 euclideanDistance(score$Comp.1,score$Comp.2)
                               },
                               "euclidean"={
                                 euclideanDistance(data.mat$c1,data.mat$c2)
                               },
                               "correlation"={
                                 correlationDistance(data.mat$c1,data.mat$c2)
                               },
                               "manhattan"={
                                 manhattanDistance(data.mat$c1,data.mat$c2)
                               },
                               "minkowski"={
                                 minkowskiDistance(data.mat$c1,data.mat$c2,2)
                               },
                               "fourier"={
                                 fourierDistance(data.mat$c1,data.mat$c2,n=(floor(length(data.mat$c1)/2)+1))
                               },
                               "mean"={
                                 abs(mean(data.mat$c1-data.mat$c2))
                               },
                               "dtw"={
                                 dtwDistance(data.mat$c1, data.mat$c2)
                               },
                               "edr"={
                                 edrDistance(data.mat$c1, data.mat$c2, epsilon=0.1)#, sigma)
                               },
                               "erp"={
                                 erpDistance(data.mat$c1, data.mat$c2, g=0)
                               },
                               "lcss"={
                                 lcssDistance(data.mat$c1, data.mat$c2, epsilon=0.1)
                               },
                               "coord"={
                                 y1 <- coordinates$long[i]
                                 x1 <- coordinates$lat[i]
                                 y2 <- coordinates$long[j]
                                 x2 <- coordinates$lat[j]
                                 sqrt(((y2-y1)^2) + ((x2-x1)^2))
                               },
                               "lm"={
                                 x <- data.frame(x=data.mat$c1)
                                 y <- data.frame(y=data.mat$c2)
                                 lx <- lm('x ~.', x)
                                 ly <- lm('y ~.', y)
                                 euclideanDistance(lx$fitted.values,ly$fitted.values)
                               },
                               "weuclid"={#distance from patel paper. weighted euclidean distance
                                 v1 <- var(data.mat$c1)
                                 v2 <- var(data.mat$c2)
                                 d <- ((data.mat$c1-data.mat$c2)^2)/(v1+v2)
                                 sum(d)
                               },
                               "maj"={
                                 d = diss.AR.MAH(data.mat$c1, data.mat$c2, dependence=FALSE, permissive=TRUE)
                                 d$statistic #check what is p-value in this
                               },
                               "shrinkage"={
                                 lambda = 0.3
                                 c = cor.shrink(data.mat,lambda)
                                 c[1,2]
                               },
                               'pdist'={
                                 #This cannot be used because- its with predicted values
                                 #and its the distance algo for hclust repeatedly
                                 ni <- 1#No.of series in a cluster
                                 nj <- 1
                                 n <- ((ni * nj)/(ni+nj))
                                 vi <- var(data.mat$c1)/(ni^2)
                                 vj <- var(data.mat$c2)/(nj^2)
                                 mi <- mean(data.mat$c1)
                                 mj <- mean(data.mat$c2)
                                 abs(n *(vi + vj - (mi-mj)^2))
                               },
                               "mahalanobis"={
                                 if(i == j){
                                   d <- 0
                                 }else{
                                   d <- mahalanobis(data.mat,
                                                    colMeans(data.mat),
                                                    cov(data.mat))
                                 }
                                 sqrt(sum(d))
                               }#,
                               #"fish"={},

      )
    }
  }
}

get.distance.mat <- function(aggr.seq){
  len <- length(aggr.seq)
  mat <- matrix(0,nrow=len, ncol=len)
  print(aggr.seq)

  for(i in seq(1,len)){
    for(j in seq(1, len)){
      i1 <- aggr.seq[i]
      j1 <- aggr.seq[j]
      mat[i,j] <- dist.mat[i1,j1]
    }
  }
  return(mat)
}

get.err.dist.data <- function(){
  avg.dist <- c()
  for(i in seq(1,length(plotdata$seq))){
    seq <- unlist(strsplit(as.character(plotdata$seq[i]),''))
    seq <- as.integer(seq[2:length(seq)]) #remove S from sequence
    seq <- replace(seq,seq==0,10)#replace 0 back to 10

    if(length(seq)>1){
      dist.mat <- get.distance.mat(seq)
      d <- as.vector(as.dist(dist.mat))
      dist <- mean(d)
      avg.dist <- c(avg.dist,dist)
    }else{
      avg.dist <- 0
    }
  }

  #combination result set
  for(i in seq(1,length(avg.dist))){
    if(avg.dist[i] <= threshold.dist && plotdata$err[i] <= threshold.rmse){
      combination.result <<- c(combination.result, plotdata$seq[i])
    }
  }

  #err <- normalizeData(as.vector(plotdata$err),type="0_1")
  data.to.plot <- data.frame(err.rmse=plotdata$err, dist=avg.dist)
  return(data.to.plot)
}

plot.for.algo <- function(){
  par(mfrow=c(5,2))
  par(mar=c(0.5, 4.5, 0.5, 0.5))
  plot.window(xlim=c(0,300),ylim=c(0,1),asp=1)

  #plot1-10
  setwd(file.path.all)
  for(i in seq(1,sites.count)){
    load.err.data(i)
    data.to.plot <- get.err.dist.data()
    #testing <<- data.to.plot##remove this line later

    y <- data.to.plot$err.rmse
    x <- data.to.plot$dist
    d <- data.frame(y,x)
    f <- formula('y~x')

    plot(f,d,
         main=cor(y,x),#method="spearman"),
         xlab=paste("distance combination ",i),
         ylab=paste("rmse combination ",i))

    dev.copy2pdf(file =plot.file)
    #dev.off()
    combination.result <<- as.character(combination.result)
    combination.result <<- gsub("0","10",combination.result)
    df <- c(paste("RMSE threshold", as.character(threshold.rmse)))
    df <- c(df,paste("Similarity Threshold",as.character(threshold.dist)))
    df <- c(df, combination.result)
    write.table(df, result.file)
  }
}

plot.all <- function(){
  typ <- 'statistical'
  load.data()
  for(folder in folder.ip.type.vec){
    for(err in err.type.vec){
      for(mes in simi.meas.vec){
        gen.dist.mat(mes)
        for(alg in algo.vec){
          setvals(alg,mes,typ,folder,err)
          plot.for.algo()
          #break
        }
        #break
      }
    }
  }
}

plot.all()
