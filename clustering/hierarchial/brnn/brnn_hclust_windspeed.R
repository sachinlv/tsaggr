require(RMySQL)
require(RSNNS)
require(TSdist)
require(ff)
require(brnn)
require(ftsa)
require(matrixStats)

history.length <- 10
sites.count <- 10
data.len <- 1440#52560
data.len.day <<- 144
window.size <- 144
start.date <- '20061222'
end.date <- '20070101'
train.data.percent <- 0.96
hidden.nodes <<- 10

simi.meas.vec <<- c(
  'euclidean',
  'minkowski',
  'manhattan',
  'fourier',
  'correlation',
  'pca',
  'coord',
  'lm',
  'weuclid',
  'maj',
  'shrinkage',
  'pdist')

plot.file.generic <- 'brnn_shortterm_windspeed_hclust_'
plot.file.path <- '/home/freak/Programming/Thesis/results/plots/history10/specific_sites/clustering/'
file.name.generic <<-'brnn_shortterm_hclust'
forecast.aggr.err.file <<- 'brnn_shortterm_hclust_aggr.csv'
filepath.generic <<- '/home/freak/Programming/Thesis/results/results/brnn_shortterm_hclust/'

table.ip.type <- "random"#c("random","specific")
powdata <<- matrix(0, nrow=data.len, ncol=sites.count)
winddata <<- matrix(0, nrow=data.len, ncol=sites.count)
dist.mat <<- matrix(0, nrow=sites.count, ncol=sites.count)

drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")

if(table.ip.type == "random"){
  t <- c('onshore_SITE_00002',
         'onshore_SITE_00003',
         'onshore_SITE_00004',
         'onshore_SITE_00005',
         'onshore_SITE_00006',
         'onshore_SITE_00007',
         'onshore_SITE_00008',
         'onshore_SITE_00012',
         'onshore_SITE_00013',
         'onshore_SITE_00014')
  tables <<- data.frame(cbind(numeric(0),t))
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
  tables <<- data.frame(cbind(numeric(0),t))
}

loaddata <- function(){
  colindx <- 1
  for(indx in seq(1,sites.count)){
    tab <- tables[indx,1]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow,spd from ", tab, " WHERE (mesdt >= ",start.date," && mesdt < ",end.date,") LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    winddata[,indx] <<- as.double(data06[,2])
  }
}

generate.distance.matrix <- function(){
  dist.mat <<- matrix(0, nrow=sites.count, ncol=sites.count)
  for(i in seq(1,sites.count)){
    for(j in seq(1, sites.count)){
      data.mat <- data.frame(c1=powdata[,i], c2=powdata[,j])

      dist.mat[i,j] <<- switch(sim.meas,
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

generate.hierarchial.cluster <- function(){
  loaddata()
  generate.distance.matrix()
  hc <<- hclust(as.dist(dist.mat))
  plot(hc)
}


generate.aggr.matrix <- function(size,vect){
  aggr.mat <<- matrix(0,nrow=data.len,ncol=size)
  col.names <- c()
  clust.sizes <- c()
  for(cut.no in seq(1,size)){
    indx.seq <- which(vect==cut.no)

    pmat <- as.matrix(powdata[,indx.seq])
    wmat <- as.matrix(winddata[,indx.seq])
    clust.sizes <- c(clust.sizes,length(indx.seq))
    col.names <- c(col.names, paste('S',paste(indx.seq, collapse=""), sep=""))
    if(length(pmat[1,])>1){
      pow.aggrdata[,cut.no] <<- rowSums(pmat)
      wind.aggrdata[,cut.no] <<- rowWeightedMeans(wmat)
    }else{
      pow.aggrdata[,cut.no] <<- pmat[,1]
      wind.aggrdata[,cut.no] <<- wmat[,1]
    }
  }
  colnames(aggr.mat) <<- col.names
  return(max(clust.sizes))
}

predict.pow <- function(aggrno) {
  pow.data.normalized <- normalizeData(as.vector(pow.aggrdata[,aggrno]),type="0_1")
  wind.data.normalized <- normalizeData(as.vector(wind.aggrdata[,aggrno]),type="0_1")
  pow.normParms <<- getNormParameters(pow.data.normalized)
  wind.normParms <<- getNormParameters(wind.data.normalized)

  pow.data.set <- as.vector(pow.data.normalized[,1])
  wind.data.set <- as.vector(wind.data.normalized[,1])

  indx.start <<- 1
  indx.end <<- indx.start + window.size - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
    print(paste("Cut size: ",cut.size," AggrNo.: ", aggrno, " Slide No.: ", count+1))
    y <- as.vector(pow.data.set[indx.start:indx.end])
    x <- as.vector(wind.data.set[indx.start:indx.end])
    dat <- data.frame(cbind(y,x))

    train.indx <- floor(window.size *  train.data.percent)
    test.indx <- train.indx + 1
    window.slide <- window.size - train.indx

    trn.data <- data.frame(dat[1:train.indx,])
    tst.x <- data.frame(x=dat$x[test.indx:window.size])
    tst.y <- data.frame(y=dat$y[test.indx:window.size])
    f = as.formula("y ~ x")

    out <<- brnn(f,
                 trn.data,
                 epochs=5,
                 cores=2,
                 mu=0.1,
                 mu_dec=0.1,
                 mu_max=1e10,
                 change = 0.001,
                 neurons=hidden.nodes,
                 normalize=FALSE,
                 verbose=FALSE,
                 Monte_Carlo = FALSE)

    data.train <<- c(data.train, trn.data$y)
    data.test <<- c(data.test, tst.y$y)

    pred <- predict.brnn(out , tst.x)
    data.out <<- c(data.out, pred)

    indx.start <<- indx.start + window.slide
    indx.end <<- indx.start + window.size
    count <- count + 1
  }

  train.data <<- cbind(train.data, data.train)
  test.data <<-  cbind(test.data, data.test)
  output <<- cbind(output, data.out)
  test.data.denorm <<- cbind(test.data.denorm, as.vector(denormalizeData(data.test,pow.normParms)))
  output.denorm <<- cbind(output.denorm, as.vector(denormalizeData(data.out, pow.normParms)))
}



prediction.error <- function(){
  parm.count <- 7
  filepath <- paste(filepath.generic,sim.meas,'/',sep="")
  file.name <- paste(file.name.generic,cut.size,'.csv',sep="")
  setwd(filepath)

  if(cut.size!=0){
    err.data <<- matrix(NA,nrow=cut.size, ncol=parm.count, byrow=TRUE)
    colnames(err.data) <<- c("Clust.Id","rmse", "mape", "mae", "mse", "sd", "cor")
    col.names <- colnames(aggr.mat)
    for(clust in seq(1:cut.size)){
      site.name <- col.names[clust]
      test <- test.data[,clust]
      pred <- output[,clust]
      err.rmse <- error(forecast=pred, true=test,method="rmse")
      err.mape <- error(forecast=pred, true=test,method="mape")
      err.mae <- error(forecast=pred, true=test,method="mae")
      err.mse <- error(forecast=pred, true=test,method="mse")
      err.sd <- sd(pred-test)# need to normalize it to installed power
      err.cor <- cor(pred,test)
      err.data[clust,] <<- c(site.name,err.rmse, err.mape, err.mae, err.mse, err.sd, err.cor)
    }
    write.csv(err.data, file=file.name)

    if(cut.size > 1){
      test.data.sum <<- rowSums(test.data.denorm[,1:cut.size])
      test.data.sum <<- normalizeData(test.data.sum, type="0_1")

      output.sum <<- rowSums(output.denorm[,1:cut.size])
      output.sum <<- normalizeData(output.sum, type="0_1")
    }else{
      test.data.sum <<- test.data.denorm[,1]
      test.data.sum <<- normalizeData(test.data.sum, type="0_1")

      output.sum <<- output.denorm[,1]
      output.sum <<- normalizeData(output.sum, type="0_1")
    }


    err.rmse <- error(forecast=pred, true=test,method="rmse")
    err.mape <- error(forecast=pred, true=test,method="mape")
    err.mae <- error(forecast=pred, true=test,method="mae")
    err.mse <- error(forecast=pred, true=test,method="mse")
    err.sd <- sd(pred-test)# need to normalize it to installed power
    err.cor <- cor(pred,test)
    err.row <- c(cut.size,max.clust.size,
                 err.rmse,err.mape,err.mae,err.mse,err.sd,err.cor)
    aggr.err[cut.size,] <<- err.row
  }
}

predict.for.aggrdata <- function(){
  train.data <<- c()
  test.data <<- c()
  test.data.denorm <<- c()
  output <<- c()
  output.denorm <<- c()

  for(aggr.nr in seq(1,cut.size)){
    predict.pow(aggr.nr)
  }
  prediction.error()
}

plot.err.aggr <- function(){
  setwd(plot.file.path)
  err.tbl <- data.frame(aggr.err)
  x <- err.tbl$max.clust.size
  y <- err.tbl$rmse
  plot(y~x)
  dev.copy2pdf(file =paste(plot.file.generic,sim.meas,'.pdf',sep=""))
}

predict.hclust.aggregates <- function(){
  for(sim in simi.meas.vec){
    sim.meas <<- sim
    generate.hierarchial.cluster()
    cut.tree.mat <- cutree(hc,k=1:sites.count)
    no.of.cuts <- length(cut.tree.mat[1,])
    aggr.err <<- matrix(0,ncol=8,nrow=sites.count,byrow=TRUE, dimnames=NULL)
    colnames(aggr.err) <- c("cutsize", "maxclustsize","rmse", "mape", "mae", "mse", "sd", "cor")

    for(cut in seq(1,no.of.cuts)){
      cut.size <<- cut
      cut.group.vec <- cut.tree.mat[,cut]
      pow.aggrdata <<- ff(NA, dim=c(data.len, cut), vmode="double")
      wind.aggrdata <<- ff(NA, dim=c(data.len, cut), vmode="double")
      max.clust.size <<- generate.aggr.matrix(cut, cut.group.vec)
      predict.for.aggrdata()
      #break
    }
    filepath <- paste(filepath.generic,sim.meas,'/',sep="")
    setwd(filepath)
    write.csv(aggr.err,forecast.aggr.err.file)
    plot.err.aggr()
  }
}

predict.hclust.aggregates()
