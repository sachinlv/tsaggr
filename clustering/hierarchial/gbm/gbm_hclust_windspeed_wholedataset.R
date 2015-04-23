require(RMySQL)
require(RSNNS)
require(TSdist)
require(ff)
require(gbm)
require(ftsa)
require(matrixStats)

history.length <- 10
sites.count <- 100
data.len.day <<- 144
data.len <- history.length * data.len.day
window.size <- 144
start.date <- '20061222'
end.date <- '20070101'
train.data.percent <- 0.96
hidden.nodes <<- 10

plot.file.generic <- 'gbm_shortterm_windspeed_hclust_'
plot.file.path <- '/home/freak/Programming/Thesis/results/plots/hclust/wholedataset/'
file.name.generic <<-'gbm_shortterm_hclust'
forecast.aggr.err.file <<- 'gbm_shortterm_hclust_aggr.csv'
filepath.generic <<- '/home/freak/Programming/Thesis/results/results/gbm_wholedataset/'
sim.meas <<- "euclidean"
powdata <<- matrix(0, nrow=data.len, ncol=sites.count)
winddata <<- matrix(0, nrow=data.len, ncol=sites.count)
dist.mat <<- matrix(0, nrow=sites.count, ncol=sites.count)

drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")

loaddata <- function(){
  tbllist_qry = paste("SELECT TABLE_NAME FROM information_schema.TABLES ",
                      "WHERE TABLE_SCHEMA = 'eastwind' AND",
                      "TABLE_NAME LIKE 'onshore_SITE_%' LIMIT 100;")
  tabls <- data.frame(dbGetQuery(con,statement=tbllist_qry), check.names=FALSE)
  tabls <- tabls[,1]

  colindx <- 1
  for(indx in seq(1,sites.count)){
    print(paste("Loading from table :: ", tabls[indx]))
    query <- paste(" select pow,spd from ", tabls[indx], " WHERE (mesdt >= ",start.date," && mesdt < ",end.date,") LIMIT ", data.len, ";")
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
      dist.mat[i,j] <<- euclideanDistance(data.mat$c1,data.mat$c2)
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
    col.names <- c(col.names, paste('Cut',cut.no, sep=""))
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

    out <<- gbm(f,
                data=trn.data,
                distribution ="gaussian",
                n.trees=1000,
                interaction.depth = 10,
                n.minobsinnode = 5,
                shrinkage =  0.008,
                n.cores=3)

    data.train <<- c(data.train, trn.data$y)
    data.test <<- c(data.test, tst.y$y)

    pred <- predict.gbm(out, tst.x, out$n.trees)
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

measure.error <- function(pred,test){
  err.rmse <- error(forecast=pred, true=test,method="rmse")
  err.mape <- error(forecast=pred, true=test,method="mape")
  err.mae <- error(forecast=pred, true=test,method="mae")
  err.mse <- error(forecast=pred, true=test,method="mse")
  err.sd <- sd(pred-test)# need to normalize it to installed power
  err.cor <- cor(pred,test)

  return(c(err.rmse, err.mape, err.mae, err.mse, err.sd, err.cor))
}


prediction.error <- function(time.taken){
  parm.count <- 7
  file.name <- paste(file.name.generic,cut.size,'.csv',sep="")
  setwd(filepath.generic)

  if(cut.size!=0){
    err.data <<- matrix(NA,nrow=cut.size, ncol=parm.count, byrow=TRUE)
    colnames(err.data) <<- c("Clust.Id","rmse", "mape", "mae", "mse", "sd", "cor")
    col.names <- colnames(aggr.mat)
    for(clust in seq(1:cut.size)){
      site.name <- col.names[clust]
      err.data[clust,] <<- c(site.name, measure.error(output[,clust], test.data[,clust]))
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
    err.row <- c(cut.size,max.clust.size, measure.error(output.sum,test.data.sum),time.taken)
    aggr.err[cut.size,] <<- err.row
  }
}

predict.for.aggrdata <- function(){
  train.data <<- c()
  test.data <<- c()
  test.data.denorm <<- c()
  output <<- c()
  output.denorm <<- c()

  time.start <- Sys.time()
  for(aggr.nr in seq(1,cut.size)){
    predict.pow(aggr.nr)
  }
  time.end <- Sys.time()
  time.taken <- as.numeric(time.end - time.start, units="secs")

  prediction.error(time.taken)
}

plot.err.aggr <- function(){
  setwd(plot.file.path)
  err.tbl <- data.frame(aggr.err)
  x <- err.tbl$maxclustsize
  y <- err.tbl$rmse
  plot(y~x)
  dev.copy2pdf(file =paste(plot.file.generic,sim.meas,'.pdf',sep=""))
}

predict.hclust.aggregates <- function(){
  #for(sim in simi.meas.vec){
    #sim.meas <<- sim
    generate.hierarchial.cluster()
    cut.tree.mat <- cutree(hc,k=1:sites.count)
    #no.of.cuts <- length(cut.tree.mat[1,])
    no.of.cuts <- c(1,10,20,30,40,50,60,70,80,90,100)
    aggr.err <<- matrix(0,ncol=9,nrow=len(no.of.cuts),byrow=TRUE, dimnames=NULL)
    colnames(aggr.err) <<- c("cutsize", "maxclustsize","rmse", "mape", "mae", "mse", "sd", "cor", "exectime")

    for(cut in no.of.cuts){
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
  #}
}

predict.hclust.aggregates()
