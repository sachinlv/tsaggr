require(RMySQL)
require(ff)
require(neuralnet)
require(brnn)
require(earth)
require(gbm)
require(TSdist)
require(forecast)
require(timeSeries)
require(TSclust)
require(corpcor)
require(MADAM)
require(Metrics)
require(ppls)
require(combinat)
require(RSNNS)
require(ftsa)
require(zoo)
require(matrixStats)
require(forecast)
require(timeSeries)

sites.count <- 10
history.length <- 10
data.len.day <<- 144
data.len <- history.length * data.len.day
window.size <- 144
train.data.percent <- 0.96
start.date <- '20061112'
end.date <- '20070101'
hidden.nodes <<- 10
powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
winddata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")

table.ip.type <- 'random'
dist.mat.file.path <-'/home/freak/Programming/Thesis/results/resultsreview/random_sites/distance/'
error.path <- '/home/freak/Programming/Thesis/results/resultsreview/random_sites/error/'

drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")

if(table.ip.type == "random"){
  t <-   c("onshore_SITE_00002",
           "onshore_SITE_00003",
           "onshore_SITE_00004",
           "onshore_SITE_00005",
           "onshore_SITE_00006",
           "onshore_SITE_00007",
           "onshore_SITE_00008",
           "onshore_SITE_00012",
           "onshore_SITE_00013",
           "onshore_SITE_00014")
  tables <<- data.frame(cbind(numeric(0),t))
}else{
  t <- c("onshore_SITE_04468",
         "onshore_SITE_04476",
         "onshore_SITE_04665",
         "onshore_SITE_04640",
         "onshore_SITE_04290",
         "onshore_SITE_04181",
         "onshore_SITE_04607",
         "onshore_SITE_04247",
         "onshore_SITE_04605",
         "onshore_SITE_04094")
  tables <<- data.frame(cbind(numeric(0),t))
}

algo.vec <<- c('neuralnet','brnn','gbm','mars')
simi.meas.vec <<- c(
  'euclidean',
  'fourier',
  'pca',
  'lm',
  'weuclid')


loaddata <- function(){
  colindx <- 1
  for(indx in seq(1,sites.count)){
    tab <- tables[indx,]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow,spd from ", tab, " WHERE (mesdt >= ",start.date," && mesdt < ",end.date,") LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    winddata[,indx] <<- as.double(data06[,2])
  }
}


gen.dist.mat <- function(measure){
  print(paste("Distance measure: ", measure))
  
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
                               "fourier"={
                                 fourierDistance(data.mat$c1,data.mat$c2,n=(floor(length(data.mat$c1)/2)+1))
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
                               }
      )
    }
  }
}


gen.all.distmat <- function(){
  for(mes in simi.meas.vec){
    dist.mat <<- matrix(0, nrow=sites.count, ncol=sites.count)
    gen.dist.mat(mes)
    result.file <- paste(dist.mat.file.path,
                         'distance_matrix_',
                         mes,'.csv',sep="")
    write.table(dist.mat, result.file)
  }
}


predict.pow <- function(algo.type, site) {
  pow.data.normalized <-normalizeData(as.vector(powdata[,site]),type="0_1")
  wind.data.normalized <- normalizeData(as.vector(winddata[,site]),type="0_1")
  
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
    print(paste(" SiteNo.: ", site, " Slide No.: ", count+1))
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
    
    pred  <-  switch(algo.type,
                     "neuralnet"={
                       out <- neuralnet(f,
                                        trn.data,
                                        hidden=hidden.nodes,
                                        rep=2,
                                        stepmax = 2000,
                                        threshold=0.2,
                                        learningrate=1,
                                        algorithm="rprop+", #'rprop-', 'sag', 'slr'
                                        startweights=NULL,
                                        lifesign="none",
                                        err.fct="sse",
                                        act.fct="logistic",
                                        exclude = NULL,
                                        constant.weights = NULL,
                                        linear.output=TRUE)
                       compute(out, tst.x)$net.result
                     },
                     "gbm"={
                       out <- gbm(f,
                                  data=trn.data,
                                  distribution ="gaussian",
                                  n.trees=1000,
                                  interaction.depth = 10,
                                  n.minobsinnode = 5,
                                  shrinkage =  0.008,
                                  n.cores=3)
                       predict.gbm(out, tst.x, out$n.trees)
                     },
                     "brnn"={
                       out <- brnn(f,
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
                       predict.brnn(out , tst.x)
                     },
                     "mars"={
                       out <<- earth(f,
                                     data=trn.data)
                       predict(out , tst.x)
                     }
    )
    
    data.test <<- c(data.test, tst.y$y)
    data.out <<- c(data.out, pred)
    
    indx.start <<- indx.start + window.slide
    indx.end <<- indx.start + window.size
    count <- count + 1
  }
  
  test.data <<-  cbind(test.data, data.test)
  output <<- cbind(output, data.out)
}


measure.error <- function(pred,test){
  err.rmse <- error(forecast=pred, true=test,method="rmse")
  err.mse <- error(forecast=pred, true=test,method="mse")
  err.sd <- sd(pred-test)
  bias.sqr <- (mean(pred) - mean(test))^2
  pred.var <- var(pred)
  
  return(c(err.rmse, err.mse, err.sd, bias.sqr, pred.var))
}


predict.all <- function(){
  for(algo in algo.vec){
    err.data <- matrix(0,ncol=5,nrow=sites.count,byrow=TRUE, dimnames=NULL)
    colnames(err.data) <- c("rmse", "mse", "sd", "bias.sqr", "pred.var")
    rownames(err.data) <- paste("S", seq(1,sites.count),sep="")
    
    for(site in seq(1,sites.count)){
      test.data <<- c()
      output <<- c()
      predict.pow(algo,site)
      err.data[site,] <- measure.error(output,test.data)
    }
    error.file <- paste(error.path, algo,".csv", sep="")
    write.table(err.data, error.file)
  }
}


cal.dist.and.predict <- function(){
  loaddata()
  gen.all.distmat()
  predict.all()
}

cal.dist.and.predict()

