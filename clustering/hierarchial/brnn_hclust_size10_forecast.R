require(RMySQL)
require(RSNNS)
require(TSdist)
require(ff)
require(brnn)
require(ftsa)

history.length <- 50
sites.count <- 10
data.len <- 52560
data.len.day <<- 144
mat.size <<- 365
window.size <- 10
sim.meas <- 'pca'
plot.file <- 'brnn_shortterm_hclust_aggr_plot.pdf'
train.data.percent <- 0.7
file.name.generic <<-'brnn_shortterm_hclust'
forecast.aggr.err.file <<- 'brnn_shortterm_hclust_aggr.csv'
filepath <<- '/home/freak/Programming/Thesis/results/results/brnn_shortterm_hclust/'
aggr.err <<- matrix(0,ncol=6,nrow=sites.count,byrow=TRUE, dimnames=NULL)
colnames(aggr.err) <- c('cut size', 'max clust size','rmse','mape','sse','mse')

table.ip.type <- "random"#"specific"
powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
dist.mat <<- matrix(0, nrow=sites.count, ncol=sites.count)
dist.mat.norm <<- matrix(0, nrow=sites.count, ncol=sites.count)


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
    tab <- as.character(tables[indx,1])
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]), type="0_1")
  }
}


generate.distance.matrix <- function(){
  for(i in seq(1,sites.count)){
    for(j in seq(1, sites.count)){
      data.mat <- data.frame(c1=powdata[,i], c2=powdata[,j])

      dist.mat[i,j] <- switch(sim.meas,
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
                              "mahalanobis"={
                                d <- matrix(NA,nrow=52560,ncol=2)
                                #d[,1] <- powdata[,i1]
                                #d[,2] <- powdata[,j1]
                                #c <- cov(d[,1:2])#,method=c("pearson"))
                                #m <- c(mean(d[,1]), mean(d[,2]))
                                #mat[i,j] <- mahalanobis(d,m,c)
                              }
      )
    }
  }
}


generate.hierarchial.cluster <- function(){
  loaddata()
  generate.distance.matrix()
  hc <<- hclust(as.dist(dist.mat))
  hc.norm <<- hclust(as.dist(dist.mat.norm))
}


generate.aggr.matrix <- function(size,vect){
  aggr.mat <<- matrix(0,nrow=52560,ncol=size)
  col.names <- c()
  clust.sizes <- c()
  for(cut.no in seq(1,size)){
    indx.seq <- which(vect==cut.no)
    mat <- as.matrix(powdata[,c(indx.seq)])
    clust.sizes <- c(clust.sizes,length(indx.seq))
    col.names <- c(col.names, paste('S',paste(indx.seq, collapse=""), sep=""))
    if(length(mat[1,])>1){
      aggr.mat[,cut.no] <<- rowSums(mat)
    }else{
      aggr.mat[,cut.no] <<- mat[,1]
    }
  }
  colnames(aggr.mat) <<- col.names
  return(max(clust.sizes))
}

predict.pow <- function(cut.size, aggrno, indx) {
  if(indx < 1 || indx >= data.len){
    print("Enter indx Greater than 0 and less than the data size")
    return
  }
  data.normalized <- normalizeData(as.vector(aggr.mat[,aggrno]),type="0_1")
  normParms <<- getNormParameters(data.normalized)
  data.set <- as.vector(data.normalized[,1])
  indx.start <<- indx
  indx.end <<- indx.start + (window.size * data.len.day) - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
    print(paste("Cut size: ",cut.size," AggrNo.: ", aggrno, " Slide No.: ", count+1))
    data.mat <- matrix(data.set[indx.start:indx.end], nrow=data.len.day, ncol=window.size, byrow=FALSE)
    colnames(data.mat) <- paste("d",c(1:window.size), sep="")

    train.dataset.indx <- floor(data.len.day *  train.data.percent)
    test.dataset.indx <- train.dataset.indx + 1
    window.slide <- data.len.day - train.dataset.indx
    data.mat.train <- data.mat[1:train.dataset.indx,]
    data.mat.test <- data.mat[test.dataset.indx:data.len.day,]

    formula.set <- colnames(data.mat)
    y = formula.set[window.size]
    x = formula.set[1:window.size-1]
    f = as.formula(paste(y, " ~ ", paste(x, collapse="+")))

    out <<- brnn(f,
                 data.mat.train,
                 epochs=1000,
                 cores=2,
                 mu=0.1,
                 mu_dec=0.1,
                 mu_max=1e10,
                 change = 0.01,
                 neurons=window.size,
                 normalize=TRUE,
                 verbose=FALSE,
                 Monte_Carlo = FALSE)


    data.train <<- c(data.train, data.mat.train[,window.size])
    data.test <<- c(data.test, data.mat.test[,window.size])

    pred <- predict.brnn(out ,data.mat.test[, 1:window.size-1])
    data.out <<- c(data.out, pred)

    indx.start <<- indx.start + window.slide
    indx.end <<- indx.start + (window.size * data.len.day)
    count <- count + 1
    #if(count == 10){
    #  break
    #}
  }

  train.data <<- cbind(train.data, data.train)
  test.data <<-  cbind(test.data, data.test)
  test.data.denorm <<- cbind(test.data.denorm, as.vector(denormalizeData(data.test, normParms)))
  output <<- cbind(output, data.out)
  output.denorm <<- cbind(output.denorm, as.vector(denormalizeData(data.out, normParms)))
}


prediction.error <- function(cut.size,max.clust.size){
  parm.count <- 5

  file.name <- paste(file.name.generic,cut.size,'.csv',sep="")

  setwd(filepath)

  if(cut.size!=0){
    err.data <<- matrix(NA,nrow=cut.size, ncol=parm.count, byrow=TRUE)
    colnames(err.data) <<- c("Clust.Id","rmse", "mape", "sse", "mse")
    col.names <- colnames(aggr.mat)
    for(clust in seq(1:cut.size)){
      site.name <- col.names[clust]
      test <- test.data[,clust]
      pred <- output[,clust]
      err.rmse <- error(forecast=pred, true=test,method="rmse")
      err.mape <- error(forecast=pred, true=test,method="mape")
      err.sse <- error(forecast=pred, true=test,method="sse")
      err.mse <- error(forecast=pred, true=test,method="mse")
      err.data[clust,] <<- c(site.name, err.rmse, err.mape, err.sse, err.mse)
      #break
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


    err.rmse <- error(forecast=output.sum, true=test.data.sum,method="rmse")
    err.mape <- error(forecast=output.sum, true=test.data.sum,method="mape")
    err.sse <- error(forecast=output.sum, true=test.data.sum,method="sse")
    err.mse <- error(forecast=output.sum, true=test.data.sum,method="mse")

    err.row <- c(cut.size,max.clust.size, err.rmse, err.mape, err.sse, err.mse)
    aggr.err[cut.size,] <<- err.row
  }
}

predict.for.aggrdata <- function(cut.size,max.clust.size){
  slide.indx <- data.len - (history.length * data.len.day) + 1
  train.data <<- c()
  test.data <<- c()
  test.data.denorm <<- c()
  output <<- c()
  output.denorm <<- c()

  for(aggr.nr in seq(1,cut.size)){
    predict.pow(cut.size,aggr.nr,slide.indx)
  }
  prediction.error(cut.size,max.clust.size)
}

plot.err.aggr <- function(){
  setwd(filepath)
  err.tbl <- data.frame(read.csv(forecast.aggr.err.file))
  x <- err.tbl$max.clust.size
  y <- err.tbl$rmse
  plot(y~x)
  dev.copy2pdf(file =plot.file)
}

predict.hclust.aggregates <- function(){
  generate.hierarchial.cluster()
  plot(hc)
  cut.tree.mat <- cutree(hc,k=1:sites.count)
  no.of.cuts <- length(cut.tree.mat[1,])

  for(cut.size in seq(1,no.of.cuts)){
    cut.group.vec <- cut.tree.mat[,cut.size]
    max.clust.size <- generate.aggr.matrix(cut.size,cut.group.vec)
    predict.for.aggrdata(cut.size,max.clust.size)
    #break
  }
  setwd(filepath)
  write.csv(aggr.err,forecast.aggr.err.file)
  plot.err.aggr()
}

predict.hclust.aggregates()
