require(gbm)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(combinat)
require(RSNNS)
require(ftsa)

#aggr.cluster.size <- 3
sites.count <- 10
#indxcombicnt <- 10 #no. of combination
#aggr.mat.size <- indxcombicnt#floor(sites.count/aggr.cluster.size)
data.len <- 52560
data.len.day <<- 144
mat.size <<- 365
window.size <- 10
train.data.percent <- 0.7
slide.count <- mat.size-window.size+1
#indxseq <- c(seq(1,sites.count))
filepath <<- '/home/freak/Programming/Thesis/results/results/gbm_shortterm_aggr/all/'

powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
#aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
#train.data <<- c()
#test.data <<- c()
#output <<- c()
#indxcombimat <<- ff(NA, dim=c(aggr.mat.size,indxcombicnt),vmode="integer")

drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")
tablelist_statement = paste("SELECT TABLE_NAME FROM information_schema.TABLES ",
                            "WHERE TABLE_SCHEMA = 'eastwind' AND",
                            "TABLE_NAME LIKE 'onshore_SITE_%' "," LIMIT ",sites.count, ";")
tables <- dbGetQuery(con, statement=tablelist_statement)
tables <- data.frame(tables)

loaddata <- function(){
  colindx <- 1
  for(indx in seq(1,sites.count)){
    tab <- tables[indx,]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
  }
}

generate.seq.matrix <- function(){
  mat <- as.matrix(combn(sites.count, aggr.cluster.size)) #generating different combination
  len <- length(mat[1,])
  indxseq <- sample(1:len, indxcombicnt)
  indxcombimat <<- mat[,indxseq]
}

gen.aggrdata <- function(){
  col.names <- c()
  if(aggr.mat.size != 0){
    for(i in seq(1 ,aggr.mat.size)){
      indx.seq <- indxcombimat[,i]
      mat <- as.matrix(powdata[,indx.seq])
      col.names <- c(col.names, paste('S',paste(indx.seq, collapse=""), sep=""))
      aggrdata[,i] <<- rowSums(mat)
    }
    colnames(aggrdata) <<- col.names
  }
  else{
    indx.seq <- seq(1,sites.count)
    print(indx.seq)
    aggrdata10 <<- rowSums(as.matrix(powdata[,indx.seq]))
  }

}


predict <- function(aggrno, indx) {
  if(indx < 1 || indx >= data.len){
    print("Enter indx Greater than 0 and less than the data size")
    return
  }

  data.normalized <- c()
  if(aggr.mat.size!=0){
    data.normalized <- normalizeData(as.vector(aggrdata[,aggrno]),type="0_1")
  }
  else{
    data.normalized <- normalizeData(as.vector(aggrdata10),type="0_1")
  }
  data.set <- as.vector(data.normalized[,1])
  #data <- c()
  #if(aggr.mat.size!=0){
  #  data <- as.vector(aggrdata[,aggrno])
  #}
  #else{
  # data <- as.vector(aggrdata10)
  #}
  #data.set <- data#as.vector(data[,1])

  indx.start <<- indx
  indx.end <<- indx.start + (window.size * data.len.day) - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
    print(paste("Cluster size: ",aggr.cluster.size," AggrNo.: ", aggrno, " Slide No.: ", count+1))
    data.mat <- matrix(data.set[indx.start:indx.end], nrow=data.len.day, ncol=window.size, byrow=FALSE)
    colnames(data.mat) <- paste("d",c(1:window.size), sep="")

    train.dataset.indx <- floor(data.len.day *  train.data.percent)
    test.dataset.indx <- train.dataset.indx + 1
    window.slide <- data.len.day - train.dataset.indx
    data.mat.train <- data.frame(data.mat[1:train.dataset.indx,])
    data.mat.test <- data.frame(data.mat[test.dataset.indx:data.len.day,])

    formula.set <- colnames(data.mat)
    y = formula.set[window.size]
    x = formula.set[1:window.size-1]
    f = as.formula(paste(y, " ~ ", paste(x, collapse="+")))

    out <<- gbm(f,
                data=data.mat.train,
                distribution ="gaussian",
                n.trees=10000,
                interaction.depth = 10,
                n.minobsinnode = 5,
                shrinkage =  0.008,
                n.cores=3)

    data.train <<- c(data.train, data.mat.train[,window.size])
    data.test <<- c(data.test, data.mat.test[,window.size])

    pred <- predict.gbm(out, data.mat.test[,1:window.size-1], out$n.trees)
    data.out <<- c(data.out, pred)

    indx.start <<- indx.start + window.slide
    indx.end <<- indx.start + (window.size * data.len.day)
    count <- count + 1
    if(count == 10){
      break
    }
  }

  train.data <<- cbind(train.data, data.train) #data.mat.train[,window.size]
  test.data <<-  cbind(test.data, data.test) #data.mat.test[,window.size]
  output <<- cbind(output, data.out)
}


predict.for.combination <- function(){
  slide.indx <- 45361
  #loaddata()
  #generate.seq.matrix()
  gen.aggrdata()
  if(aggr.mat.size != 0){
    for(aggr.indx in seq(1,aggr.mat.size)){
      predict(aggr.indx,slide.indx)
      break
    }
  }
  else{
    predict(0,slide.indx)
  }

}

prediction.error <- function(){
  parm.count <- 5
  file.name <- paste('gbm_shortterm_aggr_combi',aggr.cluster.size,'.csv',sep="")
  setwd(filepath)
  #file <- paste(filepath,file.name)
  #print(file)
  if(aggr.mat.size!=0){
    err.data <<- matrix(,nrow=indxcombicnt, ncol=parm.count, byrow=TRUE)
    colnames(err.data) <<- c("AggrNo.Seq","rmse", "mape", "sse", "mse")
    col.names <- colnames(aggrdata)
    for(site in seq(1:(aggr.mat.size))){
      site.name <- col.names[site]
      test <- test.data[,site]
      pred <- output[,site]
      err.rmse <- error(forecast=pred, true=test,method="rmse")
      err.mape <- error(forecast=pred, true=test,method="mape")
      err.sse <- error(forecast=pred, true=test,method="sse")
      err.mse <- error(forecast=pred, true=test,method="mse")
      err.data[site,] <<- c(site.name, err.rmse, err.mape, err.sse, err.mse)
      break
    }
    write.csv(err.data, file=file.name)
  }
  else{
    err.data <<- matrix(,nrow=1,ncol=parm.count, byrow=TRUE)
    colnames(err.data) <<- c("AggrNo.Seq","rmse", "mape", "sse", "mse")
    col.name <- paste('S',paste(seq(1,aggr.cluster.size), collapse=""), sep="")
    site.name <- col.name
    test <- test.data[,1]
    pred <- output[,1]
    err.rmse <- error(forecast=pred, true=test,method="rmse")
    err.mape <- error(forecast=pred, true=test,method="mape")
    err.sse <- error(forecast=pred, true=test,method="sse")
    err.mse <- error(forecast=pred, true=test,method="mse")
    err.data[1,] <<- c(site.name, err.rmse, err.mape, err.sse, err.mse)

    write.csv(err.data, file=file.name)
  }

}

predict.all.combination <- function(){
  loaddata()
  for(combi in seq(2,10)){#sites.count
    aggr.cluster.size <<- combi

    if(combi != sites.count){
      indxcombicnt <<-length(combn(sites.count,combi)[1,])
      aggr.mat.size <<- indxcombicnt
      aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
      indxcombimat <<- ff(NA, dim=c(combi,indxcombicnt),vmode="integer")
      generate.seq.matrix()
    }
    else{
      indxcombicnt <<- 0
      aggr.mat.size <<- 0
      aggrdata10 <<- ff(NA, dim=data.len, vmode="double")
    }

    train.data <<- c()
    test.data <<- c()
    output <<- c()

    predict.for.combination()
    prediction.error()
  }

}

predict.all.combination()

#plotting
length(test.data[,1])
x1 = train.data[,1]
x2 = test.data[,1]
y = output[,1]
length(y)
plot(x2, type="l")

dataToPlot = data.frame(seq(1,440),x2, y)
Line <- gvisLineChart(dataToPlot)
plot(Line)


#error measure
err <- error(forecast=y, true=x2,method="mape")
err
