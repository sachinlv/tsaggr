require(neuralnet)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(combinat)
require(RSNNS)
require(ftsa)


history.length <- 50
sites.count <- 10
hidden.nodes <<- 10#c(round(window.size/2), window.size,1)
data.len <- 52560
data.len.day <<- 144
mat.size <<- 365
window.size <- 1440
train.data.percent <- 0.7
#slide.count <- mat.size-window.size+1
filepath <<- '/home/freak/Programming/Thesis/results/results/neuralnet_shortterm_windspeed_aggr/all/'

powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
winddata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
table.ip.type <- "specific"#"random"

drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")

if(table.ip.type == "random"){
  tablelist_statement = paste("SELECT TABLE_NAME FROM information_schema.TABLES ",
                              "WHERE TABLE_SCHEMA = 'eastwind' AND",
                              "TABLE_NAME LIKE 'onshore_SITE_%' "," LIMIT ",sites.count, ";")
  tables <- dbGetQuery(con, statement=tablelist_statement)
  tables <- data.frame(tables)
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
    query <- paste(" select pow,spd from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    winddata[,indx] <<- as.double(data06[,2])
  }
}

#generate.seq.matrix <- function(){
#  mat <- as.matrix(combn(sites.count, aggr.cluster.size)) #generating different combination
#  len <- length(mat[1,])
#  indxseq <- seq(1,len)#sample(1:len, indxcombicnt)
#  indxcombimat <<- mat[,indxseq]
#}

gen.aggrdata <- function(){
  col.names <- c()
  if(aggr.mat.size != 0){
    for(i in seq(1 ,aggr.mat.size)){
      indx.seq <- indxcombimat[,i]
      pmat <- as.matrix(powdata[,indx.seq])
      wmat <- as.matrix(winddata[,indx.seq])
      col.names <- c(col.names, paste('S',paste(indx.seq, collapse=""), sep=""))

      if(length(pmat[1,]) > 1){
        pow.aggrdata[,i] <<- rowSums(pmat)
        wind.aggrdata[,i] <<- rowMeans(wmat)
      }else{
        pow.aggrdata[,i] <<- pmat[,1]
        wind.aggrdata[,i] <<- wmat[,1]
      }
    }
    colnames(pow.aggrdata) <<- col.names
  }
  else{
    indx.seq <- seq(1,sites.count)
    pow.aggrdata10 <<- rowSums(as.matrix(powdata[,indx.seq]))
    wind.aggrdata10 <<- rowMeans(as.matrix(powdata[,indx.seq]))
  }
}


predict.pow <- function(aggrno, indx) {
  if(indx < 1 || indx >= data.len){
    print("Enter indx Greater than 0 and less than the data size")
    return
  }

  pow.data.normalized <- c()
  wind.data.normalized <- c()

  if(aggr.mat.size!=0){
    pow.data.normalized <- normalizeData(as.vector(pow.aggrdata[,aggrno]),type="0_1")
    wind.data.normalized <- normalizeData(as.vector(wind.aggrdata[,aggrno]),type="0_1")
  }
  else{
    pow.data.normalized <- normalizeData(as.vector(pow.aggrdata10),type="0_1")
    wind.data.normalized <- normalizeData(as.vector(wind.aggrdata10),type="0_1")
  }

  pow.normParms <<- getNormParameters(pow.data.normalized)
  wind.normParms <<- getNormParameters(wind.data.normalized)
  pow.data.set <- as.vector(pow.data.normalized[,1])
  wind.data.set <- as.vector(wind.data.normalized[,1])

  indx.start <<- indx
  indx.end <<- indx.start + window.size - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
    print(paste("Cluster size: ",aggr.cluster.size," AggrNo.: ", aggrno, " Slide No.: ", count+1))
    y <- as.vector(pow.data.set[indx.start:indx.end])
    x <- as.vector(wind.data.set[indx.start:indx.end])
    dat <- data.frame(cbind(y,x))

    train.indx <- floor(window.size *  train.data.percent)
    test.indx <- train.indx + 1
    window.slide <- window.size - train.indx
    #y.train <- y[1:train.indx]
    #x.train <- x[1:train.indx]
    #y.test <- y[test.indx:window.size]
    #x.test <- x[test.indx:window.size]

    #trn.data <- data.frame(y.train,x.train)
    #f = as.formula("y.train ~ x.train ")
    trn.data <- data.frame(dat[1:train.indx,])
    tst.x <- data.frame(x=dat$x[test.indx:window.size])
    tst.y <- data.frame(y=dat$y[test.indx:window.size])
    f = as.formula("y ~ x")

    out <<- neuralnet(f,
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
                      linear.output=TRUE #If true, act.fct is not applied to the o/p of neuron. So it will be only integartion function
    )


    data.train <<- c(data.train, trn.data$y)
    data.test <<- c(data.test, tst.y$y)

    pred <- compute(out, tst.x)$net.result
    data.out <<- c(data.out, pred)

    indx.start <<- indx.start + window.slide
    indx.end <<- indx.start + window.size
    count <- count + 1
    #if(count == 10){
    #  break
    #}
  }

  train.data <<- cbind(train.data, data.train)
  test.data <<-  cbind(test.data, data.test)
  output <<- cbind(output, data.out)
  if(aggr.cluster.size == 1){
    test.data.denorm <<- cbind(test.data.denorm, as.vector(denormalizeData(data.test,pow.normParms)))
    output.denorm <<- cbind(output.denorm, as.vector(denormalizeData(data.out, pow.normParms)))
  }
}


predict.for.combination <- function(){
  slide.indx <- data.len - (history.length * data.len.day) + 1
  #loaddata()
  #generate.seq.matrix()
  gen.aggrdata()
  if(aggr.mat.size != 0){
    for(aggr.indx in seq(1,aggr.mat.size)){
      predict.pow(aggr.indx,slide.indx)
      #break
    }
  }
  else{
    predict.pow(0,slide.indx)
  }
}

prediction.error <- function(){
  parm.count <- 5
  file.name <- paste('neuralnet_shortterm_windspeed_aggr_combi',aggr.cluster.size,'.csv',sep="")
  setwd(filepath)

  if(aggr.cluster.size == 1){
    err.data <<- matrix(0,nrow=1,ncol=parm.count, byrow=TRUE)
    colnames(err.data) <<- c("AggrNo.Seq","rmse", "mape", "sse", "mse")
    col.name <- paste('S',paste(seq(1,aggr.cluster.size), collapse=""), sep="")
    site.name <- col.name
    test <- rowSums(test.data.denorm[,1:aggr.mat.size])
    test <- normalizeData(test, type="0_1")
    pred <- rowSums(output.denorm[,1:aggr.mat.size])
    pred <- normalizeData(pred, type="0_1")
    err.rmse <- error(forecast=pred, true=test,method="rmse")
    err.mape <- error(forecast=pred, true=test,method="mape")
    err.sse <- error(forecast=pred, true=test,method="sse")
    err.mse <- error(forecast=pred, true=test,method="mse")
    err.data[1,] <<- c(site.name, err.rmse, err.mape, err.sse, err.mse)

    write.csv(err.data, file=file.name)

  }else if(aggr.cluster.size>1 && aggr.cluster.size < sites.count){
    err.data <<- matrix(0,nrow=indxcombicnt, ncol=parm.count, byrow=TRUE)
    colnames(err.data) <<- c("AggrNo.Seq","rmse", "mape", "sse", "mse")
    col.names <- colnames(pow.aggrdata)
    for(site in seq(1:(aggr.mat.size))){
      site.name <- col.names[site]
      test <- test.data[,site]
      pred <- output[,site]
      err.rmse <- error(forecast=pred, true=test,method="rmse")
      err.mape <- error(forecast=pred, true=test,method="mape")
      err.sse <- error(forecast=pred, true=test,method="sse")
      err.mse <- error(forecast=pred, true=test,method="mse")
      err.data[site,] <<- c(site.name, err.rmse, err.mape, err.sse, err.mse)
      #break
    }
    write.csv(err.data, file=file.name)
  }else{
    err.data <<- matrix(0,nrow=1,ncol=parm.count, byrow=TRUE)
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
  for(combi in seq(1,10)){#sites.count
    aggr.cluster.size <<- combi
    if(combi != sites.count){
      indxcombicnt <<-length(combn(sites.count,combi)[1,])
      aggr.mat.size <<- indxcombicnt
      pow.aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
      wind.aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
      #indxcombimat <<- ff(NA, dim=c(combi,indxcombicnt),vmode="integer")
      indxcombimat <<- as.matrix(combn(sites.count, aggr.cluster.size))
      #generate.seq.matrix()
    }else{
      indxcombicnt <<- 0
      aggr.mat.size <<- 0
      pow.aggrdata10 <<- ff(NA, dim=data.len, vmode="double")
      wind.aggrdata10 <<- ff(NA, dim=data.len, vmode="double")
    }

    train.data <<- c()
    test.data <<- c()
    test.data.denorm <<- c()
    output <<- c()
    output.denorm <<- c()

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
