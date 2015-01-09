require(neuralnet)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(combinat)
require(RSNNS)
require(ftsa)

sites.count <- 10
indxcombicnt <- 10 #no. of combination
aggr.mat.size <- sites.count/2
hidden.nodes <<- 10#c(round(window.size/2), window.size,1)
data.len <- 52560
#data.len.day <<- 144
#mat.size <<- 365
window.size <- 1440
train.data.percent <- 0.7
indxseq <- c(seq(1,sites.count))
#slide.count <- mat.size-window.size+1

powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
winddata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
pow.aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
wind.aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")

train.data <<- c()
test.data <<- c()
output <<- c()
indxcombimat <<- ff(NA, dim=c(aggr.mat.size,indxcombicnt),vmode="integer")

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
    query <- paste(" select pow, spd from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    winddata[,indx] <<- as.double(data06[,2])
  }
}

generate.seq.matrix <- function(){
  mat <- as.matrix(combn(sites.count, aggr.mat.size)) #generating different combination
  len <- length(mat[1,])
  indxseq <- sample(1:len, indxcombicnt)
  indxcombimat <<- mat[,indxseq]
}

gen.aggrdata <- function(seq1, seq2){
  seq1 <- as.vector(seq1)
  seq2 <- as.vector(seq2)
  for(i in seq(1 ,aggr.mat.size)){
    indx1 <- seq1[i]
    indx2 <- seq2[i]
    pow.aggrdata[,i] <<- c(powdata[,indx1]) + c(powdata[,indx2])
    wind.aggrdata[,i] <<- rowMeans(cbind(winddata[,indx1], winddata[,indx2]))
  }
}


predict <- function(siteno, indx) {
  if(indx < 1 || indx >= data.len){
    print("Enter indx Greater than 0 and less than the data size")
    return
  }

  pow.data.normalized <- normalizeData(as.vector(pow.aggrdata[,siteno]),type="0_1")
  wind.data.normalized <- normalizeData(as.vector(wind.aggrdata[,siteno]),type="0_1")
  #data.set <- matrix(normalized.data, nrow=data.len.day, ncol=mat.size, byrow=FALSE)
  pow.data.set <- as.vector(pow.data.normalized[,1])
  wind.data.set <- as.vector(wind.data.normalized[,1])

  indx.start <<- indx
  indx.end <<- indx.start + window.size - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
    y <- as.vector(pow.data.set[indx.start:indx.end])
    x <- as.vector(wind.data.set[indx.start:indx.end])

    train.indx <- floor(window.size *  train.data.percent)
    test.indx <- train.indx + 1
    window.slide <- window.size - train.indx
    y.train <- y[1:train.indx]
    x.train <- x[1:train.indx]
    y.test <- y[test.indx:window.size]
    x.test <- x[test.indx:window.size]

    train.data <- data.frame(y.train,x.train)
    f = as.formula("y.train ~ x.train ")

    out <<- neuralnet(f,
                      train.data,
                      hidden=hidden.nodes,
                      rep=1,
                      stepmax = 2000,
                      threshold=0.1,
                      learningrate=1,
                      algorithm="rprop+", #'rprop-', 'sag', 'slr'
                      startweights=NULL,
                      lifesign="none",
                      err.fct="sse",
                      act.fct="logistic",
                      #exclude = NULL,
                      #constant.weights = NULL,
                      linear.output=TRUE #If true, act.fct is not applied to the o/p of neuron. So it will be only integartion function
    )


    data.train <<- c(data.train, y.train)
    data.test <<- c(data.test, y.test)

    pred <- compute(out, x.test)$net.result
    data.out <<- c(data.out, pred)

    indx.start <<- indx.start + window.slide
    indx.end <<- indx.start + window.size
    count <- count + 1
    if(count == 10){
      break
    }
  }

  train.data <<- cbind(train.data, data.train) #data.mat.train[,window.size]
  test.data <<-  cbind(test.data, data.test) #data.mat.test[,window.size]
  output <<- cbind(output, data.out)
}


predict_for_combination <- function(){
  slide.indx <- 1
  loaddata()
  generate.seq.matrix()

  for(i in seq(1:indxcombicnt)){
    mat.indx1 <- as.vector(indxcombimat[,i])
    mat.indx2 <- as.vector(indxseq[is.na(pmatch(indxseq,mat.indx1))])
    gen.aggrdata(mat.indx1, mat.indx2)

    for(aggr.indx in seq(1,aggr.mat.size)){
      predict(aggr.indx,slide.indx)
      break
    }
    break
  }
}

prediction.error <- function(){
  parm.count <- 4
  err.data <<- matrix(,nrow=sites.count, ncol=parm.count, byrow=TRUE)
  #colnames(err.data) <<- c("site.id", "rmse", "mape", "sse", "mse")
  colnames(err.data) <<- c("rmse", "mape", "sse", "mse")

  for(site in seq(1:(indxcombicnt*aggr.mat.size))){
    #site.name <- tables[site,]
    test <- test.data[,site]
    pred <- output[,site]
    err.rmse <- error(forecast=pred, true=test,method="rmse")
    err.mape <- error(forecast=pred, true=test,method="mape")
    err.sse <- error(forecast=pred, true=test,method="sse")
    err.mse <- error(forecast=pred, true=test,method="mse")
    #err.data[site,] <<- c(site.name, err.rmse, err.mape, err.sse, err.mse)
    err.data[site,] <<- c(err.rmse, err.mape, err.sse, err.mse)
    break
  }
  #write.csv(err.data, file=paste(file.path,file.name))
}

predict_for_combination()
prediction.error()

err.data
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
