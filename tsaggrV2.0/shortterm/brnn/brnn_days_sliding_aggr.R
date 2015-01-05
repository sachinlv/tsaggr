require(brnn)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(combinat)

sites.count <- 10
indxcombicnt <- 10 #no. of combination
aggr.mat.size <- sites.count/2
data.len <- 52560
data.len.day <<- 144
mat.size <<- 365
window.size <- 10
train.data.percent <- 0.7
indxseq <- c(seq(1,sites.count))
slide.count <- mat.size-window.size+1

powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
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
          query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
          data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
	      powdata[,indx] <<- as.double(data06[,1])
	    }
}

generate.seq.matrix <- function(){
	  mat <- as.matrix(combn(sites.count, aggr.mat.size)) #generating different combination
  len <- length(mat[1,])
    seq <- sample(1:len, indxcombicnt)
    indxcombimat <<- mat[,seq]
}

gen.aggrdata <- function(seq1, seq2){
  seq1 <- as.vector(seq1)
  seq2 <- as.vector(seq2)
    for(i in seq(1 ,aggr.mat.size)){
      indx1 <- seq1[i]
      indx2 <- seq2[i]
      aggrdata[,i] <<- c(powdata[,indx1]) + c(powdata[,indx2])
    }
}


predict <- function(siteno, indx) {
  if(indx < 1 || indx >= data.len){
    print("Enter indx Greater than 0 and less than the data size")
    return
  }

  #data.normalized <- normalizeData(as.vector(aggrdata[,siteno]),type="0_1")
  #data.set <- as.vector(data.normalized[,1])
  data.set <- as.vector(aggrdata[,siteno])

  indx.start <<- indx
  indx.end <<- indx.start + (window.size * data.len.day) - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  #for(indx in seq(1 ,mat.size-window.size+1)){
  while(indx.end <= data.len){
      print(paste("Slide Count: ", count+1))
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
                    mu=0.005,
                    mu_dec=0.1,
                    mu_max=1e10,
                    change = 0.001,
                    neurons=window.size,
                    normalize=TRUE,
                    verbose=TRUE,
                    Monte_Carlo = FALSE)

      data.train <<- c(data.train, data.mat.train[,window.size])
      data.test <<- c(data.test, data.mat.test[,window.size])

      pred <- predict.brnn(out ,data.mat.test[,1:window.size-1])
      data.out <<- c(data.out,as.double(pred))

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
  write.csv(err.data, file=paste(file.path,file.name))
}

predict_for_combination()
generate.seq.matrix()


#plotting
x1 = train.data
x2 = test.data
y = output

plot(y, type='l')
plot(x1, type='l')
plot(x2, type='l')
length(x2)
dataToPlot = data.frame(seq(1,430),x2,y)
Line <- gvisLineChart(dataToPlot)
plot(Line)


#error measure
err <- error(forecast=y, true=x2,method="mape")
err
