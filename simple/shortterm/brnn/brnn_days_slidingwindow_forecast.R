require(brnn)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(RSNNS)
require(ftsa)

sites.count <- 10
history.length <- 50
data.len <- 52560
data.len.day <<- 144
mat.size <<- 365
window.size <- 10
train.data.percent <- 0.7
file.name <- "brnn_shortterm_simple.csv"
file.path <- "/home/freak/Programming/Thesis/results/results/brnn_shortterm_simple/"
table.ip.type <- "specific"#"random"
powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
train.data <<- c()
test.data <<- c()
output <<- c()

drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")

if(table.ip.type == "random"){
    tablelist_statement = paste("SELECT TABLE_NAME FROM information_schema.TABLES ",
                                "WHERE TABLE_SCHEMA = 'eastwind' AND",
                                "TABLE_NAME LIKE 'onshore_SITE_%' "," LIMIT ",sites.count, ";")
    tables <<- dbGetQuery(con, statement=tablelist_statement)
    tables <<- data.frame(tables)
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
  for(indx in seq(1,sites.count)){
    tab <- tables[indx,1]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]),type="0_1")
  }
}


predict.pow <- function(siteno,indx) {
  if(indx < 1 || indx >= data.len){
    print("Enter indx Greater than 0 and less than the data size")
    return
  }

  data.set <<- as.vector(powdata[,siteno])#powdata.normalized

  indx.start <<- indx
  indx.end <<- indx.start + (window.size * data.len.day) - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
    print(paste("Site no: ", siteno, "Slide Count: ", count+1))
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
                 normalize=FALSE,
                 verbose=FALSE,
                 Monte_Carlo = FALSE)

    data.train <<- c(data.train, data.mat.train[,window.size])
    data.test <<- c(data.test, data.mat.test[,window.size])

    pred <- predict.brnn(out ,data.mat.test[, 1:window.size-1])
    data.out <<- c(data.out, pred)

    indx.start <<- indx.start + window.slide
    indx.end <<- indx.start + (window.size * data.len.day) -1
    count <- count + 1
    if(count == 10){
      break
    }
  }

  train.data <<- cbind(train.data, data.train)
  test.data <<-  cbind(test.data, data.test)
  output <<- cbind(output, data.out)
}

predict.power <- function(){
  slide.indx <- data.len - (history.length * data.len.day) + 1
  loaddata()

  for(site in seq(1:sites.count)){
    predict.pow(site,slide.indx)
    break
  }
}

prediction.error <- function(){
  parm.count <- 5
  err.data <<- matrix(,nrow=sites.count, ncol=parm.count, byrow=TRUE)
  colnames(err.data) <<- c("site.id", "rmse", "mape", "sse", "mse")
  setwd(file.path)
  for(site in seq(1:sites.count)){
    site.name <- as.character(tables[site,1])
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

predict.power()
prediction.error()

err.data

#plotting
x1 = data.train[,1]
x2 = data.test[,1]
y = data.out[,1]
length(x2)

plot(y, type='l')
plot(x1, type='l')
plot(x2, type='l')

dataToPlot = data.frame(seq(1,440),x2,y)
Line <- gvisLineChart(dataToPlot)
plot(Line)

#error measure
err <- error(forecast=y, true=x2,method="mape")
err
rmse(x2,y)
mean(abs((x2-y)/x2) * 100)
mean(y)

