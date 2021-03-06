require(gbm)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(RSNNS)
require(ftsa)

sites.count <- 10
data.len <- 52560
data.len.day <<- 144
mat.size <<- 365
history.length <- 50
hidden.nodes <<- 10#c(round(window.size/2), window.size,1)
window.size <- 1440
train.data.percent <- 0.7
file.name <- "gbm_shortterm_windspeed_simple.csv"
file.path <- "/home/freak/Programming/Thesis/results/results/gbm_shortterm_windspeed_simple/"
table.ip.type <- "specific"#"random"
powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
winddata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
winddata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")

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
    query <- paste(" select pow,spd from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    winddata[,indx] <<- as.double(data06[,2])
    powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]),type="0_1")
    winddata.normalized[,indx] <<- normalizeData(as.vector(data06[,2]),type="0_1")
  }
}


predict.pow <- function(siteno, indx) {
  if(indx < 1 || indx >= data.len){
    print("Enter indx Greater than 0 and less than the data size")
    return
  }

  data.set.pow <<- as.vector(powdata.normalized[,siteno])
  data.set.wind <<- as.vector(winddata.normalized[,siteno])

  indx.start <<- indx
  indx.end <<- indx + (window.size - 1)
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
    print(paste("Site no: ", siteno, "Slide Count: ", count+1))
    y <- as.vector(data.set.pow[indx.start:indx.end])
    x <- as.vector(data.set.wind[indx.start:indx.end])
    dat <- data.frame(cbind(y,x))

    train.indx <- floor(window.size *  train.data.percent)
    test.indx <- train.indx + 1
    window.slide <- window.size - train.indx
    #y.train <- y[1:train.indx]
    #x.train <- x[1:train.indx]
    #y.test <- y[test.indx:window.size]
    #x.test <- data.frame(x[test.indx:window.size])

    #trn.data <- data.frame(y.train,x.train)
    trn.data <- data.frame(dat[1:train.indx,])
    tst.x <- data.frame(x = dat$x[test.indx:window.size])
    tst.y <- data.frame(y=dat$y[test.indx:window.size])
    f = as.formula("y ~ x")

    out <<- gbm(f,
                data=trn.data,
                distribution ="gaussian",
                n.trees=10000,
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
x1 = train.data[,1]
x2 = test.data[,1]
y = output[,1]

length(y)
plot(y, type='l')
plot(x1, type='l')
plot(x2, type='l')

dataToPlot = data.frame(seq(1,4330),x2,y)
Line <- gvisLineChart(dataToPlot)
plot(Line)

powpred <- denormalizeData(y,getNormParameters(powdata.normalized))
plot(powpred, type='l')
summary(out)
out$generalized.weights


#Error Calculation
install.packages('ftsa')
require('ftsa')
err <- error(forecast=y, true=x2,method="mape")
err
