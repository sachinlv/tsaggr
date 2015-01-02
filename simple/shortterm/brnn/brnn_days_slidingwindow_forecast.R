require(brnn)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(ftsa)

data.len <- 52560
data.len.day <<- 144
mat.size <<- 365
window.size <- 10
train.data.percent <- 0.7
powdata <<- ff(NA, dim=data.len, vmode="double")
#train.data <<- ff(NA, dim=c(data.len.day*train.data.percent, mat.size-window.size+1), vmode="double")
#test.data <<- ff(NA, dim=c(data.len.day*(1-train.data.percent), mat.size-window.size+1), vmode="double")
#output <<- ff(NA, dim=c(data.len.day*(1-train.data.percent), mat.size-window.size+1), vmode="double")

query <- "select t.pow from onshore_SITE_07829 t where (t.mesdt >= 20060101 && t.mesdt < 20070101);"
drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")
powdata <<- data.frame(dbGetQuery(con, statement=query), check.names=FALSE)
powdata.normalized <<- normalizeData(as.vector(powdata[,1]),type="0_1")
#data.set <<- matrix(powdata.normalized[,1], nrow=data.len.day, ncol=mat.size, byrow=FALSE)
data.set <<- as.vector(powdata.normalized[,1])

predict <- function(indx) {
  if(indx < 1 || indx >= data.len){
    print("Enter indx Greater than 0 and less than the data size")
    return
  }

  indx.start <<- indx
  indx.end <<- indx.start + (window.size * data.len.day) - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
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

    out <<- brnn(f,#x,y,
                 data.mat.train,
                 epochs=1000,
                 cores=2,
                 mu=0.005,
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
    indx.end <<- indx.start + (window.size * data.len.day) -1
    count <- count + 1

    if(count == 10){
      break
    }
  }
}

predict(1)

#plotting
x1 = data.train
x2 = data.test
y = data.out
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

