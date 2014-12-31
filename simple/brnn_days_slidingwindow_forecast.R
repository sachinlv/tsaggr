require(brnn)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(ftsa)

row.size <- 52560
row.size.day <<- 144
mat.size <<- 365
window.size <- 10
train.data.percent <- 0.7
powdata <<- ff(NA, dim=row.size, vmode="double")
train.data <<- ff(NA, dim=c(row.size.day*train.data.percent, mat.size-window.size+1), vmode="double")
test.data <<- ff(NA, dim=c(row.size.day*(1-train.data.percent), mat.size-window.size+1), vmode="double")
output <<- ff(NA, dim=c(row.size.day*(1-train.data.percent), mat.size-window.size+1), vmode="double")

query <- "select t.pow from onshore_SITE_07829 t where (t.mesdt >= 20060101 && t.mesdt < 20070101);"
drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")
powdata <<- data.frame(dbGetQuery(con, statement=query), check.names=FALSE)

data.set <- matrix(powdata[,1], nrow=row.size.day, ncol=mat.size, byrow=FALSE)

predict <- function() {
  for(indx in seq(1 ,mat.size-window.size+1)){
    print(paste("window no. :", indx))
    data.set.indx <- c(seq(indx, indx+window.size-1))
    data.mat <- as.matrix(data.set[,data.set.indx])
    colnames(data.mat) <- paste("d",c(1:window.size), sep="")
    formula.set <- colnames(data.mat)
    train.dataset.indx <- row.size.day *  train.data.percent
    test.dataset.indx <- train.dataset.indx + 1
    data.mat.train <- data.mat[1:train.dataset.indx,]
    data.mat.test <- data.mat[test.dataset.indx:row.size.day,]
    #y = as.vector(data.mat.train[,window.size])
    #x = as.matrix(data.mat.train[,1:window.size-1])
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

    train.data[, indx] <<- data.mat.train[,window.size]
    test.data[, indx] <<- data.mat.test[,window.size]

    pred <- predict.brnn(out ,data.mat.test[, 1:window.size-1])
    output[,indx] <<- as.double(pred)

    #plot(pred, type='l')
    if(indx == 10){
      break
    }
  }
}

predict()


x1 = as.vector(train.data[,1:10])
x2 = as.vector(test.data[,1:10])
y = as.vector(output[,1:10])
length(x2)
plot(y, type='l')
plot(x1, type='l')
plot(x2, type='l')

dataToPlot = data.frame(seq(1,430),x2,y)
Line <- gvisLineChart(dataToPlot)
plot(Line)

#error measure
err <- error(forecast=y, true=x2,method="mape")
err
rmse(x2,y)
mean(abs((x2-y)/x2) * 100)
mean(y)

