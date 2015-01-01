require(brnn)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)

row.size <- 52560
row.size.monthly <<- 4380
mat.size <<- 12
window.size <- 4
train.data.percent <- 0.7
powdata <<- ff(NA, dim=row.size, vmode="double")
train.data <<- ff(NA, dim=c(row.size.monthly*train.data.percent, mat.size-window.size+1), vmode="double")
test.data <<- ff(NA, dim=c(row.size.monthly*(1-train.data.percent), mat.size-window.size+1), vmode="double")
output <<- ff(NA, dim=c(row.size.monthly*(1-train.data.percent), mat.size-window.size+1), vmode="double")

query <- "select t.pow from onshore_SITE_07829 t where (t.mesdt >= 20060101 && t.mesdt < 20070101);"
drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")
powdata <<- data.frame(dbGetQuery(con, statement=query), check.names=FALSE)

data.set <- matrix(powdata[,1], nrow=row.size.monthly, ncol=mat.size, byrow=FALSE)

predict <- function() {
  for(indx in seq(1 ,mat.size-window.size+1)){
    data.set.indx <- c(seq(indx, indx+window.size-1))
    data.mat <- as.matrix(data.set[,data.set.indx])
    colnames(data.mat) <- paste("d",c(1:window.size), sep="")

    train.dataset.indx <- row.size.monthly *  train.data.percent
    test.dataset.indx <- train.dataset.indx + 1
    data.mat.train <- data.mat[1:train.dataset.indx,]
    data.mat.test <- data.mat[test.dataset.indx:row.size.monthly,]
    out <<- brnn(d4 ~ d1 + d2 + d3,
                 data.mat.train,
                 epochs=1000,
                 cores=2,
                 mu=0.005,
                 mu_dec=0.1,
                 mu_max=1e10,
                 change = 0.001,
                 neurons=24,
                 normalize=FALSE,
                 verbose=TRUE,
                 Monte_Carlo = FALSE)

    train.data[, indx] <<- data.mat.train[,window.size]
    test.data[, indx] <<- data.mat.test[,window.size]
    pred <- predict.brnn(out ,data.mat.test[, 1:window.size-1])
    output[,indx] <<- as.double(pred)
    plot(pred, type='l')
    break
  }
}

predict()

x1 = train.data[1:1314,1]
x2 = test.data[,1]
y = output[,1]
length(y)
plot(y, type='l')
plot(x1, type='l')
plot(x2, type='l')

dataToPlot = data.frame(seq(1,1314),x2,y)
Line <- gvisLineChart(dataToPlot)
plot(Line)


n = normalize.vector(x)

