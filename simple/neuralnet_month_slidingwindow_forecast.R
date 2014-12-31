require(neuralnet)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)

row.size <- 52560
row.size.monthly <<- 4380
mat.size <<- 12
window.size <- 4
train.data.percent <- 0.7
powdata <<- ff(NA, dim=row.size, vmode="double")
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

      out <<- neuralnet(d4 ~ d1 + d2 + d3,
                        data.mat.train,
                        hidden=window.size,
                        rep=10,
                        stepmax = 5000,
                        threshold=0.01,
                        learningrate=1,
                        algorithm="rprop+",
                        lifesign="none",
                        err.fct="sse",
                        exclude = NULL,
                        constant.weights = NULL,
                        linear.output=TRUE)
      plot.nn(out)
      test.data[, indx] <<- data.mat.test[,window.size]
      pred <- compute(out, data.mat.test[,1:window.size-1])$net.result

      output[,indx] <<- as.double(pred)
      break
    }
}

predict()

x = test.data[,1]
y = output[,1]
length(x)
plot(y, type='l')
plot(x, type='l')

mean(y)

mean(data.set[(0.7*row.size.monthly+1):row.size.monthly,1])

