require(neuralnet)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)

row.size <- 52560
row.size.monthly <<- 4380
mat.size <<- 12
powdata <- ff(NA, dim=row.size, vmode="double")
train.data.percent <- 0.7
query <- "select t.pow from onshore_SITE_07829 t where (t.mesdt >= 20060101 && t.mesdt < 20070101);"
drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")
powdata <<- data.frame(dbGetQuery(con, statement=query), check.names=FALSE)


predict <- function(){
  d06 <- as.vector(c(unlist(powdata)))
  d06.monthly <- matrix(d06, nrow=row.size.monthly, ncol=mat.size, byrow=FALSE)
  colnames(d06.monthly) <- paste("mon",c(1:12), sep="")
  train.dataset.indx <- row.size.monthly *  train.data.percent
  test.dataset.indx <- train.dataset.indx + 1
  data.mat.train <- data.frame(d06.monthly[1:train.dataset.indx, 1:11])
  data.mat.test <- data.frame(d06.monthly[test.dataset.indx:row.size.monthly, 1:11])

  out <<- neuralnet(mon11 ~ mon1 + mon2 + mon3 + mon4 + mon5 + mon6 + mon7 + mon8 + mon9 + mon10,
                    data.mat.train,
                    threshold=0.01,
                    hidden=10,
                    err.fct="sse",
                    linear.output=FALSE
  )
  plot.nn(out)
  test <<- as.data.frame(as.vector(d06.monthly[,"mon11"]))
  pred <<- compute(out, data.mat.test[,1:10])$net.result
}
