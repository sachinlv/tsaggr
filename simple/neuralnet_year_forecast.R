require(neuralnet)
require(RMySQL)
require(ff)

row.size <- 52560
powdata <- ff(NA, dim=row.size, vmode="double")

query1 <- "select t.pow from onshore_SITE_07829 t where (t.mesdt >= 20050101 && t.mesdt < 20060101);"
query2 <- "select t.pow from onshore_SITE_07829 t where (t.mesdt >= 20060101 && t.mesdt < 20070101);"

drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")
powdata05 <- data.frame(dbGetQuery(con, statement=query1), check.names=FALSE)
powdata06 <- data.frame(dbGetQuery(con, statement=query2), check.names=FALSE)

data.year <- as.matrix(cbind(powdata05, powdata06))
colnames(data.year) <- c("c05", "c06")
data.mat <- data.frame(data.year)

tmp = as.vector(unlist(powdata05))
mean(tmp)

out <<- neuralnet(c06 ~ c05,
                  data.mat,
                  hidden=10,
                  rep=5,
                  stepmax = 1000,
                  threshold=0.01,
                  learningrate=1,
                  algorithm="rprop+",
                  lifesign="none",
                  err.fct="sse",
                  exclude = NULL,
                  constant.weights = NULL,
                  linear.output=TRUE)
plot.nn(out)

pred <- compute(out, data.mat[,"c06"])
x = as.vector(unlist(data.mat[,"c06"]))[50001:52560]
y = pred$net.result
plot(y, type="l")

debug(neuralnet)
