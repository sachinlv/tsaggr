require(nnet)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)


row.size <- 52560
row.size.monthly <<- 4380
mat.size <<- 12
powdata <- ff(NA, dim=row.size, vmode="double")

query <- "select t.pow from onshore_SITE_07829 t where (t.mesdt >= 20060101 && t.mesdt < 20070101);"
drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")
powdata <<- data.frame(dbGetQuery(con, statement=query), check.names=FALSE)

predictnn <- function(){
  d06 <- as.vector(c(unlist(powdata)))
  d06.monthly <- matrix(d06, nrow=row.size.monthly, ncol=mat.size, byrow=FALSE)
  colnames(d06.monthly) <- paste("mon",c(1:12), sep="")
  data.mat <<- as.matrix(d06.monthly[,1:10])
  data.train <<- d06.monthly[,11]
  data.test <<- d06.monthly[,12]

  out <<- nnet( data.mat,
                data.train,
                linout=FALSE,
                size = 5)

  pred <<- predict(out, data.test)#,type="raw")
}

predictnn()

head(data.train)
head(data.test)
plot(pred, type="l")


#example1
data=seq(0,2,by=0.01)
x.train = data[1:150]
x.test = data[151:201]
y=2*sin(2*pi*(x.train-1/4))
plot(y, type="l")
data.mat <- data.frame(cbind(x.train,y))
o <- nnet(x.train~., data.mat, size=2)
summary(o)
p <- predict(o, x.test)
o

plot(p, type="l")




#example2
# load data
data(longley)
head(longley)
x <- longley[,1:6]
y <- longley[,7]
# fit model
fit <- nnet(Employed~., longley, size=12, maxit=500, linout=T, decay=0.01)
# summarize the fit
summary(fit)
# make predictions
predictions <- predict(fit, x, type="raw")
# summarize accuracy
rmse <- mean((y - predictions)^2)
print(rmse)
