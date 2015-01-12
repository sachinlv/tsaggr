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
window.size <- 10
train.data.percent <- 0.7
file.name <- "gbm_shortterm_simple.csv"
file.path <- "/home/freak/Programming/Thesis/results/results/gbm_shortterm_simple/"

powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
train.data <<- c()
test.data <<- c()
output <<- c()
#train.data <<- ff(NA, dim=c(data.len.day*train.data.percent, mat.size-window.size+1), vmode="double")
#test.data <<- ff(NA, dim=c(data.len.day*(1-train.data.percent), mat.size-window.size+1), vmode="double")
#output <<- ff(NA, dim=c(data.len.day*(1-train.data.percent), mat.size-window.size+1), vmode="double")

drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")
tablelist_statement = paste("SELECT TABLE_NAME FROM information_schema.TABLES ",
                            "WHERE TABLE_SCHEMA = 'eastwind' AND",
                            "TABLE_NAME LIKE 'onshore_SITE_%' "," LIMIT ",sites.count, ";")
tables <<- dbGetQuery(con, statement=tablelist_statement)
tables <<- data.frame(tables)


loaddata <- function(){
  for(indx in seq(1,sites.count)){
    tab <- tables[indx,]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]),type="0_1")
  }
}


predict <- function(siteno, indx) {
  if(indx < 1 || indx >= data.len){
    print("Enter indx Greater than 0 and less than the data size")
    return
  }

  data.set <<- as.vector(powdata.normalized[,siteno])

  indx.start <<- indx
  indx.end <<- indx + (window.size * data.len.day) - 1
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
    data.mat.train <- data.frame(data.mat[1:train.dataset.indx,])
    data.mat.test <- data.frame(data.mat[test.dataset.indx:data.len.day,])

    formula.set <- colnames(data.mat)
    y = formula.set[window.size]
    x = formula.set[1:window.size-1]
    f = as.formula(paste(y, " ~ ", paste(x, collapse="+")))

    out <<- gbm(f,
                data=data.mat.train,
                distribution ="gaussian",
                n.trees=10000,
                interaction.depth = 10,
                n.minobsinnode = 5,
                shrinkage =  0.008,
                n.cores=3)

    data.train <<- c(data.train, data.mat.train[,window.size])
    data.test <<- c(data.test, data.mat.test[,window.size])

    pred <- predict.gbm(out, data.mat.test[,1:window.size-1], out$n.trees)
    data.out <<- c(data.out, pred)

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

predict.power <- function(){
  slide.indx <- 1
  loaddata()

  for(site in seq(1:sites.count)){
    predict(site,slide.indx)
    break
  }
}

prediction.error <- function(){
  parm.count <- 5
  err.data <<- matrix(,nrow=sites.count, ncol=parm.count, byrow=TRUE)
  colnames(err.data) <<- c("site.id", "rmse", "mape", "sse", "mse")
  setwd(file.path)
  for(site in seq(1:sites.count)){
    site.name <- tables[site,]
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

dataToPlot = data.frame(seq(1,440),x2,y)
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


#formula creation
len = 10
form = paste('d',c(1:len), sep="")
y = form[len]
x = form[1:len-1]
f = as.formula(paste(y, " ~ ", paste(x, collapse="+")))
f


m = matrix(1:12,nrow=3, ncol=4,byrow = FALSE)
colnames(m) = paste('d',c(1:4),sep='')
x = colnames(m)



#Normalizazion testing
require("ppls")
a = c(1,2,3,4,5,6,7,8,9,100,101,102)
tmp = matrix(a,nrow=4,ncol=3)

b = normalize.vector(a)
plot(b, type="l")
plot(a, type="l")

n = normalize.vector(x)

#Normalize and Denormalize
install.packages("RSNNS")
require("RSNNS")
a = c(1,2,3,4,5,6,7,8,9,100,101,102)
tmp = matrix(a,nrow=4,ncol=3)
norm <- normalizeData(a, type="0_1")
norm[,1]
denorm <- denormalizeData(norm[6:12,1], getNormParameters(norm))
denorm[,1]
plot(norm[,1], type="l")
plot(denorm[,1], type="l")
