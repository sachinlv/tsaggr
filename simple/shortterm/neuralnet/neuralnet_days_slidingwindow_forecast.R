require(neuralnet)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(RSNNS)
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
#powdata.normalized <<- normalize.vector(as.vector(powdata[,1]))
powdata.normalized <<- normalizeData(as.vector(powdata[,1]),type="0_1")
#data.set <- matrix(powdata[,1], nrow=row.size.day, ncol=mat.size, byrow=FALSE)
data.set <<- matrix(powdata.normalized[,1], nrow=row.size.day, ncol=mat.size, byrow=FALSE)
#summary(data.set)

predict <- function() {
  for(indx in seq(1 ,mat.size-window.size+1)){
    print(paste("window no. : ", indx))
    data.set.indx <- c(seq(indx, indx+window.size-1))
    data.mat <- as.matrix(data.set[,data.set.indx])
    colnames(data.mat) <- paste("d",c(1:window.size), sep="")
    formula.set <- colnames(data.mat)
    train.dataset.indx <- row.size.day *  train.data.percent
    test.dataset.indx <- train.dataset.indx + 1
    data.mat.train <- data.mat[1:train.dataset.indx,]
    data.mat.test <- data.mat[test.dataset.indx:row.size.day,]
    y = formula.set[window.size]
    x = formula.set[1:window.size-1]
    f = as.formula(paste(y, " ~ ", paste(x, collapse="+")))
    #print(f)
    out <<- neuralnet(f,
                      data.mat.train,
                      hidden=c(round(window.size/2), window.size,1),#window.size
                      rep=10,
                      stepmax = 2000,
                      threshold=0.01,
                      learningrate=1,
                      algorithm="rprop+", #'rprop-', 'sag', 'slr', 'rprop+'
                      #lifesign="none",
                      err.fct="sse",
                      act.fct="logistic",
                      #exclude = NULL,
                      #constant.weights = NULL,
                      linear.output=FALSE #If true, act.fct is not applied to the o/p of neuron. So it will be only integartion function
                    )

    train.data[, indx] <<- data.mat.train[,window.size]
    test.data[, indx] <<- data.mat.test[,window.size]

    pred <- compute(out, data.mat.test[,1:window.size-1])$net.result
    output[,indx] <<- as.double(pred)

    if(indx == 10){
      break
    }

    plot.nn(out)

  }
}

predict()

plot(test.data[,1], type='l')
#plotting
head(train.data[,1:11])
x = data.set[,window.size]
x1 = as.vector(train.data[,1:10])
x2 = as.vector(test.data[,1:10])
y = as.vector(output[,1:10])

#plotting
length(y)
plot(y, type='l')
plot(x1, type='l')
plot(x2, type='l')
length(y)

dataToPlot = data.frame(seq(1,430),x2,y)
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
