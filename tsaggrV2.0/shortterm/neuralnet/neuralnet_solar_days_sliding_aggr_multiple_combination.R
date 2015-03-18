require(neuralnet)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(combinat)
require(RSNNS)
require(ftsa)
require(zoo)
require(matrixStats)
require(forecast)
require(timeSeries)

sites.count <- 10
history.length <- 365
data.len.day <<- 24
data.len <- history.length * data.len.day
#window.size <- 24
fcstart <- 355
horizon <- 12
start.date <- '20060101'
end.date <- '20061231'
hidden.nodes <<- 10#c(round(window.size/2), window.size,1)

filepath.generic <<- '/home/freak/Programming/Thesis/results/results/neuralnet_shortterm_solar_'
file.name.generic <<- 'neuralnet_shortterm_solar_aggr_combi'
file.name.denorm.generic <<- 'neuralnet_shortterm_solar_aggr_combi_denorm'
file.name.aggr.generic <<- 'neuralnet_shortterm_solar_aggr_combi_aggr'
file.name.aggr.denorm.generic <<- 'neuralnet_shortterm_solar_aggr_combi_aggr_denorm'

table.ip.type <- "random"#"specific"
aggr.type.vec <<- c('aggr')

powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
solardata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")


drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="solar",user="sachin",pass="password")

if(table.ip.type == "random"){
  t <-   c("solar_random_site10",
           "solar_random_site11",
           "solar_random_site12",
           "solar_random_site13",
           "solar_random_site14",
           "solar_random_site15",
           "solar_random_site16",
           "solar_random_site17",
           "solar_random_site18",
           "solar_random_site19")
  tables <<- data.frame(cbind(numeric(0),t))
}else{
  t <-   c("solar_specific_site10",
           "solar_specific_site11",
           "solar_specific_site12",
           "solar_specific_site13",
           "solar_specific_site14",
           "solar_specific_site15",
           "solar_specific_site16",
           "solar_specific_site17",
           "solar_specific_site18",
           "solar_specific_site19")

  tables <<- data.frame(cbind(numeric(0),t))
}

loaddata <- function(){
  colindx <- 1
  for(indx in seq(1,sites.count)){
    tab <- tables[indx,1]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow,rad from ", tab, " WHERE (mesdt >= ",start.date," && mesdt <= ",end.date,") LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    solardata[,indx] <<- as.double(data06[,2])
  }
}

gen.aggrdata <- function(){
  col.names <- c()
  if(aggr.mat.size != 0){
    for(i in seq(1 ,aggr.mat.size)){
      indx.seq <- indxcombimat[,i]
      pmat <- as.matrix(powdata[,indx.seq])
      smat <- as.matrix(solardata[,indx.seq])
      col.names <- c(col.names, paste('S',paste(indx.seq, collapse=""), sep=""))

      if(length(pmat[1,]) > 1){
        pow.aggrdata[,i] <<- rowSums(pmat)
        solar.aggrdata[,i] <<- rowWeightedMeans(smat)#rowMeans(wmat)
      }else{
        pow.aggrdata[,i] <<- pmat[,1]
        solar.aggrdata[,i] <<- smat[,1]
      }
    }
    colnames(pow.aggrdata) <<- col.names
  }
  else{
    indx.seq <- seq(1,sites.count)
    pow.aggrdata10 <<- rowSums(as.matrix(powdata[,indx.seq]))
    solar.aggrdata10 <<- rowWeightedMeans(as.matrix(solardata[,indx.seq]))
  }
}


predict.pow <- function(aggrno) {
  pow.data.normalized <- c()
  solar.data.normalized <- c()
  if(aggr.mat.size!=0){
    pow.data.normalized <- normalizeData(as.vector(pow.aggrdata[,aggrno]),type="0_1")
    solar.data.normalized <- normalizeData(as.vector(solar.aggrdata[,aggrno]),type="0_1")
  }
  else{
    pow.data.normalized <- normalizeData(as.vector(pow.aggrdata10),type="0_1")
    solar.data.normalized <- normalizeData(as.vector(solar.aggrdata10),type="0_1")
  }

  pow.normParms <<- getNormParameters(pow.data.normalized)
  solar.normParms <<- getNormParameters(solar.data.normalized)
  pow.data.set <- as.vector(pow.data.normalized[,1])
  solar.data.set <- as.vector(solar.data.normalized[,1])

  indx.start <<- 1
  indx.end <<- fcstart * data.len.day#indx.start + window.size - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
    print(paste("Cluster size: ",aggr.cluster.size," AggrNo.: ", aggrno, " Slide No.: ", count+1))
    y <- as.vector(pow.data.set[indx.start:indx.end])
    x <- as.vector(solar.data.set[indx.start:indx.end])
    dat <- data.frame(cbind(y,x))
    current.datalen <- length(y)

    train.indx <-current.datalen - horizon
    test.indx <- train.indx + 1

    trn.data <- data.frame(dat[1:train.indx,])
    tst.x <- data.frame(x=dat$x[test.indx:current.datalen])
    tst.y <- data.frame(y=dat$y[test.indx:current.datalen])

    f = as.formula("y ~ x")
    trn.noZeros <- data.frame(subset(trn.data, trn.data[,1]>0))
    out <<- neuralnet(f,
                      trn.noZeros,
                      hidden=hidden.nodes,
                      rep=2,
                      stepmax = 2000,
                      threshold=0.2,
                      learningrate=1,
                      algorithm="rprop+", #'rprop-', 'sag', 'slr'
                      startweights=NULL,
                      lifesign="none",
                      err.fct="sse",
                      act.fct="logistic",
                      exclude = NULL,
                      constant.weights = NULL,
                      linear.output=TRUE #If true, act.fct is not applied to the o/p of neuron. So it will be only integartion function
    )


    data.test <<- c(data.test, tst.y$y)

    pred <- compute(out, tst.x)$net.result
    data.out <<- c(data.out, pred)

    indx.end <<- indx.end + horizon
    count <- count + 1
  }

  test.data <<-  cbind(test.data, data.test)
  output <<- cbind(output, data.out)

  test.data.denorm <<- cbind(test.data.denorm, as.vector(denormalizeData(data.test,pow.normParms)))
  output.denorm <<- cbind(output.denorm, as.vector(denormalizeData(data.out, pow.normParms)))

}


predict.for.combination <- function(){
  gen.aggrdata()
  if(aggr.mat.size != 0){
    for(aggr.indx in seq(1,aggr.mat.size)){
      predict.pow(aggr.indx)
    }
  }
  else{
    predict.pow(0)
  }
}

measure.error <- function(pred,test){
  err.rmse <- error(forecast=pred, true=test,method="rmse")
  err.mape <- error(forecast=pred, true=test,method="mape")
  err.mae <- error(forecast=pred, true=test,method="mae")
  err.mse <- error(forecast=pred, true=test,method="mse")
  err.sd <- sd(pred-test)# need to normalize it to installed power
  err.cor <- cor(pred,test)

  return(c(err.rmse, err.mape, err.mae, err.mse, err.sd, err.cor))
}

prediction.error <- function(){
  parm.count <- 7
  setwd(filepath)
  col.names <- c("AggrNo.Seq","rmse", "mape", "mae", "mse", "sd", "cor")
  input.data.type <- c("norm","denorm")

  for(type in input.data.type){
    if(aggr.cluster.size == 1){
      err.data <- matrix(0,nrow=1,ncol=parm.count, byrow=TRUE)
      file.name <- ''
      colnames(err.data) <- col.names
      site.name <- paste('S',paste(seq(1,aggr.cluster.size), collapse=""), sep="")
      test <- rowSums(test.data.denorm[,1:aggr.mat.size])
      pred <- rowSums(output.denorm[,1:aggr.mat.size])

      err.data[1,] <- switch(type,
                             "norm"={
                               test <- normalizeData(test, type="0_1")
                               pred <- normalizeData(pred, type="0_1")
                               err <- measure.error(pred,test)
                               file.name <- paste(file.name.generic, aggr.cluster.size,'.csv',sep="")
                               c(site.name,err)
                             },
                             "denorm"={
                               err <- measure.error(pred,test)
                               file.name <- paste(file.name.denorm.generic, aggr.cluster.size,'.csv',sep="")
                               c(site.name, err)
                             }
      )
      write.csv(err.data, file=file.name)

    }else if(aggr.cluster.size>1 && aggr.cluster.size < sites.count){
      err.data <- matrix(0,nrow=indxcombicnt, ncol=parm.count, byrow=TRUE)
      colnames(err.data) <- col.names
      site.names <- colnames(pow.aggrdata)

      for(site in seq(1:(aggr.mat.size))){
        site.name <- site.names[site]

        err.data[site,] <- switch(type,
                                  "norm"={
                                    test <- test.data[,site]
                                    pred <- output[,site]
                                    err <- measure.error(pred,test)
                                    file.name <- paste(file.name.generic, aggr.cluster.size,'.csv',sep="")
                                    c(site.name, err)
                                  },
                                  "denorm"={
                                    test.denorm <- test.data.denorm[,site]
                                    pred.denorm <- output.denorm[,site]
                                    err <- measure.error(pred.denorm, test.denorm)
                                    file.name <- paste(file.name.denorm.generic, aggr.cluster.size,'.csv',sep="")
                                    c(site.name, err)
                                  }
        )
        write.csv(err.data, file=file.name)
      }

    }else{
      err.data <- matrix(0,nrow=1,ncol=parm.count, byrow=TRUE)
      colnames(err.data) <- col.names
      site.name <- paste('S',paste(seq(1,aggr.cluster.size), collapse=""), sep="")

      err.data[1,] <- switch(type,
                             "norm"={
                               test <- test.data[,1]
                               pred <- output[,1]
                               err <- measure.error(pred, test)
                               file.name <- paste(file.name.generic, aggr.cluster.size,'.csv',sep="")
                               c(site.name, err)
                             },
                             "denorm"={
                               test.denorm <- test.data.denorm[,1]
                               pred.denorm <- output.denorm[,1]
                               err <- measure.error(pred.denorm, test.denorm)
                               file.name <- paste(file.name.denorm.generic, aggr.cluster.size,'.csv',sep="")
                               c(site.name, err)
                             })
      write.csv(err.data, file=file.name)
    }
  }

}

predict.all.combination <- function(){
  for(aggr in aggr.type.vec){
    aggr.type <<- aggr
    filepath <<- paste(filepath.generic, aggr, '/', sep="")

    loaddata()
    for(combi in seq(10,10)){#sites.count
      aggr.cluster.size <<- combi
      if(combi != sites.count){
        indxcombicnt <<-length(combn(sites.count,combi)[1,])
        aggr.mat.size <<- indxcombicnt
        pow.aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
        solar.aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
        indxcombimat <<- as.matrix(combn(sites.count, aggr.cluster.size))
      }else{
        aggr.cluster.size <<- combi
        indxcombicnt <<- 0
        aggr.mat.size <<- 0
        pow.aggrdata10 <<- ff(NA, dim=data.len, vmode="double")
        solar.aggrdata10 <<- ff(NA, dim=data.len, vmode="double")
      }

      train.data <<- c()
      test.data <<- c()
      test.data.denorm <<- c()
      output <<- c()
      output.denorm <<- c()

      predict.for.combination()
      prediction.error()
    }
  }
}

predict.all.combination()
