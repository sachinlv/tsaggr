require(gbm)
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
history.length <- 50
data.len.day <<- 144
data.len <- history.length * data.len.day
window.size <- 10
train.data.percent <- 0.7
start.date <- '20061112'
end.date <- '20070101'
#mat.size <<- 365
#slide.count <- mat.size-window.size+1
filepath.generic <<- '/home/freak/Programming/Thesis/results/history50/random_sites5/gbm_shortterm_'
file.name.generic <<- 'gbm_shortterm_aggr_combi'
file.name.denorm.generic <<- 'gbm_shortterm_aggr_combi_denorm'
file.name.aggr.generic <<- 'gbm_shortterm_aggr_combi_aggr'
file.name.aggr.denorm.generic <<- 'gbm_shortterm_aggr_combi_aggr_denorm'

table.ip.type <- "random"#c("random","specific")
preprocess <<- 'raw' #c('raw','linear', 'expsmooth', 'wma')
aggr.type.vec <<- c('aggr', 'mean','wmean')

powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")

drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")

if(table.ip.type == "random"){
  t <- c("onshore_SITE_00002",
         "onshore_SITE_00003",
         "onshore_SITE_00004",
         "onshore_SITE_00005",
         "onshore_SITE_00006",
         "onshore_SITE_00007",
         "onshore_SITE_00008",
         "onshore_SITE_00012",
         "onshore_SITE_00013",
         "onshore_SITE_00014")
  tables <<- data.frame(cbind(numeric(0),t))
}else{
  t <- c("onshore_SITE_00538",
         "onshore_SITE_00366",
         "onshore_SITE_00623",
         "onshore_SITE_00418",
         "onshore_SITE_00627",
         "onshore_SITE_00532",
         "onshore_SITE_00499",
         "onshore_SITE_00571",
         "onshore_SITE_03247",
         "onshore_SITE_00622")
  tables <<- data.frame(cbind(numeric(0),t))
}

loaddata <- function(){
  colindx <- 1
  for(indx in seq(1,sites.count)){
    tab <- tables[indx,]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= ",start.date," && mesdt < ",end.date,") LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- switch(preprocess,
                              "raw"={
                                as.double(data06[,1])
                              },
                              "linear"={
                                d <- as.double(data06[,1])
                                f <- filter(d,rep(1/mean(d),mean(d)))
                                f[is.na(f)] <- 0
                              },
                              "expsmooth"={
                                d <- as.double(data06[,1])
                                f <- tbats(d,
                                           use.trend=FALSE,
                                           use.damped.trend=FALSE,
                                           use.box.cox=FALSE,
                                           use.arma.errors=FALSE,
                                           use.parallel=TRUE,
                                           num.cores=3)
                                f$fitted.values
                              },
                              "wma"={
                                c()
                              }
    )
  }
}


aggr.timeseries <- function(vec){
  len <- length(vec)/4
  aggr.indx <-rep(c(seq(1,len)),each=4)
  x <- as.zoo(vec)
  aggr <- as.numeric(aggregate(x,by=aggr.indx, FUN="mean"))
  return(aggr)
}

gen.aggrdata <- function(){
  col.names <- c()
  if(aggr.mat.size != 0){
    for(i in seq(1 ,aggr.mat.size)){
      indx.seq <<- indxcombimat[,i]
      mat <- as.matrix(powdata[,indx.seq])
      col.names <- c(col.names, paste('S',paste(indx.seq, collapse=""), sep=""))
      if(length(mat[1,]) > 1){
        aggrdata[,i] <<- rowSums(mat)
        aggrdata.mean[,i] <<- rowMeans(mat)
        aggrdata.wmean[,i] <<- rowWeightedMeans(mat)
      }else{
        aggrdata[,i] <<- mat[,1]
        aggrdata.mean[,i] <<- mat[,1]
        aggrdata.wmean[,i] <<- mat[,1]
      }
    }
    colnames(aggrdata) <<- col.names
    colnames(aggrdata.mean) <<- col.names
    colnames(aggrdata.wmean) <<- col.names
  }
  else{
    indx.seq <- seq(1,sites.count)
    aggrdata10 <<- rowSums(as.matrix(powdata[,indx.seq]))
    aggrdata10.mean <<- rowMeans(as.matrix(powdata[,indx.seq]))
    aggrdata10.wmean <<- rowWeightedMeans(as.matrix(powdata[,indx.seq]))
  }
}


predict.pow <- function(aggrno) {
  data.normalized <- c()
  if(aggr.mat.size!=0){
    data.normalized <- switch(aggr.type,
                              'aggr'={
                                normalizeData(as.vector(aggrdata[,aggrno]),type="0_1")
                              },
                              'mean'={
                                normalizeData(as.vector(aggrdata.mean[,aggrno]),type="0_1")
                              },
                              'wmean'={
                                normalizeData(as.vector(aggrdata.wmean[,aggrno]),type="0_1")
                              }
    )
  }
  else{
    data.normalized <- switch(aggr.type,
                              'aggr'={
                                normalizeData(as.vector(aggrdata10),type="0_1")
                              },
                              'mean'={
                                normalizeData(as.vector(aggrdata10.mean),type="0_1")
                              },
                              'wmean'={
                                normalizeData(as.vector(aggrdata10.wmean),type="0_1")
                              }
    )
  }

  normParms <<- getNormParameters(data.normalized)
  data.set <- as.vector(data.normalized[,1])
  indx.start <<- 1
  indx.end <<- indx.start + (window.size * data.len.day) - 1
  data.train <<- c()
  data.test <<- c()
  data.out <<- c()
  count <- 0

  while(indx.end <= data.len){
    print(paste("Cluster size: ",aggr.cluster.size," AggrNo.: ", aggrno, " Slide No.: ", count+1))
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
    #After test, showed that not much difference of n.trees=10000(0.231001) and n.trees=1000(0.232469)
    out <<- gbm(f,
                data=data.mat.train,
                distribution ="gaussian",
                n.trees=1000,
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

  }

  train.data <<- cbind(train.data, data.train)
  test.data <<-  cbind(test.data, data.test)
  output <<- cbind(output, data.out)

  test.data.denorm <<- cbind(test.data.denorm, as.vector(denormalizeData(data.test,normParms)))
  output.denorm <<- cbind(output.denorm, as.vector(denormalizeData(data.out, normParms)))

}


predict.for.combination <- function(){
  gen.aggrdata()
  if(aggr.mat.size != 0){
    for(aggr.indx in seq(1,aggr.mat.size)){
      predict.pow(aggr.indx)
      #break
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
  input.data.type <- c("norm","denorm", "aggrnorm", "aggrdenorm")

  for(type in input.data.type){
    if(aggr.cluster.size == 1){
      err.data <- matrix(0,nrow=1,ncol=parm.count, byrow=TRUE)
      file.name <- ''
      colnames(err.data) <- col.names
      site.name <- paste('S',paste(seq(1,aggr.cluster.size), collapse=""), sep="")
      test <- switch(aggr.type,
                     'aggr'={
                       rowSums(test.data.denorm[,1:aggr.mat.size])
                     },
                     'mean'={
                       rowMeans(test.data.denorm[,1:aggr.mat.size])
                     },
                     'wmean'={
                       rowWeightedMeans(test.data.denorm[,1:aggr.mat.size])
                     }
      )
      pred <- switch(aggr.type,
                     'aggr'={
                       rowSums(output.denorm[,1:aggr.mat.size])
                     },
                     'mean'={
                       rowMeans(output.denorm[,1:aggr.mat.size])
                     },
                     'wmean'={
                       rowWeightedMeans(output.denorm[,1:aggr.mat.size])
                     }
      )

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
                             },
                             "aggrnorm"={
                               test <- normalizeData(test, type="0_1")
                               pred <- normalizeData(pred, type="0_1")
                               test <- aggr.timeseries(test)
                               pred <- aggr.timeseries(pred)
                               err <- measure.error(pred, test)
                               file.name <- paste(file.name.aggr.generic,aggr.cluster.size,'.csv',sep="")
                               c(site.name, err)
                             },
                             "aggrdenorm"={
                               test <- aggr.timeseries(test)
                               pred <- aggr.timeseries(pred)
                               err <- measure.error(pred,test)
                               file.name <- paste(file.name.aggr.denorm.generic,aggr.cluster.size,'.csv',sep="")
                               c(site.name, err)
                             }
      )
      write.csv(err.data, file=file.name)

    }else if(aggr.cluster.size>1 && aggr.cluster.size < sites.count){
      err.data <- matrix(0,nrow=indxcombicnt, ncol=parm.count, byrow=TRUE)
      colnames(err.data) <- col.names
      site.names <- colnames(aggrdata)

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
                                  },
                                  "aggrnorm"={
                                    test <- aggr.timeseries(test.data[,site])
                                    pred <- aggr.timeseries(output[,site])
                                    err <- measure.error(pred, test)
                                    file.name <- paste(file.name.aggr.generic,aggr.cluster.size,'.csv',sep="")
                                    c(site.name, err)
                                  },
                                  "aggrdenorm"={
                                    test <- aggr.timeseries(test.data.denorm[,site])
                                    pred <- aggr.timeseries(output.denorm[,site])
                                    err <- measure.error(pred,test)
                                    file.name <- paste(file.name.aggr.denorm.generic,aggr.cluster.size,'.csv',sep="")
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
                             },
                             "aggrnorm"={
                               test <- aggr.timeseries(test.data[,1])
                               pred <- aggr.timeseries(output[,1])
                               err <- measure.error(pred, test)
                               file.name <- paste(file.name.aggr.generic,aggr.cluster.size,'.csv',sep="")
                               c(site.name, err)
                             },
                             "aggrdenorm"={
                               test <- aggr.timeseries(test.data.denorm[,1])
                               pred <- aggr.timeseries(output.denorm[,1])
                               err <- measure.error(pred,test)
                               file.name <- paste(file.name.aggr.denorm.generic,aggr.cluster.size,'.csv',sep="")
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
      if(combi != sites.count){
        aggr.cluster.size <<- combi
        indxcombicnt <<-length(combn(sites.count,combi)[1,])
        aggr.mat.size <<- indxcombicnt

        aggrdata <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
        aggrdata.mean <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")
        aggrdata.wmean <<- ff(NA, dim=c(data.len, aggr.mat.size), vmode="double")

        indxcombimat <<- as.matrix(combn(sites.count, aggr.cluster.size))
      }else{
        aggr.cluster.size <<- combi
        indxcombicnt <<- 0
        aggr.mat.size <<- 0
        aggrdata10 <<- ff(NA, dim=data.len, vmode="double")
        aggrdata10.mean <<- ff(NA, dim=data.len, vmode="double")
        aggrdata10.wmean <<- ff(NA, dim=data.len, vmode="double")
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
