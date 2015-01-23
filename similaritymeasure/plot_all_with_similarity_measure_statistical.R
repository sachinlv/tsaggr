require(RMySQL)
require(ff)
require(RSNNS)
require(TSdist)
algo.vec <- c('neuralnet','brnn','gbm')
simi.meas.vec <- c('euclideanDistance',
                   'minkowskiDistance',
                   'manhattanDistance',
                   'fourierDistance',
                   'correlationDistance')

algo <- 'gbm'
sim.meas <- 'pca'
ip.type <- 'statistical'
file.path.single <- paste('/home/freak/Dropbox/results/',algo,'_shortterm_simple/',sep="")
file.path.all <- paste('/home/freak/Dropbox/results/',algo,'_shortterm_aggr/all/',sep="")
plot.file <- paste('/home/freak/Programming/Thesis/results/plots/',
                   algo,'_',
                   ip.type,'_',
                   sim.meas,'.pdf',sep="")
result.file <- paste('/home/freak/Programming/Thesis/results/plots/combination/',
                      algo,'_',
                      ip.type,'_',
                      sim.meas,'.csv',sep="")

tables <- c('onshore_SITE_00002',
            'onshore_SITE_00003',
            'onshore_SITE_00004',
            'onshore_SITE_00005',
            'onshore_SITE_00006',
            'onshore_SITE_00007',
            'onshore_SITE_00008',
            'onshore_SITE_00012',
            'onshore_SITE_00013',
            'onshore_SITE_00014')

sites.count <- 10
data.len <- 52560
file.name.single <- paste(algo,'_shortterm_simple.csv',sep="")
file.name.all <- paste(algo,'_shortterm_aggr_combi',sep="")
powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
threshold.rmse <<- 0.3
threshold.dist <<- 0.3

combination.result <<- c()

file.list <- list.files(file.path.all, pattern="*.csv")
drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")

load.data <- function(){
  for(indx in seq(1,sites.count)){
    tab <- tables[indx]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]),type="0_1")
  }
}

load.err.data <- function(file.no){
  rmse <<- c()

  if(file.no==1){
    file <- file.name.single
    err.tbl <- read.csv(file,sep = ',')
    err <- err.tbl$rmse
    rmse <<- data.frame(err)
  }
  else{
    file <- paste(file.name.all,file.no,'.csv', sep="")
    err.tbl <- read.csv(file,sep = ',')
    seq <- err.tbl$AggrNo.Seq
    err <- err.tbl$rmse
    seq <- gsub("10","0",seq)#Replace 10 with 0 in vector
    rmse <<- data.frame(seq,err)
  }
}


gen.dist.mat <- function(aggr.seq){
  len <- length(aggr.seq)
  mat <- matrix(0,nrow=len, ncol=len)
  print(aggr.seq)


  for(i in seq(1,len)){
    for(j in seq(1, len)){
      i1 <- aggr.seq[i]
      j1 <- aggr.seq[j]
      data.mat <- data.frame(c1=powdata[,i1], c2=powdata[,j1])
      #pca-start
      mat[i,j] <- switch(sim.meas,
                  "pca"={
                    f <- formula("~ c1 + c2")
                    pc <- princomp(f,data=data.mat, cor=FALSE)
                    score <- data.frame(pc$scores)
                    euclideanDistance(score$Comp.1,score$Comp.2)
                  },
                  "euclidean"={
                    euclideanDistance(data.mat$c1,data.mat$c2)
                  },
                  "correlation"={
                    correlationDistance(data.mat$c1,data.mat$c2)
                  },
                  "manhattan"={
                    manhattanDistance(data.mat$c1,data.mat$c2)
                  },
                  "minkowski"={
                    minkowskiDistance(data.mat$c1,data.mat$c2,2)
                  },
                  "fourier"={
                    fourierDistance(data.mat$c1,data.mat$c2,n=(floor(length(data.mat$c1)/2)+1))
                  },
                  "mean"={
                    abs(mean(data.mat$c1-data.mat$c2))
                  },
                  "mahalanobis"={
                    d <- matrix(NA,nrow=52560,ncol=2)
                    #d[,1] <- powdata[,i1]
                    #d[,2] <- powdata[,j1]
                    #c <- cov(d[,1:2])#,method=c("pearson"))
                    #m <- c(mean(d[,1]), mean(d[,2]))
                    #mat[i,j] <- mahalanobis(d,m,c)
                  }
             )
    }
  }
  #print(mat)
  return(mat)
}


get.err.dist.data <- function(){
  avg.dist <- c()
  for(i in seq(1,length(rmse$seq))){
    seq <- unlist(strsplit(as.character(rmse$seq[i]),''))
    seq <- as.integer(seq[2:length(seq)]) #remove S from sequence
    seq <- replace(seq,seq==0,10)#replace 0 back to 10
    diag.len <- row.len <- col.len <- length(seq)
    dist.mat <- gen.dist.mat(seq)
    #dist.mat[upper.tri(dist.mat,diag=TRUE)] <- 0
    #lower.tri.cnt <- ((row.len*col.len)-diag.len)/2
    #dist <- sum(dist.mat)/lower.tri.cnt
    dist <- sum(dist.mat)/((row.len*col.len)-diag.len)
    avg.dist <- c(avg.dist,dist)
  }

  if(length(avg.dist) == 1){
    while(avg.dist>1){avg.dist <- avg.dist/10}
  }
  else{
    avg.dist <- normalizeData(as.vector(avg.dist),type="0_1")
  }
  #combination result set
  for(i in seq(1,length(avg.dist))){
    if(avg.dist[i] <= threshold.dist && rmse$err[i] <= threshold.rmse){
      combination.result <<- c(combination.result, rmse$seq[i])
    }
  }

  data.to.plot <- data.frame(err.rmse=rmse$err, dist=avg.dist)
  return(data.to.plot)
}



plot.all <- function(){
  load.data()
  #x11()
  par(mfrow=c(5,2))
  par(mar=c(0.5, 4.5, 0.5, 0.5))
  plot.window(xlim=c(0,300),ylim=c(0,1),asp=1)
  #plot singles
  setwd(file.path.single)
  load.err.data(1)
  data.to.plot <- rmse$err
  plot(data.to.plot,
       type="p",
       ylim=c(0,1),
       col="red",
       xlab=paste("Individuals"),
       ylab="Rmse vals",
       main="Individuals")
  legend("topright",
         legend=c("rmse"),
         col=c("red"),
         text.col=c("red"))

  #plot2-10
  setwd(file.path.all)
  for(i in seq(2,10)){
    load.err.data(i)
    data.to.plot <- get.err.dist.data()
    #typ <- if(i == 10) 'p' else 'l'
    y <- data.to.plot$err.rmse
    x <- data.to.plot$dist
    d <- data.frame(y,x)
    f <- formula('y~x')
    plot(f,d,
         main=cor(data.to.plot$err.rmse,data.to.plot$dist),
         #type="p",
         #col=c("red","blue"),
         xlab=paste("rmse combination ",i),
         ylab=paste("distance combination ",i))
         #main="Combination")neuralnet_statistical_euclidean.pdf
         #plot(data.to.plot$dist)
         #if(i<10){
         #  lines(data.to.plot$dist,col="blue")
         #}
         #else{
         #  print(data.to.plot$dist)
         #points(data.to.plot$dist,,col="blue")
         #}
         #legend("topright", legend=c("rmse","correlation"),
         #           col=c("red","blue"), text.col=c("red","blue"))
  }
  dev.copy2pdf(file =plot.file)
  #dev.off()
  combination.result <- as.character(combination.result)
  combination.result <- gsub("0","10",combination.result)
  df <- c(paste("RMSE threshold", as.character(threshold.rmse)))
  df <- c(df,paste("Similarity Threshold",as.character(threshold.dist)))
  #combination.result <- cbind(numeric(0),combination.result)
  df <- c(df, combination.result)
  write.table(df, result.file)
}

plot.all()


#PCA
mat <- as.matrix(powdata[,1:10])
colnames(mat) <- c(paste('c',seq(1,10),sep=""))
f <- formula("~ c1 + c2")


a <- prcomp(f,data=data.frame(mat))
p <- a$x
#or
b <- princomp( f,data=data.frame(mat), cor=FALSE)
p <- b$scores

o <- p[,1:2]
o <- as.matrix(o)
head(o)
head(mat[,1:2])
dist.mat <- matrix(NA,nrow=2,ncol=2)
for(i in seq(1,2)){
  for(j in seq(1,2)){
    dist.mat[i,j] <- euclideanDistance(o[,i] , o[,j])
  }
}
dist.mat
euclideanDistance(o[,1],o[,5])
?princomp
?euclideanDistance
?dist

su <- (o[,1]-o[,2])^2
sum(su)
s <- 0
for(i in su){
  s <- s + i
  print(s)
}
