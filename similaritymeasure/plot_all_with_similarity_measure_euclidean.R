require(RMySQL)
require(ff)
require(RSNNS)
require(TSdist)
algo <- 'gbm'
file.path.single <- paste('/home/freak/Dropbox/results/',algo,'_shortterm_simple/',sep="")
file.path.all <- paste('/home/freak/Dropbox/results/',algo,'_shortterm_aggr/all/',sep="")

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
      mat[i,j] <- euclideanDistance(powdata[,i1],powdata[,j1])#,2)#, n=(floor(length(powdata[,i1])/2)+1))
    }
  }

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
    dist <- sum(dist.mat)/((row.len+col.len)-diag.len)

    avg.dist <- c(avg.dist,dist)
  }

  if(length(avg.dist) == 1){
    while(avg.dist>1){avg.dist <- avg.dist/10}
  }
  else{
    avg.dist <- normalizeData(as.vector(avg.dist),type="0_1")
  }

  data.to.plot <- data.frame(err.rmse=rmse$err, eucd.dist=avg.dist)
  return(data.to.plot)
}



plot.all <- function(){
  load.data()
  par(mfrow=c(5,2))
  par(mar=c(0.5, 4.5, 0.5, 0.5))
  plot.window(xlim=c(0,300),ylim=c(0,1),asp=1)
  #plot singles
  setwd(file.path.single)
  load.err.data(1)
  data.to.plot <- rmse$err
  plot(data.to.plot,type="l",
       ylim=c(0,1),col="red",
       xlab=paste("Individuals"),ylab="Rmse vals",
       main="Individuals")
  legend("topright", legend=c("rmse"),
         col=c("red"), text.col=c("red"))

  #plot2-10
  setwd(file.path.all)
  for(i in seq(2,10)){
    load.err.data(i)
    data.to.plot <- get.err.dist.data()
    typ <- if(i == 10) 'p' else 'l'
    plot(data.to.plot$err.rmse,type=typ,
         ylim=c(0,1),col="red",
         xlab=paste("combination",i),ylab="Rmse val~distance",
         main="Combination")
    #if(i<10){
    #  lines(data.to.plot$eucd.dist,col="blue")
    #}
    #else{
    #  print(data.to.plot$eucd.dist)
      points(data.to.plot$eucd.dist,,col="blue")
    #}
    legend("topright", legend=c("rmse","euclideanDistance"),
                      col=c("red","blue"), text.col=c("red","blue"))
  }

}

plot.all()
