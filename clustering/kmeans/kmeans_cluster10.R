require(RMySQL)
require(RSNNS)
require(NbClust)
require(ff)

data.len <- 52560
sites.count <- 10
table.ip.type <- "random"#"specific"

powdata <<- matrix(NA, nrow=data.len, ncol=sites.count)
powdata.normalized <<- matrix(NA, nrow=data.len, ncol=sites.count)
#powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
#powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")

drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")

if(table.ip.type == "random"){
  tablelist_statement = paste("SELECT TABLE_NAME FROM information_schema.TABLES ",
                              "WHERE TABLE_SCHEMA = 'eastwind' AND",
                              "TABLE_NAME LIKE 'onshore_SITE_%' "," LIMIT ",sites.count, ";")
  tables <<- dbGetQuery(con, statement=tablelist_statement)
  tables <<- data.frame(tables)
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
  for(indx in seq(1,sites.count)){
    tab <- tables[indx,1]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]),type="0_1")
  }
}

generate.cluster <- function(){
  loaddata()
  d <- t(powdata)
  d2 <- powdata.normalized[,1:3]
  #d2 <- t(d2)

  c <- NbClust(d2[1:2000,1:3],
                distance="euclidean",
                min.nc=2, max.nc=4,
                method="kmeans")

  #kmeans(powdata,2,iter.max=10)
  c2 <- kmeans(d[1:10,],2,iter.max=10)
  c3 <- kmeans(d[1:10,],3,iter.max=10)
  c4 <- kmeans(d[1:10,],4,iter.max=10)
  c5 <- kmeans(d[1:10,],5,iter.max=10)
  c6 <- kmeans(d[1:10,],6,iter.max=10)
  c7 <- kmeans(d[1:10,],7,iter.max=10)
  c8 <- kmeans(d[1:10,],8,iter.max=20)
  c9 <- kmeans(d[1:10,],9,iter.max=20)
}

dim(fitted(c2))
clust <- c2$cluster
clust
c2$withinss
c3$size
c4$size
c5$size
c6$size
c7$size
c8$size
c9$size




set.seed(1)
x<-rbind(matrix(rnorm(150,sd=0.3),ncol=5),
         matrix(rnorm(150,mean=3,sd=0.2),ncol=5),
         matrix(rnorm(150,mean=1,sd=0.1),ncol=5),
         matrix(rnorm(150,mean=6,sd=0.3),ncol=5),
         matrix(rnorm(150,mean=9,sd=0.3),ncol=5))
head(x)

res<-NbClust(x, distance = "euclidean", min.nc=2, max.nc=10,
             method = "ward.D", index = "all")

res$All.index
res$Best.nc
res$All.CriticalValues
res$Best.partition
