require(RMySQL)
require(RSNNS)
require(TSdist)
require(ff)

sites.count <- 10
data.len <- 52560
table.ip.type <- "random"#"specific"
powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
#ff(NA, dim=c(data.len, sites.count), vmode="complex")
#ff(NA, dim=c(data.len, sites.count), vmode="complex")
dist.mat <<- matrix(0, nrow=sites.count, ncol=sites.count)
dist.mat.norm <<- matrix(0, nrow=sites.count, ncol=sites.count)


drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")

if(table.ip.type == "random"){
  t <- c('onshore_SITE_00002',
              'onshore_SITE_00003',
              'onshore_SITE_00004',
              'onshore_SITE_00005',
              'onshore_SITE_00006',
              'onshore_SITE_00007',
              'onshore_SITE_00008',
              'onshore_SITE_00012',
              'onshore_SITE_00013',
              'onshore_SITE_00014')
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
    tab <- as.character(tables[indx,1])
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]), type="0_1")
  }
}


generate.distance.matrix <- function(){
  for(i in seq(1,sites.count)){
    for(j in seq(1, sites.count)){
      dist.mat[i,j] <<- euclideanDistance(powdata[,i],powdata[,j])
      dist.mat.norm[i,j] <<- euclideanDistance(powdata.normalized[,i],
                                               powdata.normalized[,j])
    }
  }
}


generate.hierarchial.cluster <- function(){
  loaddata()
  generate.distance.matrix()
  hc <<- hclust(as.dist(dist.mat))
  hc.norm <<- hclust(as.dist(dist.mat.norm))
}


generate.hierarchial.cluster()

plot(hc)
hc$order
hc$merge
cutree(hc, k=4)
