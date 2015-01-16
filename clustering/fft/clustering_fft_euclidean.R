require(RMySQL)
require(ff)
require(RSNNS)
require(TSdist)

sites.count <- 10
data.len <- 52560
powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.fft <<- matrix(, data.len, sites.count) #ff(NA, dim=c(data.len, sites.count), vmode="complex")
powdata.normalized.fft <<-  matrix(, data.len, sites.count)#ff(NA, dim=c(data.len, sites.count), vmode="complex")
dist.mat.real <<- matrix(, nrow=sites.count, ncol=sites.count)
dist.mat.img <<- matrix(, nrow=sites.count, ncol=sites.count)
dist.mat.real.norm <<- matrix(, nrow=sites.count, ncol=sites.count)
dist.mat.img.norm <<- matrix(, nrow=sites.count, ncol=sites.count)


drv <- dbDriver("MySQL")
con <- dbConnect(drv, host="localhost", dbname="eastwind", user="sachin", pass="password")
tablelist_statement = paste("SELECT TABLE_NAME FROM information_schema.TABLES ",
                            "WHERE TABLE_SCHEMA = 'eastwind' AND",
                            "TABLE_NAME LIKE 'onshore_SITE_%' "," LIMIT ",sites.count, ";")
tables <- dbGetQuery(con, statement=tablelist_statement)
tables <- data.frame(tables)


loaddata <- function(){
  colindx <- 1
  for(indx in seq(1,sites.count)){
    tab <- tables[indx,]
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]), type="0_1")
  }
}

generate.fft <- function(){
  for(indx in seq(1,sites.count)){
    powdata.fft[,indx] <<- fft(powdata[,indx])
    powdata.normalized.fft[,indx] <<- fft(powdata.normalized[,indx])
  }
}

generate.distance.matrix <- function(){
  for(i in seq(1,sites.count)){
    for(j in seq(1, sites.count)){
      dist.mat.real.norm[i,j] <<- euclideanDistance(Re(powdata.fft[,i]),Re(powdata.fft[,j]))
      dist.mat.img.norm[i,j] <<- euclideanDistance(Im(powdata.fft[,i]),Im(powdata.fft[,j]))
    }
  }
}

generate.distance.matrix.norm <- function(){
  for(i in seq(1,sites.count)){
    for(j in seq(1, sites.count)){
      dist.mat.real[i,j] <<- euclideanDistance(Re(powdata.normalized.fft[,i]),Re(powdata.normalized.fft[,j]))
      dist.mat.img[i,j] <<- euclideanDistance(Im(powdata.normalized.fft[,i]),Im(powdata.normalized.fft[,j]))
    }
  }
}


generate.hierarchial.cluster <- function(){
  loaddata()
  generate.fft()
  generate.distance.matrix()
  generate.distance.matrix.norm()

  hc.real <<- hclust(as.dist(dist.mat.real))
  hc.img <<- hclust(as.dist(dist.mat.img))

  hc.real.norm <<- hclust(as.dist(dist.mat.real.norm))
  hc.img.norm <<- hclust(as.dist(dist.mat.img.norm))
}

generate.hierarchial.cluster()

plot(hc.real.norm)
summary(hc.real.norm)

groups <- cutree(hc.real.norm,k=5)
rect.hclust(hc.real.norm, k=5, border="red")
