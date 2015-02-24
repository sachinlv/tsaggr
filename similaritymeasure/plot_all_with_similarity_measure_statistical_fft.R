require(RMySQL)
require(ff)
require(RSNNS)
require(TSdist)

table.ip.type <- "random"#"specific"
#results.file.path <-'/home/freak/Dropbox/results/specific_sites/texas/'
#plot.file.path <- '/home/freak/Programming/Thesis/results/plots/specific_sites/texas/'
#combination.file.path <- '/home/freak/Programming/Thesis/results/plots/specific_sites/texas/combination/'
results.file.path <-'/home/freak/Dropbox/results/random_sites/'
plot.file.path <- '/home/freak/Programming/Thesis/results/plots/random_sites/'
combination.file.path <- '/home/freak/Programming/Thesis/results/plots/random_sites/combination/'


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


sites.count <- 10
data.len <-7200 #52560
powdata <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
powdata.normalized <<- ff(NA, dim=c(data.len, sites.count), vmode="double")
threshold.rmse <<- 0.3
threshold.dist <<- 0.3
drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")


load.data <- function(){
  for(indx in seq(1,sites.count)){
    tab <- as.character(tables[indx,1])
    print(paste("Loading from table :: ", tab))
    query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20061112 && mesdt < 20070101) LIMIT ", data.len, ";")
    data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
    powdata[,indx] <<- as.double(data06[,1])
    powdata.normalized[,indx] <<- normalizeData(as.vector(data06[,1]),type="0_1")
  }
}

load.err.data <- function(file.no){
  rmse <<- c()

  file <- paste(file.name.all,file.no,'.csv', sep="")
  err.tbl <- read.csv(file,sep = ',')
  seq <- err.tbl$AggrNo.Seq
  err <- err.tbl$rmse
  seq <- gsub("10","0",seq)#Replace 10 with 0 in vector
  rmse <<- data.frame(seq,err)
}

load.data()

head(powdata[,1:10])
f <- fft(as.matrix(powdata[,1:10]))
f.res <- f[1:3600,1:10]

f.Re <- Re(f.res)
f.Im <- Im(f.res)
f.Re

k <- kmeans(f.Re, 3)
k
summary(k)
k$centers
quantile(f.Re)
library('TSclust')
x <- powdata[1:7200,1]
y <- powdata[1:7200,2]
w <- 20
alpha <- 4
n = 100
x <- (x - mean(x)) /sd(x)
y <- (y - mean(y)) /sd(y)
paax <- PAA(x, w) #generate PAA reductions
paay <- PAA(y, w)
plot(x, type="l", main="PAA reduction of series x") #plot an example of PAA reduction
p <- rep(paax,each=length(x)/length(paax)) #just for plotting the PAA
lines(p, col="red")
#repeat the example with y
plot(y, type="l", main="PAA reduction of series y")
py <- rep(paay,each=length(y)/length(paay))
lines(py, col="blue")
#convert to SAX representation
SAXx <- convert.to.SAX.symbol( paax, alpha)
SAXy <- convert.to.SAX.symbol( paay, alpha)
#CALC THE SAX DISTANCE
MINDIST.SAX(SAXx, SAXy, alpha, n)
#this whole process can be computed using diss.MINDIST.SAX
diss.MINDIST.SAX(x, y, w, alpha, plot=TRUE)
z <- rnorm(n)^2
diss(rbind(x,y,z), "MINDIST.SAX", w, alpha)
SAX.plot( as.ts(cbind(x,y,z)), w=w, alpha=alpha)



