require(brnn)
require(RMySQL)
require(ff)
require(googleVis)
require(Metrics)
require(ppls)
require(combinat)

sites.count <- 10
indxcombicnt <<- 10 #no. of combination
indxseq <<- c(seq(1,sites.count))
aggr.mat.size <- sites.count/2
row.size <- 52560
row.size.day <<- 144
mat.size <<- 365
window.size <- 30
train.data.percent <- 0.7
indxseq <- c(seq(1,sites.count))

powdata <<- ff(NA, dim=c(row.size, sites.count), vmode="double")
aggrdata <<- ff(NA, dim=c(row.size, aggr.mat.size), vmode="double")
train.data <<- ff(NA, dim=c(row.size.day*train.data.percent, mat.size-window.size+1), vmode="double")
test.data <<- ff(NA, dim=c(row.size.day*(1-train.data.percent), mat.size-window.size+1), vmode="double")
output <<- ff(NA, dim=c(row.size.day*(1-train.data.percent), mat.size-window.size+1), vmode="double")
indxcombimat <<- ff(NA, dim=c(aggr.mat.size,indxcombicnt),vmode="integer")

drv = dbDriver("MySQL")
con = dbConnect(drv,host="localhost",dbname="eastwind",user="sachin",pass="password")
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
          query <- paste(" select pow from ", tab, " WHERE (mesdt >= 20060101 && mesdt < 20070101) LIMIT ", row.size, ";")
          data06 <- data.frame(dbGetQuery(con,statement=query), check.names=FALSE)
	      powdata[,indx] <<- as.double(data06[,1])
	    }
}

generate.seq.matrix <- function(){
	  mat <- as.matrix(combn(sites.count, aggr.mat.size)) #generating different combination
  len <- length(mat[1,])
    seq <- sample(1:len, indxcombicnt)
    indxcombimat <<- mat[,seq]
}

gen.aggrdata <- function(seq1, seq2){
	  seq1 <- as.vector(seq1)
  seq2 <- as.vector(seq2)
    for(i in seq(1 ,aggr.mat.size)){
	        indx1 <- seq1[i]
      indx2 <- seq2[i]
          aggrdata[,i] <<- c(powdata[,indx1]) + c(powdata[,indx2])
        }
}


predict <- function(indx) {
	  data.set <- matrix(aggrdata[,indx], nrow=row.size.day, ncol=mat.size, byrow=FALSE)
  print(head(data.set[,1:10]))
    for(indx in seq(1 ,mat.size-window.size+1)){
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
		          print(f)
			      out <<- brnn(f,
					                    data.mat.train,
							                     epochs=1000,
									                      cores=2,
											                       mu=0.005,
													                        mu_dec=0.1,
																                 mu_max=1e10,
																		                  change = 0.001,
																				                   neurons=window.size,
																						                    normalize=TRUE,
																								                     verbose=TRUE,
																										                      Monte_Carlo = FALSE)


			      train.data[, indx] <<- data.mat.train[,window.size]
			          test.data[, indx] <<- data.mat.test[,window.size]
			          pred <- predict.brnn(out ,data.mat.test[,1:window.size-1])
				      output[,indx] <<- as.double(pred)
				      plot(pred, type='l')
				          if(indx == 10){
						        break
				          }
				        }
}


predict_for_combination <- function(){
	  loaddata()
  generate.seq.matrix()

    for(i in seq(1:indxcombicnt)){
	        mat.indx1 <- as.vector(indxcombimat[,i])
      mat.indx2 <- as.vector(indxseq[is.na(pmatch(indxseq,mat.indx1))])
          gen.aggrdata(mat.indx1, mat.indx2)

          predict(i)
	      break
	    }
}

predict_for_combination()
generate.seq.matrix()


x1 = as.vector(train.data[,1:10])
x2 = as.vector(test.data[,1:10])
y = as.vector(output[,1:10])

plot(y, type='l')
plot(x1, type='l')
plot(x2, type='l')
length(x)
dataToPlot = data.frame(seq(1,430),x2,y)
Line <- gvisLineChart(dataToPlot)
plot(Line)
