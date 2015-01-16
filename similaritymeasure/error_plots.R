

file.path <-'/home/freak/Dropbox/results/neuralnet_shortterm_aggr/all/'

setwd(file.path)
file.list <- list.files(file.path, pattern="*.csv")
rmse <- c()
file.list

for(file in file.list){
  err.tbl <- read.csv(file,sep = ',')
  tmp <- err.tbl$rmse
  length(tmp) <- 252
  rmse <<- cbind(rmse,tmp)
}
col.names <- paste('combi', seq(2,10), sep="")
colnames(rmse) <- col.names

#matplot
matplot(rmse, type="l", col=1:9)

#lineplot
plot(rmse[,1],type="l",
     ylim=c(0,0.5),col="black",
     xlab="Count",ylab="Rmse val",
     main="Combination")
lines(rmse[,2],col="red",type="l")
lines(rmse[,3],col="orange",type="l")
lines(rmse[,4],col="yellow",type="l")
lines(rmse[,5],col="green",type="l")
lines(rmse[,6],col="blue",type="l")
lines(rmse[,7],col="purple",type="l")
lines(rmse[,8],col="pink",type="l")
lines(rmse[,9],col="darkgreen",type="l")

legend("topright",legend=c(paste("combi",c(10,2,3,4,5,6,7,8,9),sep="")),
       lty=1,lwd=2,pch=21,col=c("black","red","orange","yellow","green","blue","purple","pink","darkgreen"),
       ncol=2,bty="n",cex=0.8,
       text.col=c("black","red","orange","yellow","green","blue","purple","pink","darkgreen"),
       inset=0.01)
#max rmse
for(col in seq(1,9)){
  r <- rmse[,col]
  print(max(r,na.rm=TRUE))
}

#min rmse
for(col in seq(1,9)){
  r <- rmse[,col]
  print(min(r,na.rm=TRUE))
}

