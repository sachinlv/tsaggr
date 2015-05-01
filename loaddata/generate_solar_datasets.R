details.path <- '/usr/claudio_solar_data/Arizona/powr_rad_site_combi.csv'
pow.data.path <- '/usr/claudio_solar_data/Arizona/DA/DPV/'
rad.data.path <- '/usr/claudio_solar_data/Arizona/radiation/'
output.path <- '/usr/claudio_solar_data/Arizona/prepared/'
output.file.generic <- 'SITE'

details = read.csv(details.path)

for i in seq(1,length(dat[,1])){
 vals <- c(dat[i,])
 setwd(pow.data.path)
 pow.data.file <- vals[1]
 pow.data <- read.csv(pow.data.file)

 setwd(rad.data.path)
 rad.data.file <- paste(vals[4],'_2006_solar.csv', ,sep="")
 rad.data <- read.csv(rad.data.file)
 rad.data <- rad.data[,1:length(rad.data[3,])]

 data <- cbind(pow.data, rad.data)
 setwd(output.path)
 output.file<- paste(output.file.generic, '_', vals[2], '_', vals[3],'.csv',sep="" )
 write.csv(data, output.file)
}
