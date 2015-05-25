library(xts)
par(mfrow=c(1,2),
    cex=0.6, cex.axis=1.0, cex.lab=1.0, cex.main=1.5, cex.sub=1.0,
    font.lab=8, lab=c(10,10,10),
    oma=c(2,2,2,2),
    mar=c(5.1,4.1,4.1,2.1))

#1. Wind Forecast plot
cols = c("blue", "red")
dwind = data.frame(Observed=test.data[1:1000,1],Forecast=output[1:1000,1])
ts.plot(dwind,
        main= "Wind Forecast",
        xlab="time",
        ylab="Power Output", #"Output",
        col=cols, ylim=c(0,1.2))
legend("topright", colnames(dwind), col=cols, lty=1, cex=1)

#2. Solar Forecast plot
cols = c("blue", "red")
dsolr = data.frame(Observed=test.data[1:160,1],Forecast=output[1:160,1])
ts.plot(dsolr,
        main= "Solar Forecast",
        xlab="time",
        ylab="Power Output", #"Output",
        col=cols, ylim=c(0,1))#160
legend("topright", colnames(dsolr), col=cols, lty=1, cex=1)


#3. Chapter-4 plots
cols = c("blue", "green")
#cols = c("blue", "red")
d = data.frame(Power=pow.data.set,Wind=wind.data.set)
#d = data.frame(Power=pow.data.set[0:1000],Radiation=solar.data.set[0:1000])

ts.plot(d,
     main= "Data",
     xlab="time",
     ylab="Output",
     col=cols, ylim=c(0,1))
legend("topright", colnames(d), col=cols, lty=1, cex=1, bty="n")
