cynmix=read.csv("~/Desktop//autoparts/plots/cynmix.csv")

cyn=list(cynmix[1:7,],cynmix[8:14,],cynmix[9:21,],cynmix[22:28,])


layout(matrix(1:(4*7),4,7))
par(mar=c(2,2,2,2))
sapply(1:7, function(x){
	lapply(cyn,function(df){
		plot(1:11,df[x,3:13],pch=c(1,1,1,2,2,2,3,3,3,4,5),col=c(1,1,1,2,2,2,3,3,3,4,5), xaxt="n",cex=1.3
		#, ylim=c(1,3)
		)
		legend("topleft",legend=df[x,14],cex=1.3,bty="n")
	})
})