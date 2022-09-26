hummer=read.csv("~/Desktop//autoparts/plots/hummer.csv")

hb=list(hummer[1:6,],hummer[7:12,],hummer[13:18,],hummer[19:24,])


layout(matrix(1:(4*6),4,6))
par(mar=c(2,2,2,2))
sapply(1:6, function(x){
	lapply(hb,function(df){
		plot(1:9,df[x,3:11],pch=c(1,2,2,2,3,4,5,5,5),col=c(1,2,3,4,5,6,2,3,4), xaxt="n",cex=1.3
		#, ylim=c(1,3)
		)
		legend("topleft",legend=df[x,10],cex=1.3,bty="n")
	})
})

