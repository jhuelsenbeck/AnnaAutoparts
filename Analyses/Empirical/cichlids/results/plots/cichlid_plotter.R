cichlid=read.csv("~/Desktop//autoparts/plots/cichlid.csv")

cic=list(cichlid[1:5,],cichlid[6:10,],cichlid[11:15,],cichlid[16:20,])


layout(matrix(1:(4*5),4,5))
par(mar=c(2,2,2,2))
sapply(1:5, function(x){
	lapply(cic,function(df){
		plot(1:7,df[x,3:9],pch=c(1,2,2,2,3,3,3),col=c(1,2,3,4,2,3,4), xaxt="n",cex=1.3
		#, ylim=c(1,3)
		)
		legend("topleft",legend=df[x,8],cex=1.3,bty="n")
	})
})