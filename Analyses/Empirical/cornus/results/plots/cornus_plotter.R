cornus=read.csv("~/Desktop//autoparts/plots/cornus.csv")

cor=list(cornus[1:3,],cornus[4:6,],cornus[7:9,],cornus[10:12,])


layout(matrix(1:(4*5),4,5))
par(mar=c(2,2,2,2))
sapply(1:3, function(x){
	lapply(cor,function(df){
		plot(1:4,df[x,3:6],pch=c(1,2,2,2),col=c(1,2,3,4), xaxt="n",cex=1.3
		#, ylim=c(1,3)
		)
		legend("topleft",legend=df[x,7],cex=1.3,bty="n")
	})
})