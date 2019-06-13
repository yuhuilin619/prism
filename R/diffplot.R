diffplot<-function(l,dat,hvar,dot.size=1,legend.size=0.75){
  l2<-l$prism
  nodes.l<-l2
  nodes.l$frame$yval=as.numeric(rownames(nodes.l$frame))
  testnodes<-predict(nodes.l,dat,type="vector")
  l.coef<-l2$frame$yval2[,-1]
  testnodes.dat<-cbind(testnodes,l.coef[match(testnodes,row.names(l2$frame)),])
  dat.temp<-dat
  dat.temp$node<-testnodes.dat[,1]
  dat.temp$int0<-testnodes.dat[,2]
  dat.temp$int1<-testnodes.dat[,3]
  dat.temp$slp0<-testnodes.dat[,4]
  dat.temp$slp1<-testnodes.dat[,5]
  dat.temp$survival0<-dat.temp$int0+dat.temp$slp0*dat.temp[,which(names(dat.temp)==hvar)]
  dat.temp$survival1<-dat.temp$int0+dat.temp$int1+dat.temp$slp1*dat.temp[,which(names(dat.temp)==hvar)]
  dat.temp$surv.diff<-dat.temp$survival0-dat.temp$survival1
  n<-length(unique(dat.temp$node))
  myls<-vector("list",length=n)
  for (i in 1:n){
    myls[[i]]<-dat.temp[which(dat.temp$node==sort(unique(dat.temp$node))[i]),]
    myls[[i]]<-myls[[i]][,which((names(myls[[i]])==hvar)|(names(myls[[i]])=="surv.diff"))]
    myls[[i]]<-unique(myls[[i]])
  }
  cl<-rainbow(n)
  plot(myls[[1]][,1],myls[[1]][,2],cex=dot.size,col=cl[1],xlab=hvar,ylab="diff in log survival (log(day))",ylim=range(dat.temp$surv.diff))
  if (n>1){
    for (j in 2:n){
      points(myls[[j]][,1],myls[[j]][,2],cex=dot.size,col=cl[j])
    }
  }
  legend("topright",title="PRISM terminal node number", as.character(unique(dat.temp$node)) , lty=1, col=cl, bty='n', cex=legend.size)
}
