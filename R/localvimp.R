localvimp<-function(l,hvar=NULL,M=30,dat,cp.prism=0){
  l2<-l$prism
  f<-l$formula
  cpprism<-cp.prism
  temp<-all.vars(f)
  n<-length(temp)-3
  myls<-vector("list",length=n)
  names(myls)<-temp[3:(n+3-1)]
  k<-length(unique(l2$where))
  tt<-cbind(l2$where,l$fitted,dat[,which(names(dat)==temp[1])])
  res<-rep(NA,k)
  t<-1
  for (s in unique(l2$where)){
    d<-tt[which(tt[,1]==s),]
    res[t]<-mean((d[,2]-d[,3])^2)
    t<-t+1
  }
  if (is.null(hvar)){
    for (i in 1:n){
      cat(names(myls)[i],":","\n")
      myls[[i]]<-matrix(NA,nrow=M,ncol=k)
      for (j in 1:M){
        cat("Iteration: ",j,"\n")
        dnew<-dat
        dnew[,which(names(dnew)==names(myls)[i])]<-sample(dnew[,which(names(dnew)==names(myls)[i])])
        l.new<-prism(f,dat=dnew,cp.prism=cpprism)
        nfitted<-l.new$fitted
        tt.new<-cbind(l2$where,nfitted,dnew[,which(names(dnew)==temp[1])],dat[,which(names(dat)==temp[2])])
        m<-1
        for (c in unique(l2$where)){
          d.new<-tt.new[which(tt.new[,1]==c),]
          res<-d.new[,3]-d.new[,2]
          res.nc<-res[which(d.new[, 4]==1)]
          res.c<-res[which(d.new[, 4]==0)]
          res[which(d[, 4]==0)]<-sapply(res.c,function(rr){mean(res.nc[which(res.nc>rr)])})
          myls[[i]][j,m]<-mean(res^2,na.rm=T)
          m<-m+1
        }
      }
    }
  }else{
    h.var<-hvar
    for (i in 1:n){
      cat(names(myls)[i],":","\n")
      myls[[i]]<-matrix(NA,nrow=M,ncol=k)
      for (j in 1:M){
        cat("Iteration: ",j,"\n")
        dnew<-dat
        dnew[,which(names(dnew)==names(myls)[i])]<-sample(dnew[,which(names(dnew)==names(myls)[i])])
        l.new<-prism(f,dat=dnew,cp.prism=cpprism,hvar=h.var)
        nfitted<-l.new$fitted
        tt.new<-cbind(l2$where,nfitted,dnew[,which(names(dnew)==temp[1])],dat[,which(names(dat)==temp[2])])
        m<-1
        for (c in unique(l2$where)){
          d.new<-tt.new[which(tt.new[,1]==c),]
          res<-d.new[,3]-d.new[,2]
          res.nc<-res[which(d.new[, 4]==1)]
          res.c<-res[which(d.new[, 4]==0)]
          res[which(d[, 4]==0)]<-sapply(res.c,function(rr){mean(res.nc[which(res.nc>rr)])})
          myls[[i]][j,m]<-mean(res^2,na.rm=T)
          m<-m+1
        }
      }
    }
  }
  #myls.m<-myls.raw<-myls.rank<-vector("list",length=n)
  #names(myls.m)<-names(myls.raw)<-names(myls.rank)<-names(myls)
  myls.m<-vector("list",length=n)
  names(myls.m)<-names(myls)
  myls.rawvalue<-myls.rank<-matrix(NA,nrow=n,ncol=k)
  for (i in 1:n){
    myls.m[[i]]<-colMeans(myls[[i]])
    #myls.raw[[i]]<-myls.m[[i]]-res
    myls.rawvalue[i,]<-myls.m[[i]]-res
    #myls.rank[[i]]<-rank(-(myls.m[[i]]-res),ties.method="min")
  }
  for (v in 1:k){
    myls.rank[,v]<-rank(-(myls.rawvalue[,v]),ties.method="min")
  }
  rownames(myls.rank)<-rownames(myls.rawvalue)<-names(myls.m)
  colnames(myls.rank)<-colnames(myls.rawvalue)<-rownames(l2$frame[unique(l2$where),])
  lvimp<-NULL
  lvimp$rawvalue<-myls.rawvalue
  lvimp$rank<-myls.rank
  #lvimp$mean<-myls.m
  return(lvimp)
}
