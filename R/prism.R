prism<-function(formula,dat,hvar=NULL,cp.prism=0){
  f<-as.formula(formula)
  prism.result<-NULL
  if (is.null(hvar)){
    #prism
    l<-rpart(f,data=dat,method=ulist.prism,cp=cp.prism)
    nodes.l<-l
    nodes.l$frame$yval=as.numeric(rownames(nodes.l$frame))
    testnodes<-predict(nodes.l,dat,type="vector")
    l.coef<-l$frame$yval2[,-1]
    testnodes.dat<-cbind(testnodes,l.coef[match(testnodes,row.names(l$frame)),])
    design.mat<-cbind(rep(1,nrow(dat)),dat[,which(names(dat)==tail(all.vars(f),1))])
    l.test.hadamard<-hadamard.prod(design.mat,testnodes.dat[,-1])
    l.test.fitted<-apply(l.test.hadamard,1,sum)
    prism.result$formula<-f
    prism.result$h.variable<-NULL
    prism.result$prism<-l
    prism.result$prism.plot<-l
    prism.result$fitted<-l.test.fitted
  } else{
    #hprism
    l<-rpart(f,data=dat,method=ulist.hprism,parms=dat[,which(names(dat)==hvar)],cp=cp.prism)
    l.plot<-rpart(f,data=dat,method=ulist.hprismplot,parms=dat[,which(names(dat)==hvar)],cp=cp.prism)
    nodes.l<-l
    nodes.l$frame$yval=as.numeric(rownames(nodes.l$frame))
    testnodes<-predict(nodes.l,dat,type="vector")#find which terminal node each data point falls into
    l.coef<-l$frame$yval2[,-1]
    testnodes.dat<-cbind(testnodes,l.coef[match(testnodes,row.names(l$frame)),])#find the node specific regression coefficients for each data point
    design.mat<-cbind(rep(1,nrow(dat)),dat[,which(names(dat)==tail(all.vars(f),1))],dat[,which(names(dat)==hvar)])
    M<-data.frame(design.mat,testnodes.dat)
    M$int<-ifelse(M[,2]==0,M[,5],M[,5]+M[,6])
    M$slp<-ifelse(M[,2]==0,M[,7],M[,8])
    Design<-cbind(M[,1],M[,3])#design matrix
    Coef<-cbind(M$int,M$slp)#coef matrix
    l.test.hadamard<-hadamard.prod(Design,Coef)
    l.test.fitted<-apply(l.test.hadamard,1,sum)#fitted value
    prism.result$formula<-f
    prism.result$h.variable<-hvar
    prism.result$prism<-l
    prism.result$prism.plot<-l.plot
    prism.result$fitted<-l.test.fitted
  }
  return(prism.result)
}
