prismpredict<-function(l,newdat,hvar=NULL){
  l2<-l$prism
  #prism
  if (is.null(hvar)){
    #if the tree only have root node
    if (dim(l2$frame)[1]==1){
      temp<-all.vars(l$formula)
      w<-kmweight.fix(newdat[,which(names(newdat)==temp[1])],newdat[,which(names(newdat)==temp[2])])
      l.new<-lm(newdat[,which(names(newdat)==temp[1])]~as.factor(newdat[,which(names(newdat)==tail(temp,1))]),weights=w)
      testnodes.dat<-cbind(rep(coef(l.new)[1],nrow(newdat)),rep(coef(l.new)[2],nrow(newdat)))
      design.mat<-cbind(rep(1,nrow(newdat)),newdat[,which(names(newdat)==tail(temp,1))])
      l2.test.hadamard<-hadamard.prod(design.mat,testnodes.dat)
      l2.test.fitted<-apply(l2.test.hadamard,1,sum)
    }else{
      #if the tree does not have root node only
      nodes.l2<-l2
      nodes.l2$frame$yval=as.numeric(rownames(nodes.l2$frame))
      testnodes<-predict(nodes.l2,newdat,type="vector")#find which terminal node each new data point falls into
      l2.coef<-l2$frame$yval2[,-1]
      testnodes.dat<-cbind(testnodes,l2.coef[match(testnodes,row.names(l2$frame)),])#find the node specific regression coefficients for each new data point
      design.mat<-cbind(rep(1,nrow(newdat)),newdat[,which(names(newdat)==tail(all.vars(l$formula),1))])
      l2.test.hadamard<-hadamard.prod(design.mat,testnodes.dat[,-1])
      l2.test.fitted<-apply(l2.test.hadamard,1,sum)#predicted value
    }
  } else{
    #hprism
    #if the tree only have root node
    if (dim(l2$frame)[1]==1){
      temp<-all.vars(l$formula)
      w<-kmweight.fix(newdat[,which(names(newdat)==temp[1])],newdat[,which(names(newdat)==temp[2])])
      l.new<-lm(newdat[,which(names(newdat)==temp[1])]~as.factor(newdat[,which(names(newdat)==tail(temp,1))])+newdat[,which(names(newdat)==hvar)]:(as.factor(newdat[,which(names(newdat)==tail(temp,1))])),weights=w)
      l2.test.fitted<-fitted(l.new)
    }else{
      #if the tree does not have root node only
      nodes.l2<-l2
      nodes.l2$frame$yval=as.numeric(rownames(nodes.l2$frame))
      testnodes<-predict(nodes.l2,newdat,type="vector")#find which terminal node each new data point falls into
      l2.coef<-l2$frame$yval2[,-1]
      testnodes.dat<-cbind(testnodes,l2.coef[match(testnodes,row.names(l2$frame)),])#find the node specific regression coefficients for each new data point
      tt<-all.vars(l$formula)
      design.mat<-cbind(rep(1,nrow(newdat)),newdat[,which(names(newdat)==tail(all.vars(l$formula),1))],newdat[,which(names(newdat)==hvar)])
      M<-data.frame(design.mat,testnodes.dat)
      M$int<-ifelse(M[,2]==0,M[,5],M[,5]+M[,6])
      M$slp<-ifelse(M[,2]==0,M[,7],M[,8])
      Design<-cbind(M[,1],M[,3])#design matrix
      Coef<-cbind(M$int,M$slp)#coef matrix
      l2.test.hadamard<-hadamard.prod(Design,Coef)
      l2.test.fitted<-apply(l2.test.hadamard,1,sum)#predicted value
    }}
  return(l2.test.fitted)
}
