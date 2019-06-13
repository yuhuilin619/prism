#evaluation, splitting and initialization function for prism with hierarchical interaction
#for generating tree plot
hploteval<-function (y, wt,parms)
{
  y0<-which((y[,2]==1)&(y[,4]==0))
  y1<-which((y[,2]==1)&(y[,4]==1))
  w<-kmweight.fix(y[, 1],y[, 2])
  tfit<-lm(y[,1]~as.factor(y[,4])+y[,3]:as.factor(y[,4]),weights=w)
  tfit.coef <- coefficients(tfit)
  tfit.res<-rep(NA,dim(y)[1])
  res.nc<-residuals(tfit)[which(y[, 2]==1)]
  res.c<-residuals(tfit)[which(y[, 2]==0)]
  tfit.res[which(y[, 2]==1)]<-res.nc
  tfit.res[which(y[, 2]==0)]<-sapply(res.c,function(rr){mean(res.nc[which(res.nc>rr)])})
  if (sum(y[, 2]==0)==0) {tfit.res<-res.nc}
  tfit.gof <-sum(tfit.res^2,na.rm=T)
  list(label = c( tfit.coef[1], tfit.coef[2],tfit.coef[3],tfit.coef[4]), deviance = tfit.gof)
}

hplotsplit<-function(y, wt, x, parms, continuous)
{
  n <- length(y)
  w<-kmweight.fix(y[, 1],y[, 2])
  lasso.parent <- lm(y[,1]~as.factor(y[,4])+y[,3]:as.factor(y[,4]),weights=w)
  tfit.res<-rep(NA,dim(y)[1])
  res.nc<-residuals(lasso.parent)[which(y[, 2]==1)]
  res.c<-residuals(lasso.parent)[which(y[, 2]==0)]
  tfit.res[which(y[, 2]==1)]<-res.nc
  tfit.res[which(y[, 2]==0)]<-sapply(res.c,function(rr){mean(res.nc[which(res.nc>rr)])})
  if (sum(y[, 2]==0)==0) {tfit.res<-res.nc}
  lasso.p.gof <- sum(tfit.res^2,na.rm=T)
  #continous splitting variable
  if (continuous) {
    n <- nrow(y)
    goodness <- double(n - 1)
    direction <- goodness
    temp <- rep(0, n)
    for (i in 1:(n - 1)) {
      if (x[i] != x[i + 1]) {
        Left<-which(x<=x[i])
        #examine inclusion and exclusion criteria for a potential split
        if ((length(Left)>1)&(length(Left)<n-1)){
          yLeft<-y[which(x<=x[i]),]
          yRight<-y[which(x>x[i]),]
          if ((length(unique(yLeft[,4]))>1)&(length(unique(yRight[,4]))>1)){
            yLeft.ind0<-which((yLeft[,2]==1)&(yLeft[,4]==0))
            yLeft.ind1<-which((yLeft[,2]==1)&(yLeft[,4]==1))
            yRight.ind0<-which((yRight[,2]==1)&(yRight[,4]==0))
            yRight.ind1<-which((yRight[,2]==1)&(yRight[,4]==1))
            if ((length(yLeft.ind0)>1)&(length(yLeft.ind1)>1)&(length(yRight.ind0)>1)&(length(yRight.ind1)>1)){
              if ((length(table(yLeft[yLeft.ind0,3]))>1)&(length(table(yLeft[yLeft.ind1,3]))>1)&(length(table(yRight[yRight.ind0,3]))>1)&(length(table(yRight[yRight.ind1,3]))>1)){
                wLeft<-kmweight.fix(yLeft[, 1],yLeft[, 2])
                wRight<-kmweight.fix(yRight[, 1],yRight[, 2])
                tfit1<-lm(yLeft[,1]~as.factor(yLeft[,4])+yLeft[,3]:as.factor(yLeft[,4]),weights=wLeft)
                tfit2<-lm(yRight[,1]~as.factor(yRight[,4])+yRight[,3]:as.factor(yRight[,4]),weights=wRight)
                tfit1.res<-rep(NA,dim(yLeft)[1])
 	            tfit2.res<-rep(NA,dim(yRight)[1])
 	            res1.nc<-residuals(tfit1)[which(yLeft[, 2]==1)]
 	            res1.c<-residuals(tfit1)[which(yLeft[, 2]==0)]
 	            res2.nc<-residuals(tfit2)[which(yRight[, 2]==1)]
 	            res2.c<-residuals(tfit2)[which(yRight[, 2]==0)]
 	            tfit1.res[which(yLeft[, 2]==1)]<-res1.nc
 	            tfit1.res[which(yLeft[, 2]==0)]<-sapply(res1.c,function(rr){mean(res1.nc[which(res1.nc>rr)])})
 	            tfit2.res[which(yRight[, 2]==1)]<-res2.nc
 	            tfit2.res[which(yRight[, 2]==0)]<-sapply(res2.c,function(rr){mean(res2.nc[which(res2.nc>rr)])})
 	            if (sum(yLeft[, 2]==0)==0) {tfit1.res<-res1.nc}
 	            if (sum(yRight[, 2]==0)==0) {tfit2.res<-res2.nc}
                tfit.gof <- sum(tfit1.res^2,na.rm=T)+sum(tfit2.res^2,na.rm=T)
                goodness[i] <- max(lasso.p.gof - tfit.gof,0)
              }}}}else{goodness[i] <- 0}
        direction[i] <- -1
      }
    }
  }
  #categorical splitting variable
  else {
    w<-kmweight.fix(y[, 1],y[, 2])
    tfit.cp <- lm(y[,1]~as.factor(y[,4])+y[,3]:as.factor(y[,4]),weights=w)
    tfit.res<-rep(NA,dim(y)[1])
    res.nc<-residuals(tfit.cp)[which(y[, 2]==1)]
 	res.c<-residuals(tfit.cp)[which(y[, 2]==0)]
 	tfit.res[which(y[, 2]==1)]<-res.nc
 	tfit.res[which(y[, 2]==0)]<-sapply(res.c,function(rr){mean(res.nc[which(res.nc>rr)])})
 	if (sum(y[, 2]==0)==0) {tfit.res<-res.nc}
    tfit.cp.gof<-sum(tfit.res^2,na.rm=T)
    tfit <- lm(y[,1]~factor(x) - 1)
    ngrp <- length(tfit$coef)
    direction <- levels(factor(x))
    xx <- direction[match(x, sort(unique(x)))]
    goodness <- double(length(direction) - 1)
    xp<-sort(unique(xx))
    for (i in 1:length(goodness)) {
      Left<-which(xx<=xp[i])
      Right<-which(xx>xp[i])
      yLeft<-y[Left,]
      yRight<-y[Right,]
      #examine inclusion and exclusion criteria for a potential split
      if ((length(Left)>1)&(length(Right)>1)){
        if ((length(unique(yLeft[,4]))>1)&(length(unique(yRight[,4]))>1)){
          yLeft.ind0<-which((yLeft[,2]==1)&(yLeft[,4]==0))
          yLeft.ind1<-which((yLeft[,2]==1)&(yLeft[,4]==1))
          yRight.ind0<-which((yRight[,2]==1)&(yRight[,4]==0))
          yRight.ind1<-which((yRight[,2]==1)&(yRight[,4]==1))
          if ((length(yLeft.ind0)>1)&(length(yLeft.ind1)>1)&(length(yRight.ind0)>1)&(length(yRight.ind1)>1)){
            if ((length(table(yLeft[yLeft.ind0,3]))>1)&(length(table(yLeft[yLeft.ind1,3]))>1)&(length(table(yRight[yRight.ind0,3]))>1)&(length(table(yRight[yRight.ind1,3]))>1)){
            wLeft<-kmweight.fix(yLeft[, 1],yLeft[, 2])
            wRight<-kmweight.fix(yRight[, 1],yRight[, 2])
            tfit1<-lm(yLeft[,1]~as.factor(yLeft[,4])+yLeft[,3]:as.factor(yLeft[,4]),weights=wLeft)
            tfit2<-lm(yRight[,1]~as.factor(yRight[,4])+yRight[,3]:as.factor(yRight[,4]),weights=wRight)
            res1.nc<-residuals(tfit1)[which(yLeft[, 2]==1)]
 	        res1.c<-residuals(tfit1)[which(yLeft[, 2]==0)]
 	        res2.nc<-residuals(tfit2)[which(yRight[, 2]==1)]
 	        res2.c<-residuals(tfit2)[which(yRight[, 2]==0)]
 	        tfit1.res[which(yLeft[, 2]==1)]<-res1.nc
 	        tfit1.res[which(yLeft[, 2]==0)]<-sapply(res1.c,function(rr){mean(res1.nc[which(res1.nc>rr)])})
 	        tfit2.res[which(yRight[, 2]==1)]<-res2.nc
 	        tfit2.res[which(yRight[, 2]==0)]<-sapply(res2.c,function(rr){mean(res2.nc[which(res2.nc>rr)])})
 	        if (sum(yLeft[, 2]==0)==0) {tfit1.res<-res1.nc}
 	        if (sum(yRight[, 2]==0)==0) {tfit2.res<-res2.nc}
            tfit.gof <- sum(tfit1.res^2,na.rm=T)+sum(tfit2.res^2,na.rm=T)
            goodness[i]<-max(tfit.cp.gof - tfit.gof,0)
            }}}}else{goodness[i] <- 0}
    }
  }
  list(goodness = goodness, direction = direction)
}


hplotinit<-function (y, offset, parms=0, wt)
{
  if (is.null(offset))
    offset <- 0
  sfun <- function(yval, dev, wt, ylevel, digits) {
    paste("fint= ", format(signif(yval[,1], digits)), "\n",
          "fslope=", format(signif(yval[, 2],digits)), "\n",
          "finteraction0=", format(signif(yval[, 3],digits)), "\n",
          "finteraction1=", format(signif(yval[, 4],digits)), "\n",
          "deviance=", format(signif(dev, digits)), sep = "")
  }
  tfun <- function(yval, dev, wt, ylevel, digits, n, use.n) {
    if (use.n) {
      paste("fint= ", format(signif(yval[, 1], digits)),  "\n",
            "fslope=", format(yval[, 2], digits),"\n",
            ",finteraction0=", format(signif(yval[, 3],digits)), "\n",
            "finteraction1=", format(signif(yval[, 4],digits)), "\n",
            "Deviance=", format(dev, digits), "\nn=", n, sep = "")
    }
    else {
      paste("fint= ", format(signif(yval[, 1], digits)), "\n",
            "fslope=", format(yval[, 2], digits), "\n",
            "finteraction0=", format(signif(yval[, 3],digits)),"\n",
            "finteraction1=", format(signif(yval[, 4],digits)), "\n",
            "Deviance=", format(dev, digits))
    }
  }
  environment(sfun) <- .GlobalEnv
  environment(tfun) <- .GlobalEnv
  list(y = cbind(y, parms,offset),parms=parms,numresp = 4, numy = 4,
       summary = sfun, text = tfun)
}
ulist.hprismplot<-list(eval=hploteval,split=hplotsplit,init=hplotinit)

