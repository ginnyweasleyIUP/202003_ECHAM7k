simpleawmean<-function(fld, lats = seq(from = -90, to = 90, length.out = 73)){
  zonmean<-apply(fld,2,mean,na.rm=TRUE)
  w.lats<-cos(lats*pi/180)/sum(cos(lats*pi/180))
  #if ((sum(is.na(zonmean))) > 2) {return(NA)}
  tavg<-sum(w.lats*zonmean,na.rm=TRUE)
  return(tavg)
}

simpleawsd<-function(fld, lats = seq(from = -90, to = 90, length.out = 73)){
  zonmean<-apply(fld,2,sd,na.rm=TRUE)
  w.lats<-cos(lats*pi/180)/sum(cos(lats*pi/180))
  #if ((sum(is.na(zonmean))) > 2) {return(NA)}
  tavg<-sum(w.lats*zonmean,na.rm=TRUE)
  return(tavg)
}