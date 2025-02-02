library(xgboost)
library(caret)
library(scales)
forecast.xg <- function (trD, h = h, lags = freq*2){
  n<-ncol(trD)
  td_len<-nrow(trD)
  scale_para<-matrix(NA,nrow=2,ncol=n);scaled_trD<-trD*0
  for (i in 1:n){
    scale_para[1,i]<-mean(trD[,i])
    scale_para[2,i]<-sd(trD[,i])
    scaled_trD[,i]<-(trD[,i]-scale_para[1,i]) / scale_para[2,i]
  }
  X<-matrix(NA,ncol=lags,nrow=n*(td_len-lags))
  Y<-matrix(NA,ncol=1   ,nrow=n*(td_len-lags))
  for (i in 1:n){
    for (t in 1:(td_len-lags)){
      #cat(t,' ',t+lags-1,' ',t+lags,'\n')
      X[(i-1)*(td_len-lags)+t,]<-scaled_trD[t:(t+lags-1),i]
      Y[(i-1)*(td_len-lags)+t,]<-scaled_trD[t+lags,i]
    }
  }
  data_matrix <- xgb.DMatrix(data = X, label = Y)
  
  params <- list(
    objective = "reg:squarederror",
    eval_metric = "rmse"
  )
  model <- xgb.train(params = params, data = data_matrix, nrounds = 100)
  scaled_fit<-predict(model,data_matrix)
  scaled_fit<-matrix(scaled_fit,ncol=n)
  fit<-scaled_fit
  for (i in 1:n){
    fit[,i]<-scaled_fit[,i]*scale_para[2,i]+scale_para[1,i]
  }
  realvalue<-trD[(lags+1):td_len,]
  res<-realvalue-fit
  
  predict<-matrix(NA,nrow=0,ncol=n)#predict value
  XforPredict<-t(scaled_trD[(td_len-lags+1):td_len,])
  for (i in 1:h){
    Xi<-xgb.DMatrix(data = XforPredict )
    Xnext<-predict(model,Xi)
    predict<-rbind(predict,Xnext)
    XforPredict<-unname(cbind(XforPredict[,-1],Xnext))
  }
  for (i in 1:n){
    predict[,i]<-predict[,i]*scale_para[2,i]+scale_para[1,i]
  }
  
  return (list(fit=fit, residual=res, predict=predict))
}
