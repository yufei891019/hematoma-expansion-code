library(car)
library(ggplot2)
library(glmnet)
library(MASS)
library(survival)
library(rms)
library(survminer)
library(ggridges)
library(pROC)
library(plotROC)
library(riskRegression)
library(foreign)
library(rms)
library(epiDisplay)
library(dplyr)
library(caret)

getwd()
setwd("E:/Python05/FZ168/ER/返回20240130/Code and Data Disclosure")
data<-read.csv("./data.csv")
data.train <- data[data$train_validation.test==1,]
data.validation <- data[data$train_validation.test==2,]
data.tset <- data[data$train_validation.test==3,]

Y.train <- data.train$HE
Y.validation <- data.validation$HE
Y.test <- data.tset$HE

#radiomics data
radiomics.train <- data.train[,c(2:31)]
radiomics.validation <- data.validation[,c(2:31)]
radiomics.test <- data.tset[,c(2:31)]
#Clinical data
clinical.train <- data.train[,c(2,32:37)]
clinical.validation <- data.validation[,c(2,32:37)]
clinical.test <- data.tset[,c(2,32:37)]

###radiomics model##
LG.ML<-glm(HE~.,family="binomial",data = radiomics.train)
#radiomics.train
radiomics.train_LG_p<-predict(LG.ML,data=radiomics.train,type = "response")
roc_radiomics.train_LG<-roc(radiomics.train$HE,radiomics.train_LG_p,ci=T,auc=T)
roc_radiomics.train_LG
plot(roc_radiomics.train_LG,legacy.axes=TRUE)
#radiomics.validation
radiomics.validation_LG_p<-predict(LG.ML,newdata =radiomics.validation,type = "response")
roc_radiomics.validation_LG<-roc(radiomics.validation$HE,radiomics.validation_LG_p,ci=T,auc=T)
roc_radiomics.validation_LG
plot(roc_radiomics.validation_LG,legacy.axes=TRUE)
#radiomics.test
radiomics.test_LG_p<-predict(LG.ML,newdata =radiomics.test,type = "response")
roc_radiomics.test_LG<-roc(radiomics.test$HE,radiomics.test_LG_p,ci=T,auc=T)
roc_radiomics.test_LG
plot(roc_radiomics.test_LG,legacy.axes=TRUE)
#radiomics score
RadScore.train <- log(radiomics.train_LG_p/(1-radiomics.train_LG_p))
RadScore.valldation <- log(radiomics.validation_LG_p/(1-radiomics.validation_LG_p))
RadScore.test <- log(radiomics.test_LG_p/(1-radiomics.test_LG_p))
#################################### clinical imaging sign model  ###########################
#logistics
fit.log <- glm(HE~ Time+Volume+Irregular + BlendSign  +SmokeHistory+ GCS,family = binomial(),data = clinical.train)
summary(fit.log)
coefficients(fit.log)
exp(coefficients(fit.log))
exp (confint(fit.log))
#clinical.train
train_LC.P <- predict(fit.log, clinical.train,type="response")
train_roc_LC<-roc(clinical.train$HE,train_LC.P,ci=T,auc=T)
train_roc_LC
plot(train_roc_LC,legacy.axes=TRUE)
#clinical.validation
val_LC.P <- predict(fit.log, clinical.validation,type="response")
val_roc_LC<-roc(clinical.validation$HE,val_LC.P,ci=T,auc=T)
val_roc_LC
plot(val_roc_LC,legacy.axes=TRUE)
#clinical.test
test_LC.P <- predict(fit.log, clinical.test,type="response")
test_roc_LC<-roc(clinical.test$HE,test_LC.P,ci=T,auc=T)
test_roc_LC
plot(test_roc_LC,legacy.axes=TRUE)


train <- cbind(clinical.train,RadScore.train)
names(train)[8] <- "Rad.score"
validation <- cbind(clinical.validation,RadScore.valldation)
names(validation)[8] <- "Rad.score"
test <- cbind(clinical.test,RadScore.test)
names(test)[8] <- "Rad.score"
###################### hybrid model #######################
fit.hybrid.1 <- glm(HE~ Time+Volume+Irregular + BlendSign  +SmokeHistory+ GCS+Rad.score,family = binomial(),data = train)
summary(fit.hybrid.1)
### The p-values of Time, BlendSign, and SmokeHistory are greater than 0.05, indicating a small contribution to the model. These three variables are excluded.
fit.hybrid<-glm(HE~GCS+Volume+Irregular + Rad.score,family = binomial(),data = train)
summary(fit.hybrid)
coefficients(fit.hybrid)
exp(coefficients(fit.hybrid))
exp (confint(fit.hybrid))

### nomograph ###
dd<-datadist(train)
options(datadist='dd')
fit_nomograph<- lrm(HE~GCS+Volume+Irregular + Rad.score,data = train,x=T,y=T)
fit_nomograph

nom1 <- nomogram(fit_nomograph, fun=plogis,fun.at=c(.001, .01,0.1,seq(.2,.8, by=.2), .9, .99,.999), 
                 lp=F, funlabel="Risk of developing HE")
plot(nom1,cex.axis=1.0,cex.var=1.5,xfrac=0.35)

##radiomics model##
fit_R.score<-glm(HE ~ Rad.score,family = binomial(),data = train)
summary(fit_R.score)

################################Calculate probability################################################
#clinical imaging sign model
train$predvalue_clinic_train<-predict(fit.log, train,type="response")
validation$predvalue_clinic_validation <- predict(fit.log, validation,type="response")
test$predvalue_clinic_test<-predict(fit.log, test,type="response")
##radiomics model###
train$predvalue_R.score_train<-predict(fit_R.score, train,type="response")
validation$predvalue_R.score_val<-predict(fit_R.score, validation,type="response")
test$predvalue_R.score_test<-predict(fit_R.score, test,type="response")
##hybrid model###
train$predvalue_hybrid_train<-predict(fit.hybrid,train,type="response")
validation$predvalue_hybrid_val<-predict(fit.hybrid,validation,type="response")
test$predvalue_hybrid_test<-predict(fit.hybrid,test,type="response")



####AUC####

AUC_rad.score_train<-roc(train$HE,train$Rad.score,data=train,ci=T,auc=T)
AUC_rad.score_train
AUC_rad.score_val<-roc(HE~Rad.score,data=validation,ci=T,auc=T)
AUC_rad.score_val
AUC_rad.score_test<-roc(HE~Rad.score,data=test,ci=T,auc=T)
AUC_rad.score_test


AUC_clinic_train<-roc(HE~predvalue_clinic_train,data=train,ci=T,auc=T)
AUC_clinic_train
AUC_clinic_val<-roc(HE~predvalue_clinic_validation,data=validation,ci=T,auc=T)
AUC_clinic_val
AUC_clinic_test<-roc(HE~predvalue_clinic_test,data=test,ci=T,auc=T)
AUC_clinic_test


AUC_hybrid_train<-roc(HE~predvalue_hybrid_train,data=train,ci=T,auc=T)
AUC_hybrid_train
AUC_hybrid_val<-roc(HE~predvalue_hybrid_val,data=validation,ci=T,auc=T)
AUC_hybrid_val
AUC_hybrid_test<-roc(HE~predvalue_hybrid_test,data=test,ci=T,auc=T)
AUC_hybrid_test

###confusion matrix####
library(pROC)
Actual_test<-factor(test$HE,levels = c(1,0),labels = c("TRUE","FALSE"))
Actual_val<-factor(validation$HE,levels = c(1,0),labels = c("TRUE","FALSE"))


#clinical imaging sign model#
# train
roc_clinic_train<-roc(train$HE,train$predvalue_clinic_train)
ROC.clinic_train<-coords(roc_clinic_train,"best",ret="all",transpose=FALSE)
as.matrix(ROC.clinic_train)
#validation#
Predict_clinic_val_1_0<-ifelse(validation$predvalue_clinic_validation>0.3664429,1,0)
Predict_clinic_val_1_0<-factor(Predict_clinic_val_1_0,levels = c(1,0),labels = c("TRUE","FALSE"))
xtab_clinic_val<-table(Predict_clinic_val_1_0,Actual_val)
xtab_clinic_val
confusionMatrix(xtab_clinic_val)
#test#
Predict_clinic_test_1_0<-ifelse(test$predvalue_clinic_test>0.3664429,1,0)
Predict_clinic_test_1_0<-factor(Predict_clinic_test_1_0,levels = c(1,0),labels = c("TRUE","FALSE"))
xtab_clinic_test<-table(Predict_clinic_test_1_0,Actual_test)
xtab_clinic_test
confusionMatrix(xtab_clinic_test)

##radiomics model
#train#
roc_rad.score_train<-roc(train$HE,train$predvalue_R.score_train)
ROC.rad.score_train<-coords(roc_rad.score_train,"best",ret="all",transpose=FALSE)
as.matrix(ROC.rad.score_train)
#val#
Predict_rad.score_val_1_0<-ifelse(validation$predvalue_R.score_val>0.3004192,1,0)
Predict_rad.score_val_1_0<-factor(Predict_rad.score_val_1_0,levels = c(1,0),labels = c("TRUE","FALSE"))
xtab_rad.score_val<-table(Predict_rad.score_val_1_0,Actual_val)
xtab_rad.score_val
confusionMatrix(xtab_rad.score_val)
#test#
Predict_rad.score_test_1_0<-ifelse(test$predvalue_R.score_test>0.3004192,1,0)
Predict_rad.score_test_1_0<-factor(Predict_rad.score_test_1_0,levels = c(1,0),labels = c("TRUE","FALSE"))
xtab_rad.score_test<-table(Predict_rad.score_test_1_0,Actual_test)
xtab_rad.score_test
confusionMatrix(xtab_rad.score_test)

##hybrid model
#train#
roc_hybrid_train<-roc(train$HE,train$predvalue_hybrid_train)
ROC.hybrid_train<-coords(roc_hybrid_train,"best",ret="all",transpose=FALSE)
as.matrix(ROC.hybrid_train)
#val#
Predict_hybrid_val_1_0<-ifelse(validation$predvalue_hybrid_val>0.3861082,1,0)
Predict_hybrid_val_1_0<-factor(Predict_hybrid_val_1_0,levels = c(1,0),labels = c("TRUE","FALSE"))
xtab_hybrid_val<-table(Predict_hybrid_val_1_0,Actual_val)
xtab_hybrid_val
confusionMatrix(xtab_hybrid_val)
#test#
Predict_hybrid_test_1_0<-ifelse(test$predvalue_hybrid_test>0.3861082,1,0)
Predict_hybrid_test_1_0<-factor(Predict_hybrid_test_1_0,levels = c(1,0),labels = c("TRUE","FALSE"))
xtab_hybrid_test<-table(Predict_hybrid_test_1_0,Actual_test)
xtab_hybrid_test
confusionMatrix(xtab_hybrid_test)
############################################   ROC ###################################
####train
library(ROCR)
pred_R.score_train<-prediction(train$predvalue_R.score_train,train$HE)
perf_R.score_train<-performance(pred_R.score_train,"tpr","fpr")
plot(perf_R.score_train,xlab="1-Specificity",ylab="Sensitivity",col="green",lwd=2)
abline(0,1,col="black",lty=2,lwd=2)

pred_clinic_train<-prediction(train$predvalue_clinic_train,train$HE)
perf_clinic_train<-performance(pred_clinic_train,"tpr","fpr")
plot(perf_clinic_train,col="blue",add=TRUE,lwd=2)

pred_hybrid_train<-prediction(train$predvalue_hybrid_train,train$HE)
perf_hybrid_train<-performance(pred_hybrid_train,"tpr","fpr")
plot(perf_hybrid_train,col="red",add=TRUE,lwd=2)

legend(0.38,0.12,c("Radiomics model(AUC=0.885)","Clinical imaging sign model(AUC=0.759)","Hybrid model(AUC=0.901)"),lty = c(1,1,1),lwd=c(2,2,2,2),col=c("green","blue","red"),bty = "n")

####validation
pred_R.score_val<-prediction(validation$predvalue_R.score_val,validation$HE)
perf_R.score_val<-performance(pred_R.score_val,"tpr","fpr")
plot(perf_R.score_val,xlab="1-Specificity",ylab="Sensitivity",col="green",lwd=2)
abline(0,1,col=1,lty=2,lwd=2)

pred_clinic_val<-prediction(validation$predvalue_clinic_val,validation$HE)
perf_clinic_val<-performance(pred_clinic_val,"tpr","fpr")
plot(perf_clinic_val,col="blue",add=TRUE,lwd=2)

pred_hybrid_val<-prediction(validation$predvalue_hybrid_val,validation$HE)
perf_hybrid_val<-performance(pred_hybrid_val,"tpr","fpr")
plot(perf_hybrid_val,col="red",add=TRUE,lwd=2)

legend(0.38,0.12,c("Radiomics model(AUC=0.827)","Clinical imaging sign model(AUC=0.725)","Hybrid model(AUC=0.838)"),lty = c(1,1,1),lwd=c(2,2,2,2),col=c("green","blue","red"),bty = "n")

####test
pred_R.score_test<-prediction(test$predvalue_R.score_test,test$HE)
perf_R.score_test<-performance(pred_R.score_test,"tpr","fpr")
plot(perf_R.score_test,xlab="1-Specificity",ylab="Sensitivity",col="green",lwd=2)
abline(0,1,col=1,lty=2,lwd=2)

pred_clinic_test<-prediction(test$predvalue_clinic_test,test$HE)
perf_clinic_test<-performance(pred_clinic_test,"tpr","fpr")
plot(perf_clinic_test,col="blue",add=TRUE,lwd=2)

pred_hybrid_test<-prediction(test$predvalue_hybrid_test,test$HE)
perf_hybrid_test<-performance(pred_hybrid_test,"tpr","fpr")
plot(perf_hybrid_test,col="red",add=TRUE,lwd=2)
legend(0.38,0.12,c("Radiomics model(AUC=0.894)","Clinical imaging sign model(AUC=0.765)","Hybrid model(AUC=0.917)"),lty = c(1,1,1),lwd=c(2,2,2,2),col=c("green","blue","red"),bty = "n")


###################################calibration curve##########################################################

library(rms)
library(riskRegression)

dd<-datadist(train)
options(datadist='dd')

#train
fit_clinic_train_lrm<-lrm(HE~Time+Volume+Irregular + BlendSign  +SmokeHistory+ GCS,data = train,x=T,y=T)
cal_clinic_train_lrm<- calibrate(fit_clinic_train_lrm, cmethod='hare',method='boot',B=1000)

fit_r.score_train_lrm<-lrm(HE ~ Rad.score,data = train,x=T,y=T)
cal_r.score_train_lrm <- calibrate(fit_r.score_train_lrm, method='boot',B=1000)

fit_hybrid_train_lrm<-lrm(HE~ GCS+Volume+Irregular + Rad.score,data = train,x=T,y=T)
cal_hybrid_train_lrm <- calibrate(fit_hybrid_train_lrm, method='boot',B=1000,data=train)

plot(1,type="n",xlim=c(0,1),ylim=c(0,1),xlab = "Predicted probability",ylab = "Actual probability",legend=FALSE,subtitles=FALSE)
abline(0,1,col="black",lty=2,lwd=2)
lines(cal_r.score_train_lrm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="green",pch=16)
lines(cal_clinic_train_lrm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="blue",pch=16)
lines(cal_hybrid_train_lrm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="red",pch=16)
legend(0.55,0.18,c("Radiomics model","Clinical imaging sign model","Hybrid model","Ideal model"),lty = c(1,1,1,2),lwd=c(2,2,2,2),col=c("green","blue","red","black"),bty = "n")


#val
cal_clinic_val_pro<-predict(object=fit_clinic_train_lrm,type="fitted",newdata=validation)
fit_clinic_val<-lrm(HE~cal_clinic_val_pro,data=validation,x=T,y=T)
cal_clinic_val_lrm<- calibrate(fit_clinic_val, method='boot',B=1000,data=validation)

cal_R.score_val_pro<-predict(object=fit_r.score_train_lrm,type="fitted",newdata=validation)
fit_R.score_val<-lrm(HE~cal_R.score_val_pro,data=validation,x=T,y=T)
cal_R.score_val_lrm<- calibrate(fit_R.score_val, method='boot',B=1000,data=validation)

cal_hybrid_val_pro<-predict(object=fit_hybrid_train_lrm,type="fitted",newdata=validation)
fit_hybrid_val<-lrm(HE~cal_hybrid_val_pro,data=validation,x=T,y=T)
cal_hybrid_val_lrm<- calibrate(fit_hybrid_val, method='boot',B=1000,data=val)

plot(1,type="n",xlim=c(0,1),ylim=c(0,1),xlab = "Predicted probability",ylab = "Actual probability",legend=FALSE,subtitles=FALSE)
abline(0,1,col="black",lty=2,lwd=2)
lines(cal_R.score_val_lrm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="green",pch=16)
lines(cal_clinic_val_lrm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="blue",pch=16)
lines(cal_hybrid_val_lrm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="red",pch=16)
legend(0.55,0.18,c("Radiomics model","Clinical imaging sign model","Hybrid model","Ideal model"),lty = c(1,1,1,2),lwd=c(2,2,2,2),col=c("green","blue","red","black"),bty = "n")

#test 
cal_clinic_test_pro<-predict(object=fit_clinic_train_lrm,type="fitted",newdata=test)
fit_clinic_test<-lrm(HE~cal_clinic_test_pro,data=test,x=T,y=T)
cal_clinic_test_lrm<- calibrate(fit_clinic_test, method='boot',B=1000,data=test)

cal_R.score_test_pro<-predict(object=fit_r.score_train_lrm,type="fitted",newdata=test)
fit_R.score_test<-lrm(HE~cal_R.score_test_pro,data=test,x=T,y=T)
cal_R.score_test_lrm<- calibrate(fit_R.score_test, method='boot',B=1000,data=test)

cal_hybrid_test_pro<-predict(object=fit_hybrid_train_lrm,type="fitted",newdata=test)
fit_hybrid_test<-lrm(HE~cal_hybrid_test_pro,data=test,x=T,y=T)
cal_hybrid_test_lrm<- calibrate(fit_hybrid_test, method='boot',B=1000,data=test)

plot(1,type="n",xlim=c(0,1),ylim=c(0,1),xlab = "Predicted probability",ylab = "Actual probability",legend=FALSE,subtitles=FALSE)
abline(0,1,col="black",lty=2,lwd=2)
lines(cal_R.score_test_lrm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="green",pch=16)
lines(cal_clinic_test_lrm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="blue",pch=16)
lines(cal_hybrid_test_lrm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="red",pch=16)
legend(0.55,0.18,c("Radiomics model","Clinical imaging sign model","Hybrid model","Ideal model"),lty = c(1,1,1,2),lwd=c(2,2,2,2),col=c("green","blue","red","black"),bty = "n")


#############################################   DCA     ##########################################################
library(rmda)
#train
dca_train_r.score<-decision_curve(HE~Rad.score,data = train,family = binomial(link='logit'),thresholds=seq(0,1,by=0.01),confidence.intervals=0.95,study.design='case-control',population.prevalence=0.3)
dca_train_clinic<-decision_curve(HE~Time+Volume+Irregular + BlendSign  +SmokeHistory+ GCS,data = train,family = binomial(link='logit'),thresholds=seq(0,1,by=0.01),confidence.intervals=0.95,study.design='case-control',population.prevalence=0.3)
dca_train_hybrid<-decision_curve(HE ~ GCS+Volume+Irregular + Rad.score,data = train,family = binomial(link='logit'),thresholds=seq(0,1,by=0.01),confidence.intervals=0.95,study.design='case-control',population.prevalence=0.3)
list_dca_train<-list(dca_train_r.score,dca_train_clinic,dca_train_hybrid)
plot_decision_curve(list_dca_train,curve.names = c("Radiomics model","Clinical imaging sign model","Hybrid model"),xlab = "Threshold probability",ylab = "Net benefit",cost.benefit.axis = FALSE,col = c('green','blue','red'),confidence.intervals = FALSE,standardize = FALSE)
#validation
dca_val_r.score<-decision_curve(HE~predvalue_R.score_val,data = validation,family = binomial(link='logit'),thresholds=seq(0,1,by=0.01),confidence.intervals=0.95,study.design='case-control',population.prevalence=0.3)
dca_val_clinic<-decision_curve(HE~predvalue_clinic_validation,data = validation,family = binomial(link='logit'),thresholds=seq(0,1,by=0.01),confidence.intervals=0.95,study.design='case-control',population.prevalence=0.3)
dca_val_hybrid<-decision_curve(HE ~ predvalue_hybrid_val,data = validation,family = binomial(link='logit'),thresholds=seq(0,1,by=0.01),confidence.intervals=0.95,study.design='case-control',population.prevalence=0.3)
list_dca_val<-list(dca_val_r.score,dca_val_clinic,dca_val_hybrid)
plot_decision_curve(list_dca_val,curve.names = c("Radiomics model","Clinical imaging sign model","Hybrid model"),xlab = "Threshold probability",ylab = "Net benefit",cost.benefit.axis = FALSE,col = c('green','blue','red'),confidence.intervals = FALSE,standardize = FALSE)
#test
dca_test_r.score<-decision_curve(HE~predvalue_R.score_test,data = test,family = binomial(link='logit'),thresholds=seq(0,1,by=0.01),confidence.intervals=0.95,study.design='case-control',population.prevalence=0.3)
dca_test_clinic<-decision_curve(HE~predvalue_clinic_test,data = test,family = binomial(link='logit'),thresholds=seq(0,1,by=0.01),confidence.intervals=0.95,study.design='case-control',population.prevalence=0.3)
dca_test_hybrid<-decision_curve(HE ~ predvalue_hybrid_test,data = test,family = binomial(link='logit'),thresholds=seq(0,1,by=0.01),confidence.intervals=0.95,study.design='case-control',population.prevalence=0.3)
list_dca_test<-list(dca_test_r.score,dca_test_clinic,dca_test_hybrid)
plot_decision_curve(list_dca_test,curve.names = c("Radiomics model","Clinical imaging sign model","Hybrid model"),xlab = "Threshold probability",ylab = "Net benefit",cost.benefit.axis = FALSE,col = c('green','blue','red'),confidence.intervals = FALSE,standardize = FALSE)



############################################################################################
################IDI##################

library(PredictABEL)

#train
reclassification(data = train,cOutcome = 1,predrisk1 = train$predvalue_R.score_train, predrisk2 = train$predvalue_hybrid_train,cutoff = c(0,0.5,1))
reclassification(data = train,cOutcome = 1,predrisk1 = train$predvalue_clinic_train, predrisk2 = train$predvalue_hybrid_train,cutoff = c(0,0.5,1))
##validation
reclassification(data = validation,cOutcome = 1,predrisk1 = validation$predvalue_R.score_val, predrisk2 = validation$predvalue_hybrid_val,cutoff = c(0,0.5,1))
reclassification(data = validation,cOutcome = 1,predrisk1 = validation$predvalue_clinic_validation, predrisk2 = validation$predvalue_hybrid_val,cutoff = c(0,0.5,1))
##test
reclassification(data = test,cOutcome = 1,predrisk1 = test$predvalue_R.score_test, predrisk2 = test$predvalue_hybrid_test,cutoff = c(0,0.5,1))
reclassification(data = test,cOutcome = 1,predrisk1 = test$predvalue_clinic_test, predrisk2 = test$predvalue_hybrid_test,cutoff = c(0,0.5,1))


