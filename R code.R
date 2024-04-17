library(readr)
library(irr)
library(lpSolve)
library(reportROC)
library(pROC)

####################################   ICC      ###############################


feature_1 <- read.csv("./ICC1.csv")
feature_3 <- read.csv("./ICC2.csv")

data_name=names(feature_1)[2:ncol(feature_1)]
icc_val<-vector(length=1219)
selected <- feature_1[feature_1$ID %in% feature_3$ID,]
for (i in data_name){
  ratings <- cbind(selected[,i],feature_3[,i])
  icc <- icc(ratings, model = "twoway", 
             type = "agreement", 
             unit = "single", r0 = 0, conf.level = 0.95)
  icc_val[i] <- icc$value
}
icc_val["wavelet.LLL_glcm_JointEntropy"]
thr <- 0.75
Index <- as.data.frame(which(icc_val > thr))
Index_name <- row.names(Index)
Index_ICC <- which(icc_val > thr)

ICC_after_radiomics_x<- data_Radiomics[,names(Index_ICC )]
ICC_after_radiomics <- cbind(data_Radiomics$ID,cbind(data_Radiomics$FZ,cbind(data_Radiomics$HE,ICC_after_radiomics_x)))
names(ICC_after_radiomics)[1] <- "ID"
names(ICC_after_radiomics)[2] <- "FZ"
names(ICC_after_radiomics)[3] <- "HE"


#####  Standardization  #######
DataTrain_Radiomics <- data_train[,names(ICC_after_radiomics )]
DataTest_Radiomics <- data_test[,names(ICC_after_radiomics )]
DataVal_Radiomics <- data_gy[,names(ICC_after_radiomics )]

data_train_reduce <- DataTrain_Radiomics[,-(1:3)]
data_test_reduce <- DataTest_Radiomics[,-(1:3)]
data_val_reduce <- DataVal_Radiomics[,-(1:3)]


library(caret)
normal_para <- preProcess(x=data_train_reduce,mothod=c("center","scale"))
df_train_normal <- predict(object = normal_para,newdata = data_train_reduce)
df_test_normal <- predict(object = normal_para,newdata = data_test_reduce)
df_val_normal <- predict(object = normal_para,newdata = data_val_reduce)

data_train_stantard<- cbind(data_train$ID,cbind(data_train$FZ,cbind(data_train$HE,df_train_normal)))
names(data_train_stantard)[1] <- "ID"
names(data_train_stantard)[2] <- "FZ"
names(data_train_stantard)[3] <- "HE"

data_test_stantard<- cbind(data_test$ID,cbind(data_test$FZ,cbind(data_test$HE ,df_test_normal)))
names(data_test_stantard)[1] <- "ID"
names(data_test_stantard)[2] <- "FZ"
names(data_test_stantard)[3] <- "HE"

data_val_stantard<- cbind(data_gy$ID,cbind(data_gy$FZ,cbind(data_gy$HE ,df_val_normal)))
names(data_val_stantard)[1] <- "ID"
names(data_val_stantard)[2] <- "FZ"
names(data_val_stantard)[3] <- "HE"



################
##### T  test#######
################


data_t1=names(data_train_stantard)[4:ncol(data_train_stantard)]
pval=c()

for(i in data_t1){
  p=t.test(data_train_stantard[,i]~data_train_stantard$HE,var.equal=T)  
  pval[i]=p$p.value  
}
pval

thr <- 0.05
Index_t1 <- which(pval < thr)



names(Index_t1)
data_train_t_after_reduce <- data_train_stantard[,names(Index_t1)]
data_train_t_after_reduce
data_test_t_after_reduce <- data_test_stantard[,names(Index_t1)]
train_t_after<- cbind(data_train$ID,cbind(data_train$FZ,cbind(data_train$HE,data_train_t_after_reduce)))

names(train_t_after)[1] <- "ID"
names(train_t_after)[2] <- "FZ"
names(train_t_after)[3] <- "HE"

test_t_after<- cbind(data_test$ID,cbind(data_test$FZ,cbind(data_test$HE,data_test_t_after_reduce)))
names(test_t_after)[1] <- "ID"
names(test_t_after)[2] <- "FZ"
names(test_t_after)[3] <- "HE"

##################################
#############LASSO###############
##################################
library(glmnet)
library(rms)
library(foreign)
library(car)
library(ggplot2)
library(MASS)
library(survival)
library(survminer)
library(ggridges)
library(pROC)
library(plotROC)
library(riskRegression)
library(epiDisplay)
library(dplyr)
library(caret)
dev01 <- train_t_after
dev <- dev01[,-(1:2)]
ddist<-datadist(dev)
options(datadist="ddist")
y<-as.matrix(dev[,1])
x<-as.matrix(dev[,2:ncol(dev)])

fit3<-glmnet(x,y,alpha=1,family='binomial')
plot(fit3, xvar = "lambda", label = TRUE)

library(ggplot2)
set.seed(3699)
cv.fit3 <- cv.glmnet(x,y,family='binomial',type.measure = "deviance",nfolds =10)
plot(cv.fit3)
plot(fit3, xvar = "lambda", label = TRUE)
abline(v=log(c(cv.fit3$lambda.min,cv.fit3$lambda.1se)),lty=2)

cv.fit3$lambda.1se
log(cv.fit3$lambda.1se)
Coefficients <- coef(fit3, s = cv.fit3$lambda.1se)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Index
Active.Coefficients
row.names(Coefficients)[Active.Index]
lasso_name <- row.names(Coefficients)[Active.Index][-1]
lasso_name

data_train_lasso <- data.frame(train_t_after)[lasso_name]
data_test_lasso <- data_test_stantard[,names(data_train_lasso)]

train_lasso_after<- cbind(data_train$ID,cbind(data_train$FZ,cbind(data_train$HE,data_train_lasso)))
names(train_lasso_after)[1] <- "ID"
names(train_lasso_after)[2] <- "FZ"
names(train_lasso_after)[3] <- "HE"

test_lasso_after<- cbind(data_test$ID,cbind(data_test$FZ,cbind(data_test$HE,data_test_lasso)))
names(test_lasso_after)[1] <- "ID"
names(test_lasso_after)[2] <- "FZ"
names(test_lasso_after)[3] <- "HE"

##########################
########## pearson #######
##########################
norm_result <- apply(data_train_lasso, 2, function(x) shapiro.test(x)$p.value)
norm_feature <- data_train_lasso[which(norm_result >= 0.05)]

cor_nor <- cor(norm_feature, method = "pearson")
cor_all <- cor(data_train_lasso, method = "spearman")

num_nor <- dim(cor_nor)[1]

cor_all[upper.tri(cor_all)] <- 0
diag(cor_all) <- 0
data_reduce = data_train_lasso[, !apply(cor_all, 2, function(x) any(abs(x) > 0.8))]
dim(data_reduce)

data_train_model<- data_train_lasso[,names(data_reduce)]
data_test_model <- data_test_lasso[,names(data_reduce)]
train_model<- cbind(data_train$ID,cbind(data_train$FZ,cbind(data_train$HE,data_train_model)))
names(train_model)[1] <- "ID"
names(train_model)[2] <- "FZ"
names(train_model)[3] <- "HE"
test_model<- cbind(data_test$ID,cbind(data_test$FZ,cbind(data_test$HE,data_test_model)))
names(test_model)[1] <- "ID"
names(test_model)[2] <- "FZ"
names(test_model)[3] <- "HE"
val_model<- data_val_stantard[,names(train_model)]

write.csv(train_model,"train_model.csv",row.names = F)
write.csv(test_model,"test_model.csv",row.names = F)
write.csv(val_model,"val_model.csv",row.names = F)

