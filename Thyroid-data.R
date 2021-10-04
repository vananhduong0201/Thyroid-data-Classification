rm(list=ls())

library(dplyr)
library(ggplot2)
library(class)
library(caTools)
library(rpart)
library(rpart.plot)
library(scales)
library(caret)
library(e1071)
# Load the data set
Thyroid<-read.csv('Thyroid_data.csv',header=TRUE,sep=',')
str(Thyroid)
summary(Thyroid)

# check if there is any missing/null/NA/duplicated value in the data
any(is.na(Thyroid))
any(is.null(Thyroid))
duplicated(Thyroid)

# Checking for outliers
boxplot(Thyroid,main="Box Plot Analysis",xlab="", ylab="value")

#Basic EDA & plotting
  table(Thyroid$CLASS)

  par(mfrow=c(2,3))
  barplot(table(Thyroid$CLASS), main="Distribution of Class")
  hist(Thyroid$T3, main='T3 (Resin uptake test results)')
  hist(Thyroid$TST, main="TST")
  hist(Thyroid$TSTR, main="TSTR")
  hist(Thyroid$TSH, main="TSH")
  hist(Thyroid$MAD.TSH, main="MAD.TSH")
  
  plot(Thyroid, col=as.factor(Thyroid$CLASS))

  # Checking correlation
  corrplot(cor(Thyroid), "ellipse")
  
  #Normalized data using min-max method
  normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
  }
  Thyroid_2<-normalize(Thyroid)
  str(Thyroid_2)
  
  Thyroid_2<-Thyroid_2 %>%
    mutate(CLASS=Thyroid$CLASS)
  
  # Spliting data
  sample=sample.split(Thyroid_2$CLASS,SplitRatio=0.7)
  train=subset(Thyroid_2,sample==TRUE)
  test=subset(Thyroid_2,sample==FALSE)
  
  ThyroidVis<-Thyroid_2%>%mutate(train=sample)
  ggplot(ThyroidVis,aes(x=T3, y=TST, col=as.factor(CLASS)))+
    geom_point()+
    geom_jitter()+
    facet_wrap(~train)

  # Conducting classification model: decision trees
  
  tree_model<-rpart(CLASS~., data=train, method="class",minbucket=20)
  rpart.plot(tree_model)
  
  rpart.rules(tree_model)
  
  # Calcualte the training and test prediction
  treeTrainpred=predict(tree_model,train,type="class")
  treeTestpred=predict(tree_model,test,type="class")
  
  # Construct a confusion matrix (true class vs predicted class) 
  tabTrain=table(train$CLASS,treeTrainpred)  
  tabTest=table(test$CLASS,treeTestpred)  

  # Calculate for the performance measures
  
  accuracyTrain=(tabTrain[1,1]+tabTrain[2,2]+tabTrain[3,3])/sum(tabTrain)
  accuracyTest=(tabTest[1,1]+tabTest[2,2]+tabTest[3,3])/sum(tabTest)
  
  #Checking if accuracy rates are correct
  confusionMatrix(tabTrain)
  confusionMatrix(tabTest)
  
  # Construct a loop to see how the number of leaves
  # impacts the training and test performance
  treeTrain=matrix(0,nrow=50,ncol=200)
  treeTest=matrix(0,nrow=50,ncol=200)
  
  for(i in 1:50){
    for(j in 1:200){
      sample=sample.split(Thyroid_2$CLASS,SplitRatio=0.7)
      train=subset(Thyroid_2,sample==TRUE)
      test=subset(Thyroid_2,sample==FALSE)
      tree_model<-rpart(CLASS~., data=train, method='class',minbucket=i)
      treeTrainpred = predict(tree_model,train,type="class")
      treeTestpred = predict(tree_model,test,type="class")
      treeTrain[i,j]=1-mean(train$CLASS==treeTrainpred)
      treeTest[i,j]=1-mean(test$CLASS==treeTestpred)
    }
  }
  treeTrain=rowMeans(treeTrain)
  treeTest=rowMeans(treeTest)
  
  #Looking for appropriate leaves in the Tree model
  plot(c(1:50),treeTrain,type='l',main='Training(blue) vs Testing(red) error - Tree model',col='blue', xlab='Minimum leave size',ylab='Error')
  lines(treeTest, col="red")
  abline(v=20, col='green', lty=2)
  abline(v=25, col='green', lty=2)
  
  # Construct an overfit model 
  tree_model_2<-rpart(CLASS~., data=train, method='class', minbucket=1) # Overfit model
  rpart.plot(tree_model_2)
  
  # Calcualte the training and test accuracies
  treeTrain_2=predict(tree_model_2,train,type="class")
  treeTest_2=predict(tree_model_2,test,type="class")
  
  confusionMatrix(table(train$CLASS,treeTrain_2))
  confusionMatrix(table(test$CLASS,treeTest_2))
  
  
  ########################**************
  # KNN Model
  
  KnnTrain_pred=knn(train[,2:6],train[,2:6],as.factor(train$CLASS),12)
  table(train$CLASS,KnnTrain_pred)
  
  accuracyKTrain=sum(train$CLASS==KnnTrain_pred)/NROW(train$CLASS) 
  confusionMatrix(table(train$CLASS,KnnTrain_pred))
  
  KnnTest_pred=knn(train[,2:6],test[,2:6],as.factor(train$CLASS),12)
  table(test$CLASS,KnnTest_pred)
  
  accuracyKTest=sum(test$CLASS==KnnTest_pred)/NROW(test$CLASS) 
  confusionMatrix(table(test$CLASS,KnnTest_pred))
  
  ### Construct an underfit model
  KnnTrain_pred_Un=knn(train[,2:6],train[,2:6],as.factor(train$CLASS),2)
  KnnTest_pred_Un=knn(train[,2:6],test[,2:6],as.factor(train$CLASS),2)
  
  confusionMatrix(table(train$CLASS,KnnTrain_pred_Un))
  confusionMatrix(table(test$CLASS,KnnTest_pred_Un))
  
  # Construct a loop to see how the number of k
  # impacts the training and test performance
  KNNTrain=matrix(0,nrow=50,ncol=200)
  KNNTest=matrix(0,nrow=50,ncol=200)
  
  for(i in 1:50){
    for(j in 1:200){
      sample=sample.split(Thyroid_2$CLASS,SplitRatio=0.7)
      train=subset(Thyroid_2,sample==TRUE)
      test=subset(Thyroid_2,sample==FALSE)
      KnnTrain_pred=knn(train[,2:6],train[,2:6],as.factor(train$CLASS),i)
      KnnTest_pred=knn(train[,2:6],test[,2:6],as.factor(train$CLASS),i)
      
      KNNTrain[i,j]=1-mean(train$CLASS== KnnTrain_pred)
      KNNTest[i,j]=1-mean(test$CLASS== KnnTest_pred)
    }
  }
  KNNTrain=rowMeans(KNNTrain)
  KNNTest=rowMeans(KNNTest)
  ### Results: k from 10-20
  
  #Looking for appropriate k in the KNN model
  par(mfrow=c(1,1))
   plot(c(1:50),KNNTrain,type='l',main='Training(blue) vs Testing(red) error-KNN model',col='blue', xlab='Number of k',ylab='Error')
  lines(KNNTest, col="red")
  abline(v=9, col='green', lty=2)
  abline(v=15, col='green', lty=2)
 ###result: k from 9-15
 
####################***************
  #Redo the KNN without scaling
  sample_N<-sample.split(Thyroid$CLASS,SplitRatio=0.7)
  train_N<-subset(Thyroid,sample==TRUE)
  test_N<-subset(Thyroid,sample==FALSE)
  
  KnnTrain_pred_N = knn(train_N[,2:6],train_N[,2:6],as.factor(train_N$CLASS),10)
  KnnTest_pred_N = knn(train_N[,2:6],test_N[,2:6],as.factor(train_N$CLASS),10)
  
  confusionMatrix(table(train_N$CLASS,KnnTrain_pred_N)) 
  confusionMatrix(table(test_N$CLASS,KnnTest_pred_N)) 
  
  # Construct a loop to see how the number of k
  # impacts the training and test performance
  KNNTrain_N=matrix(0,nrow=50,ncol=200)
  KNNTest_N=matrix(0,nrow=50,ncol=200)
  
  for(i in 1:50){
    for(j in 1:200){
      sample_N=sample.split(Thyroid$CLASS,SplitRatio=0.7)
      train_N=subset(Thyroid,sample==TRUE)
      test_N=subset(Thyroid,sample==FALSE)
      KnnTrain_pred_N=knn(train_N[,2:6],train_N[,2:6],as.factor(train_N$CLASS),i)
      KnnTest_pred_N=knn(train_N[,2:6],test_N[,2:6],as.factor(train_N$CLASS),i)
      
      KNNTrain_N[i,j]=1-mean(train_N$CLASS== KnnTrain_pred_N)
      KNNTest_N[i,j]=1-mean(test_N$CLASS== KnnTest_pred_N)
    }
  }
  KNNTrain_N=rowMeans(KNNTrain_N)
  KNNTest_N=rowMeans(KNNTest_N)
  
  #Looking for appropriate k in the KNN model
  par(mfrow=c(1,1))
  plot(c(1:50),KNNTrain_N,type='l',main='Training(blue) vs Testing(red) error
       KNN model - no scaling',col='blue', xlab='Number of k',ylab='Error')
  lines(KNNTest_N, col="red")

  
  
  
  
  
  
  
