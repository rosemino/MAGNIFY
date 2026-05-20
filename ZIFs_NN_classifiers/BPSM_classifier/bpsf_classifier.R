# ANALISIS COLECTIVE VARIABLES PLUMED

# IMPORT LIBRARIES
library(readr)
library(dplyr)
library(ggplot2)
library(forcats)
library(stringr)
library(tidymodels)
library(discrim)
library(LiblineaR)
library(neuralnet)
library(keras)
library(RColorBrewer)

################################################################
# LOAD DATA

plumed = read.table(file='fullset.dat' ,sep='')
#######################################################

# NORMALIZE DATA
# I WANT EACH VARIABLE TO BE BETWEEN 0 AND 1

xmin  = 1:(dim(plumed)[2]-1)
xmax  = 1:(dim(plumed)[2]-1)
delta = 1:(dim(plumed)[2]-1)

for (i in {1:(dim(plumed)[2]-1)})
{
  xmin[i] = min(plumed[i])
  xmax[i] = max(plumed[i])
  delta[i] = xmax[i]-xmin[i]
  plumed[i] = (plumed[i]-xmin[i])/delta[i]
  }

#SPLIT DATA INTO TRAIN AND TEST SETS

set.seed(421)
split <- initial_split(plumed, prop = 0.8, strata = V13)
train = training(split)
test = testing(split)

###############################################################
###############################################################
###############################################################

# NEURAL NETWORK (NON LINEAR CLASSIFIER)
## NN TENSORFLOW

# DEFINE NET
model <- keras_model_sequential() %>%
  layer_dense(units = 24, activation = "relu", input_shape = c(dim(train)[2]-1)) %>%
  layer_dense(units = 7, activation = "softmax")

#DEFINE ALGORITHM
compile(model,
  loss='categorical_crossentropy',
  optimizer = optimizer_adam(),
  metrics = c('accuracy')     
)

# DATA PREPARATION

x_train = as.matrix(train[1:(dim(train)[2]-1)])
y_train = as.matrix(to_categorical(train$V13))
x_test = as.matrix(test[1:(dim(test)[2]-1)])
y_test = as.matrix(to_categorical(test$V13))

# FITTING
model %>% fit(x_train, y_train, epochs = 500, batch_size = 128)

#SAVE MODEL 
save_model_weights_tf(model, '7states_2.0') 

scaling = t(rbind(xmin,xmax))
write.table(scaling,file='scaling_7states_2.0.dat')

#COMPUTE TEST SET
score <- model %>% evaluate(x_test, y_test)
predict  = predict(model, x_test)

#confusion matrix
predict_factor = max.col(predict)-1
y_test_factor = max.col(y_test)-1
table(predict_factor,y_test_factor)














