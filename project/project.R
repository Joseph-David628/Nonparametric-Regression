#Stat527 Project

library(kdensity)
library(logKDE)
library(rpart)

data = read_excel('/Users/joseph/Courses/stat527/project/data.xlsx')

#Only take data from 2021-2022
data = data[data$EVENT_DATE > '2021-01-01',]
data = data[data$COUNTRY=='Haiti',]
n = length(data$ISO) #number of observations

#plot histogram of fatalities
hist(sqrt(data$FATALITIES), freq = FALSE, ylim=c(0,2), xlim=c(0,7),
     main='Histogram and KDE of sqrt(Fatalities)', xlab = 'sqrt(Fatalities)')

kde = kdensity(sqrt(data$FATALITIES), kernel = "gaussian")
lines(kde, main = 'KDE of Fatalities', col = 'blue')

kde.civilians = kdensity(c(sqrt(data$FATALITIES[data$INTER1==7]), 
     sqrt(data$FATALITIES[data$INTER2==7])), kernel='gaussian')
plot(kde.civilians, main = 'KDE of Fatalities With/Without Civilians',
     col = 'red', xlab = 'sqrt(Fatalities)')
noncivil1 = which(data$INTER1!=7)
noncivil2 = which(data$INTER2!=7)
kde.noncivil = kdensity(sqrt(data$FATALITIES[intersect(noncivil1,noncivil2)]), 
                        kernel = "gaussian")
lines(kde.noncivil, col = 'blue')
legend('topright', legend=c('Civilians Present','Civilians Not Present'),
       col = c('red','blue'), lty=c(1,1))

#NW Regression on Lat/Long, lower bandwidth
bandwidth = 0.2
new_long <- seq(-77,-70, length = 30)
new_lat <- seq(16,22, length = 15)
reg = outer(new_lat, new_long, NW_estimator)
persp(new_long, new_lat, t(sqrt(reg)), main='NW Estimator with Gaussian Kernel', 
      col = 'white', theta = 0, phi = 25, shade = 0.75,lwd=0.1,
      xlab = 'Longitude', ylab = 'Latitude',zlab='Density Estimate',
      axes = TRUE, ticktype = 'detailed') #plot regression function
legend('topleft', legend=c('Bandwidth = .5'))

NW_estimator = function(lat,long){
  estimate = 0
  denom = 0
  for(i in 1:n){
    dist = sqrt((lat - data$LATITUDE[i])^2 + (long - data$LONGITUDE[i])^2)
    estimate = estimate + data$FATALITIES[i]*dnorm(dist, mean=0,sd = bandwidth)
    denom = denom + dnorm(dist, mean=0,sd = bandwidth)
  }  
  estimate = estimate / denom 
}

########################################################################
#Fit Tree Decision model on whether civilians present
data.tree = data
#remove when actors are inter 1 (only 7 observations)
civilians.inter1 = which(data$INTER1==7) 
#Set actor variable as type factor
data.tree$INTER1 = as.factor(data$INTER1)

fit.tree = rpart(CIVILIANS ~ LATITUDE + LONGITUDE + FATALITIES + INTER1, 
                    subset = setdiff(1:1251,civilians.inter1), method = 'class', 
                    data = data.tree, 
                    control=rpart.control(minsplit=20, minbucket=10, cp=0.01))
printcp(fit.tree)
plotcp(fit.tree)
summary(fit.tree)
plot(fit.tree, uniform=TRUE, 
     main="Classification Tree for Civilians Present") #plot tree
text(fit.tree, use.n=TRUE, all=TRUE, cex=.8)

#Classification Trees with Bagging
bs.sizes = seq(5,100, by=5)
bs.sizes = c(c(1,2,3,4), bs.sizes)
bs.train.errors = numeric(length(bs.sizes))
bs.test.errors = numeric(length(bs.sizes))

testData = data.tree[c("LATITUDE", "LONGITUDE", "FATALITIES", 
                       "INTER1", "CIVILIANS")]
trainData = data.tree[c("LATITUDE", "LONGITUDE", "FATALITIES", 
                        "INTER1", "CIVILIANS")]
indices = sample(1:n, 100)
testData = testData[indices,]
trainData = trainData[-indices,]
y.test.true = data.tree$CIVILIANS[indices]
y.train.true = data.tree$CIVILIANS[-indices]
n.train = length(trainData$LATITUDE)
n.test = length(testData$LATITUDE)

for(k in 1:length(bs.sizes)){
  train.predictions = numeric(n.train)
  test.predictions = numeric(n.test)
  for(i in 1:bs.sizes[k]){
    #Sample rows and fit a classification tree
    indices = sample(1:n.train, 100)
    fit.tree = rpart(CIVILIANS ~ LATITUDE + LONGITUDE + FATALITIES + INTER1, 
                     subset = indices, method = 'class', data = trainData, 
                     control=rpart.control(minsplit=5, minbucket=3, cp=0.1))
    #Generate predictions with this tree
    test.predictions = test.predictions + (predict(fit.tree, 
                                                testData, type='vector') - 1)
    train.predictions = train.predictions + (predict(fit.tree, 
                                                trainData, type='vector') - 1)
  }
  #classify by most common pred
  test.predictions = round(test.predictions / bs.sizes[k]) 
  train.predictions = round(train.predictions / bs.sizes[k])
  bs.test.errors[k] = sum(abs(test.predictions - y.test.true))/n.test
  bs.train.errors[k] = sum(abs(train.predictions - y.train.true))/n.train
}

plot(bs.sizes, bs.train.errors, xlab='Number of Bootstrap Samples',
     ylab='Classification Error', col = 'blue',ylim =c(.19, .3),
     main='Test vs Training Error')
lines(bs.sizes, bs.train.errors, col = 'blue')
lines(bs.sizes, bs.test.errors, col = 'green')
legend('topright', legend=c('Train','Test'), col = c('blue','green'),lty=1)

################################################################
#Boosting
data.boost = data.tree[c("LATITUDE", "LONGITUDE", "FATALITIES", 
                       "INTER1", "CIVILIANS")]
data.boost$CIVILIANS = 2*data.boost$CIVILIANS - 1
indices.xval = sample(1:n, 100)
indices.test = sample(setdiff(1:n, indices.xval), 100)
indices.train = setdiff(1:n, union(indices.test, indices.xval))

testData = data.boost[indices.test,]
xvalData = data.boost[indices.xval,]

#X-Validate to determine shrinkage param and no. of trees
shrinkage = seq(.001, .02, length = 10)
bag_sizes = seq(100, 2000, length = 20)

trainError = xvalError = array(0, dim = c(length(shrinkage), length(bag_sizes)))
n.train = length(indices.train)
n.xval = length(indices.xval)
n.test = length(indices.test)

y.train = data.boost$CIVILIANS[indices.train]
y.xval = data.boost$CIVILIANS[indices.xval]
y.test = data.boost$CIVILIANS[indices.test]

for(i in 1:length(shrinkage)){
  for(j in 1:length(bag_sizes)){
    pred.train = numeric(n.train)
    pred.xval = numeric(n.xval)
    trainData = data.boost[indices.train,]
    for(k in 1:bag_sizes[j]){
      fit.tree = rpart(CIVILIANS ~ LATITUDE + LONGITUDE + FATALITIES + INTER1, 
                        data = trainData, control=rpart.control(minsplit=5, 
                        minbucket=3, cp=0.07))
      
      new.pred = predict(fit.tree, trainData, type='vector')
      pred.train = pred.train + shrinkage[i]*new.pred
      new.xval = predict(fit.tree, xvalData, type='vector')
      pred.xval = pred.xval + shrinkage[i]*new.xval
      
      trainData$CIVILIANS = trainData$CIVILIANS - shrinkage[i]*new.pred
    }
    #pred.train = pred.train/abs(pred.train)
    #pred.xval = pred.xval/abs(pred.xval)
    
    trainError[i,j] = 0.5 * mean((pred.train - y.train)^2)
    xvalError[i,j] = 0.5 * mean((pred.xval - y.xval)^2)
  }
}

plot(bag_sizes, xvalError[1,], main = 'Cross-Validated Error', xlab='No. of Bags',
     ylab='Error', ylim = c(0.3, 0.5), type = 'l')
for(i in 1:length(shrinkage)){
  lines(bag_sizes, xvalError[i,], col=i)
}
legend('topright', legend = trunc(shrinkage*10^3)/10^3, col=1:10, lty=1)

plot(bag_sizes, trainError[1,], main = 'Training Error', 
     xlab='No. of Bags',
     ylab='Error', ylim = c(0.3, 0.6), type = 'l')
for(i in 1:length(shrinkage)){
  lines(bag_sizes, xvalError[i,], col=i)
}
legend('topright', legend = trunc(shrinkage*10^3)/10^3, col=1:10, lty=1)

#Fit tree with best hyperparameters
shrinkage = .012
num_bags = 750
pred.test = numeric(n.test)
indices.train = union(indices.train, indices.xval)
trainData = data.boost[indices.train,]
for(k in num_bags){
  fit.tree = rpart(CIVILIANS ~ LATITUDE + LONGITUDE + FATALITIES + INTER1, 
                   data = trainData, control=rpart.control(minsplit=5, 
                                        minbucket=3, cp=0.07))
  
  new.pred = predict(fit.tree, trainData, type='vector')
  new.pred.test = predict(fit.tree, testData, type='vector')
  pred.test = pred.test + new.pred.test*shrinkage
  trainData$CIVILIANS = trainData$CIVILIANS - shrinkage*new.pred
}
pred.test = pred.test/abs(pred.test)
error = 0.5 * mean(abs(pred.test - y.test))
plot(0.5*abs(pred.test - y.test))




