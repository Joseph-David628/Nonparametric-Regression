#Assignment 1

#Problem 2
n=200
x = sort(runif(n))
y1 = 2*x + rnorm(n)
y2 = sin(2*pi*x) + rnorm(n)
y3 = sin(30*x) + rnorm(n)

#polynomial with degrees 1-5
plot(x,y3)
mse1 = numeric(5)
for(d in 1:5){
  reg = lm(y3 ~ poly(x,degree=d))
  lines(x,fitted(reg),col=d)
  mse1[d] = (1/n)*sum((fitted(reg) - y3)^2)
}

#NW box kernel
y3_nw_box = numeric(n)
box_width = .04
for(i in 1:n){  #create predictions vector
  y3_nw_box[i] = mean(y3[abs(x-x[i])<box_width])
}
mse2 = mean((y3 - y3_nw_box)^2)
lines(x,y3_nw_box,col=6)

#NW gaussian kernel
y3_nw_gauss = ksmooth(x,y3,bandwidth = 0.04,kernel='normal')
lines(y3_nw_gauss,col=7)
mse3 = mean((y3 - y3_nw_gauss$y)^2)
legend('topleft',legend=c('deg 1', 'deg 2', 'deg 3', 'deg 4', 'deg 5', 'NW Box',
                          'NW Gauss'),col=1:7,lty=rep(1,7),cex=0.7)

#Problem 4
#Read data
fev = read.table('fev', header = TRUE, sep = "", dec = ".")
x = fev[,5]
y = fev[,4]
plot(x,y)
bandwidths = 1:10

#Box Kernel LOOCV
mse_box_loo = rep(0,length(bandwidths))
for(b in bandwidths){
  for(i in 1:length(x)){
    y_box_loo = ksmooth(x[-i],y[-i],bandwidth=b,kernel='box',x.points=c(x[i]))$y
    mse_box_loo[b] = mse_box_loo[b] + (y[i] - y_box_loo)^2
  }
}
mse_box_loo = mse_box_loo / length(x)

#Box Kernel 5-fold
mse_box_5 = rep(0,length(bandwidths))
for(b in bandwidths){
  indices = 1:length(x)
  for(i in 1:5){
    new_i = sample(indices, 130, replace=F)
    print(new_i)
    indices = setdiff(indices, new_i)
    y_box_5 = ksmooth(x[-new_i],y[-new_i],bandwidth=b,kernel='box',
                        x.points=x[new_i])$y
    mse_box_5[b] = mse_box_5[b] + mean((y[new_i] - y_box_5)^2)
  }
}
mse_box_5 = mse_box_5 / 5

#Epanechikov LOOCV
mse_ep_loo = rep(0,10)
for(b in bandwidths){
  for(i in 1:length(x)){
    y_ep_loo = sum((3/4)*(1-(abs(x[i]-x[-i])/(b/20))^2) * y[-i])/sum((3/4)*(1-(abs(x[i]-x[-i])/(b/20))^2))
    mse_ep_loo[b] = mse_ep_loo[b] + (y[i] - y_ep_loo)^2
  }
}
mse_ep_loo = mse_ep_loo / length(x)

#Epanechikov 5-fold
mse_ep_5 = rep(0,10)
for(b in bandwidths){
  indices = 1:654
  for(i in 1:5){
    new_i = sample(indices, 130, replace=F)
    indices = setdiff(indices, new_i)
    y_ep_5 = numeric(130)
    for(j in 1:130){
      y_ep_5[j] = sum((3/4)*(1-(abs(x[new_i[j]]-x[-new_i])/(b/10))^2) * y[-new_i])/sum((3/4)*(1-(abs(x[new_i[j]]-x[-new_i])/(b/10))^2))
    }
    mse_ep_5[b] = mse_ep_5[b] + mean((y[new_i] - y_ep_5)^2)
  }
}
mse_ep_5 = mse_ep_5 / 5


