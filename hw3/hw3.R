#HW3 Stat527
#Problem 2

n_s = c(20,50,100,200,500,1000)
mse = array(0, dim=c(10,8,length(n_s)))
for(i in 1:length(n_s)){
 #print(i)
   for(j in 1:10){
     for(d in 1:5){
       x = sort(runif(n_s[i]))
       y = sin(30*x) + rnorm(n_s[i]) #change to compare using other truth functions
     
         reg = lm(y ~ poly(x,degree=d))
         mse[j,d,i] = (1/n_s[i])*sum((fitted(reg) - y)^2)
      }
     
     reg = ksmooth(x,y,bandwidth = 0.05,kernel='normal')
     mse[j,6,i] = mean((y-reg$y)^2)
       
     reg = numeric(n_s[i])
     box_width = .05
     for(k in 1:n_s[i]){  #create predictions vector
         reg[k] = mean(y[abs(x-x[k])<box_width])
      }
     mse[j,7,i] = mean((y - reg)^2)
       
     reg = loess(y ~ x)
     mse[j,8,i] = mean((y - fitted(reg))^2)
   }
 }
average_mse = array(0, dim = c(8,length(n_s)))
for(i in 1:8){
   for(j in 1:length(n_s)){
       average_mse[i,j] = mean(mse[,i,j])
     }
 }
plot(n_s, average_mse[1,],xlim=c(0,1000),ylim=c(0.7,2),xlab='Sample Size', ylab = 'Average MSE',main='y=sin(30*x): Avg MSE vs Sample Size')
lines(n_s,average_mse[1,],lty=1)
for(i in 2:8){
   lines(n_s, average_mse[i,],col=i,lty=i)
 }
legend('topright',legend=c('deg1','deg2','deg3','deg4','deg5','NWBox','NWGauss','LocPoly'),col=1:8,lty=1:8,cex=0.7)
 
 