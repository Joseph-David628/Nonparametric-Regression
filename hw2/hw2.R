#Stat527 HW2

#Problem 3
local_poly <- function(x, y, deg, bandwidth, kernel = 1, new_x = x){
  n = length(new_x)
  f_hat = numeric(n)
  for(k in 1:n){
    x0 = x[k]
    beta_hat = get_beta(x0, x, y, deg, bandwidth, kernel)
    zi = numeric(deg)
    for(j in 1:(deg+1)){zi[j] = x0^(j-1)}
    f_hat[k] = zi %*% beta_hat
  }
  return(f_hat)
}

uniform <- function(x, bandwidth){
  if(abs(x) > bandwidth){
    return(0)
  }
  else{return(1)}
}

get_beta = function(x0,x,y,deg, bandwidth, kernel = 1){
  n = length(x)
  Z = matrix(0,n,deg+1)
  W = matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:(deg+1)){
      Z[i,j] = x[i]^(j-1)
    }
  }
  if(kernel == 1){
  for(i in 1:n){
    W[i,i] = uniform(x[i] - x0, bandwidth)
  }
  }
  if(kernel == 2){
    for(i in 1:n){
      W[i,i] = dnorm(x[i] - x0, mean = 0, sd=bandwidth)
    }
  }
  if(kernel == 3){
    for(i in 1:n){
      if(abs(x0-x[i]) > bandwidth){W[i,i]=0}
      else{W[i,i] = (3/4)*(1-(abs(x[i]-x0)/bandwidth)^2)}
    }
  }
  return(solve(t(Z)%*%W%*%Z)%*%t(Z)%*%W%*%y)
}