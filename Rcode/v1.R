ls()
rm(list=ls())
generateA = function(n, a, k) {
  flag = 0
  
  while(flag == 0) {
    # first sample generate features
    nones = rpois(1, a)
    while(nones == 0) nones = rpois(1, a)
    A = matrix(rep(1, nones), nrow = 1)
    
    # subsequent samples
    for(i in 2 : n) {
      # update existing features
      A = rbind(A, runif(ncol(A)) < colSums(A) / i)
      # generate new features
      nones = rpois(1, a / i)
      if(nones != 0) {
        A = cbind(A, matrix(0, nrow(A), nones))
        A[i, (ncol(A) - nones) : ncol(A)] = 1
      }
    }
    
    # accepting criteria
    if(ncol(A) == k) flag = 1
  }
  return(A)
}

require(simFrame)


generatedata_overlap <- function(seednum){
  d <- c(15,200,15)
  set.seed(2022)
  y <- array(NA,dim = d)
  r <- 3
  m <- 1
  rho <- 0.3
  nac<-NAControl(NArate=0.5)
  gamma <- c(rho/2,1-rho,rho/2)
  mu <- 0
  sigma2 <- 1
  v <- 7
  c1o <- generateA(d[1],m,r)
  #c2o <- matrix(rbinom(d[2]*r,1,rho),nrow = d[2])
  c2o = matrix(as.numeric(runif(d[2] * r) < rho), d[2], r)
  #c3o <- matrix(sample(c(-1,0,1),d[3]*r,TRUE,gamma),nrow = d[3])
  c3o = matrix(as.numeric(runif(d[3] * r) < rho), d[3], r)
  c3o[c3o == 1] = 2 * (runif(sum(c3o == 1)) < 0.5) - 1
  l1 <- seq(6.0,6+0.5*(r-1),length.out = r)
  #l2 <- t(seq(2.0,4.5,by = 0.5)%o%rep(1,d[3]))
  l2 <- l1
  
  lambda1 <- matrix(rep(l1,d[3]),nrow = d[3],byrow = TRUE)
  lambda2 <- lambda1
  #generate b1,b2
  #b1 <- rnorm(d[2]*d[3])
  b1 <- log(0.1)
  #b2 <- rnorm(d[2]*d[3])
  b2 <- b1
  b1m <- matrix(b1,nrow = d[2],ncol = d[3])
  b2m <- b1m
  
  #generate z & y
  theta1 <- array(0,dim =d)
  theta2 <- array(0,dim =d)
  for(i in 1:r){
    theta1 <- theta1+l1[i]*c1o[,i]%o%c2o[,i]%o%(c3o==-1)[,i]
    theta2 <- theta2+l2[i]*c1o[,i]%o%c2o[,i]%o%(c3o==1)[,i]
  }
  theta1 <- theta1+b1
  theta2 <- theta2+b2
  
  theta1 <- exp(theta1)
  theta2 <- exp(theta2)
  set.seed(seednum)
  # generate latent ternary indicator z
  U0 = array(runif(d[1]*d[2]*d[3]),dim = d)
  
  z0 = (U0 > (theta1 + 1) / (theta1 + theta2 + 1)) - (U0 < theta1 / (theta1 + theta2 + 1))
  
  #y[z0==1] <- rnorm(sum(z0==1),mu+10,5)
  y[z0==1] <- runif(sum(z0==1),mu,mu+5)
  
  y[z0==0] <- rnorm(sum(z0==0),mu,sigma2)
  
  #y[z0==-1] <- rnorm(sum(z0==-1),mu-10,5)
  y[z0==-1] <- runif(sum(z0==-1),mu-5,mu)
  yd <- data.frame(as.vector(y))
  yna <- setNA(yd,nac)
  yna <- yna[,1]
  yna <- array(yna,dim = d)
  return(list(y,yna,c1o,c2o,c3o))
}

generate_nonoverlap <- function(seednum){
  set.seed(2)
  d <- c(15,200,15)
  y <- array(NA,dim = d)
  r <- 3
  rho <- 0.3
  gamma <- c(rho/2,1-rho,rho/2)
  
  mu <- 0
  sigma2 <- 1
  v <- 7
  
  # generate feature <-> latent feature matrix C (ternary): number of features p = 50
  #pi <- rdirichlet(1,c(2,2,6))
  c1o <- t(rmultinom(d[1],1,c(3,3,4)))
  #c2o <- matrix(rbinom(d[2]*r,1,rho),nrow = d[2])
  c2o <- t(rmultinom(d[2],1,c(3,3,4)))
  #c3o <- matrix(sample(c(-1,0,1),d[3]*r,TRUE,gamma),nrow = d[3])
  c3o <- t(rmultinom(d[3],1,c(3,3,4)))
  c3o[c3o == 1] = 2 * (runif(sum(c3o == 1)) < 0.5) - 1
  l1 <- seq(6.0,6+0.5*(r-1),length.out = r)
  #l2 <- t(seq(2.0,4.5,by = 0.5)%o%rep(1,d[3]))
  l2 <- l1
  
  lambda1 <- matrix(rep(l1,d[3]),nrow = d[3],byrow = TRUE)
  lambda2 <- lambda1
  #generate b1,b2
  #b1 <- rnorm(d[2]*d[3])
  b1 <- log(0.1)
  #b2 <- rnorm(d[2]*d[3])
  b2 <- b1
  b1m <- matrix(b1,nrow = d[2],ncol = d[3])
  b2m <- b1m
  
  #generate z & y
  theta1 <- array(0,dim =d)
  theta2 <- array(0,dim =d)
  for(i in 1:r){
    theta1 <- theta1+l1[i]*c1o[,i]%o%c2o[,i]%o%(c3o==-1)[,i]
    theta2 <- theta2+l2[i]*c1o[,i]%o%c2o[,i]%o%(c3o==1)[,i]
  }
  theta1 <- theta1+b1
  theta2 <- theta2+b2
  
  theta1 <- exp(theta1)
  theta2 <- exp(theta2)
  set.seed(seednum)
  # generate latent ternary indicator z
  U0 = array(runif(d[1]*d[2]*d[3]),dim = d)
  
  z0 = (U0 > (theta1 + 1) / (theta1 + theta2 + 1)) - (U0 < theta1 / (theta1 + theta2 + 1))
  
  #y[z0==1] <- rnorm(sum(z0==1),mu+10,5)
  y[z0==1] <- runif(sum(z0==1),mu,mu+5)
  
  y[z0==0] <- rnorm(sum(z0==0),mu,sigma2)
  
  #y[z0==-1] <- rnorm(sum(z0==-1),mu-10,5)
  y[z0==-1] <- runif(sum(z0==-1),mu-5,mu)
  return(list(y,c1o,c2o,c3o))
}

generatedata_na <- function(seednum,nac){
  d <- c(15,200,15)
  set.seed(2022)
  y <- array(NA,dim = d)
  r <- 3
  m <- 1
  rho <- 0.3
  gamma <- c(rho/2,1-rho,rho/2)
  mu <- 0
  sigma2 <- 1
  v <- 7
  c1o <- generateA(d[1],m,r)
  #c2o <- matrix(rbinom(d[2]*r,1,rho),nrow = d[2])
  c2o = matrix(as.numeric(runif(d[2] * r) < rho), d[2], r)
  #c3o <- matrix(sample(c(-1,0,1),d[3]*r,TRUE,gamma),nrow = d[3])
  c3o = matrix(as.numeric(runif(d[3] * r) < rho), d[3], r)
  c3o[c3o == 1] = 2 * (runif(sum(c3o == 1)) < 0.5) - 1
  l1 <- seq(6.0,6+0.5*(r-1),length.out = r)
  #l2 <- t(seq(2.0,4.5,by = 0.5)%o%rep(1,d[3]))
  l2 <- l1
  
  lambda1 <- matrix(rep(l1,d[3]),nrow = d[3],byrow = TRUE)
  lambda2 <- lambda1
  #generate b1,b2
  #b1 <- rnorm(d[2]*d[3])
  b1 <- log(0.1)
  #b2 <- rnorm(d[2]*d[3])
  b2 <- b1
  b1m <- matrix(b1,nrow = d[2],ncol = d[3])
  b2m <- b1m
  
  #generate z & y
  theta1 <- array(0,dim =d)
  theta2 <- array(0,dim =d)
  for(i in 1:r){
    theta1 <- theta1+l1[i]*c1o[,i]%o%c2o[,i]%o%(c3o==-1)[,i]
    theta2 <- theta2+l2[i]*c1o[,i]%o%c2o[,i]%o%(c3o==1)[,i]
  }
  theta1 <- theta1+b1
  theta2 <- theta2+b2
  
  theta1 <- exp(theta1)
  theta2 <- exp(theta2)
  set.seed(seednum)
  # generate latent ternary indicator z
  U0 = array(runif(d[1]*d[2]*d[3]),dim = d)
  
  z0 = (U0 > (theta1 + 1) / (theta1 + theta2 + 1)) - (U0 < theta1 / (theta1 + theta2 + 1))
  
  #y[z0==1] <- rnorm(sum(z0==1),mu+10,5)
  y[z0==1] <- runif(sum(z0==1),mu,mu+5)
  
  y[z0==0] <- rnorm(sum(z0==0),mu,sigma2)
  
  #y[z0==-1] <- rnorm(sum(z0==-1),mu-10,5)
  y[z0==-1] <- runif(sum(z0==-1),mu-5,mu)
  yna <- y
  #nac <- sample(1:(d[1]*d[2]),0.5*d[1]*d[2])
  for (g in 1:d[3]) {
    #for (t in 1:d[1]) {
    yna[,,g][nac] <- NA
    #}
  }
  return(list(y,yna,c1o,c2o,c3o))
}

require(e1071)
min_ham = function(A, B) {
  ## list all permutations
  permutation = e1071::permutations(ncol(A))
  
  DHamming = rep(NA, nrow(permutation))
  for(l in 1 : nrow(permutation)) {
    NB = as.vector(B[, permutation[l, ]])
    DHamming[l] = sum(abs(as.vector(A) - NB))
  }
  return(min(DHamming))
}

modev <- function(x){
  return(as.numeric(names(table(x))[table(x) == max(table(x))]))}
tomembershipmatrix <- function(z){
  d = length(z)
  Z = list()
  for (i in 1:d){
    zi = z[[i]]
    pi = length(zi)
    ri = length(unique(zi))
    Zi = matrix(0, pi, ri)
    for (k in 1:ri){
      zi.subset = which(zi==k)
      Zi[zi.subset,k] = 1
    }
    Z = c(Z, list(Zi))
  }
  return(Z)
}
perE = function(A, B) {
  ## list all permutations
  permutation = e1071::permutations(ncol(A))
  
  DHamming = rep(NA, nrow(permutation))
  for(l in 1 : nrow(permutation)) {
    NB = as.vector(B[, permutation[l, ]])
    DHamming[l] = sum(abs(as.vector(A) - NB))
  }
  
  return(permutation[which.min(DHamming), ])
}
require("doParallel")   
require("foreach") 
cl<- makeCluster(8) 
registerDoParallel(cl)
bayes_multi_n <- function(y,c1,c2,c3){
  set.seed(1)
  {
    c1track <- list()
    c2track <- list()
    c3track <- list()
    al=1;bl=3;av=0.5;bv=5;am=5;bm=1;sigmal=1;mub=-1
    ;sigmab=5;arho=1;brho=1;psi_1=1;psi0=8;psi1=1;as=1;bs=1;sigmamu=10;
    
    #//dimension
    d=dim(y);
    
    #//class number
    r=1;
    c1 = matrix(as.numeric(runif(d[1] * r) < 0.5), d[1], r)
    c2 = matrix(as.numeric(runif(d[2] * r) < 0.5), d[2], r)
    c3 = matrix(as.numeric(runif(d[3] * r) < 0.5), d[3], r)
    lambda1 <- matrix(rgamma(d[3]*r, al, 1/bl), d[3],r)
    lambda2 <- lambda1
    #//paramters
    mu1 <- rep(0,d[1])
    mu2 <- rep(0,d[2])
    mu3 <- rep(0,d[3])
    sigma2 <- rep(1,d[1])
    v1 <- 5*sigma2
    v2 <- v1
    #b1 <- matrix(rnorm(d[2]*d[3]),nrow = d[2],ncol = d[3])
    b1 <- matrix(-2,nrow = d[2],ncol = d[3])
    b2 <-b1
    m=3;rho=0.3;g=1;
    
    #vec gamma(3,fill::randu);
    gamma <- c(rho/2,1-rho,rho/2)
    #z <- z0
    z <- array(0,dim = d)
    z[which(is.na(y))] <- NA
  }
  #//MCMC update
  pb <- progress_bar$new(format = "  complete [:bar] :percent eta: :eta",
                         total = 5000, clear = FALSE, width= 60)
  for(it in 1:5000){
    #ztrack[[it]] <- z
    
    
    z <- update_z(c1, c2, c3, lambda1, lambda2, d, z, b1, b2, v1, v2, av, bv, y, mu1, mu2, mu3, sigma2);
    v1 <- update_v1(z,y,mu1,mu2,mu3,v1,v2,d,av,bv,sigma2)
    v2 <- update_v2(z,y,mu1,mu2,mu3,v1,v2,d,av,bv,sigma2)
    li1 <- update_c1(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,gamma, m,g,rho,al,bl)
    c1 <- li1[[1]]
    c2 <- li1[[2]]
    c3 <- li1[[3]]
    lambda1 <- li1[[4]]
    lambda2 <- li1[[5]]
    r <- ncol(c1)
    c2 <- update_c2(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,rho);
    li2 <-  update_c3(c1,c2,c3,lambda1,lambda2,d,z,b1,b2,gamma,sigmal,al,bl)
    c3 <-li2[[1]]
    lambda1 <- li2[[2]]
    lambda2 <- li2[[3]]
    
    #m <- rgamma(1,am+r,bm+har(d[1]))
    #//update $\rho$
    rho <-  rbeta(1,arho+sum(c2),brho+d[2]*r-sum(c2))
    #//update $\gamma$
    probs <- c(psi_1+sum(c3==-1),psi0+sum(c3==0),psi1+sum(c3==1))
    gamma <- rdirichlet(1,probs)
    li5 <- update_bunit(d,mub,sigmab,b1,b2,c1,c2,c3,lambda1,lambda2,z)
    b1 <- li5[[1]]
    b2 <- li5[[2]]
    sigma2 <- update_sigma2(z,y,mu1,mu2,mu3,v1,v2)
    mu1 <- update_mu1(z,y,sigmamu,sigma2,mu1,mu2,mu3)
    #mu2 <- update_mu2(z,y,sigmamu,sigma2,mu1,mu2,mu3)
    mu3 <- update_mu3(z,y,sigmamu,sigma2,mu1,mu2,mu3)
    c1track[[it]] <- c1
    c2track[[it]] <- c2
    c3track[[it]] <- c3
    pb$tick()
    Sys.sleep(1/100)
    print(it)
  }
  c1track <- c1track[2501:5000]
  c2track <- c2track[2501:5000]
  c3track <- c3track[2501:5000]
  rvec <- unlist(lapply(c1track, function(x) return(ncol(x))))
  rp <- as.numeric(names(table(rvec))[table(rvec) == max(table(rvec))])
  c1track <- c1track[which(rvec==rp)]
  c2track <- c2track[which(rvec==rp)]
  c3track <- c3track[which(rvec==rp)]
  c1p <- matrix(unlist(c1track),nrow = d[1]*rp)
  c2p <- matrix(unlist(c2track),nrow = d[2]*rp)
  c3p <- matrix(unlist(c3track),nrow = d[3]*rp)
  c1p <- matrix(rowMeans(c1p),nrow=d[1])
  c2p <- matrix(rowMeans(c2p),nrow=d[2])
  c3p <- matrix(rowMeans(c3p),nrow=d[3])
  return(list(c1p,c2p,c3p))
}
ob2= foreach(i=1:20,.packages=c("simFrame","progress","testarma","gtools")) %dopar% 
  {
    data <- generatedata_overlap(i)
    y <- data[[1]]
    #yna <- data[[2]]
    c1 <- data[[3]]
    c2 <- data[[4]]
    c3 <- data[[5]]
    bayes_multi_n(y,c1,c2,c3)
  }
save.image("sen(av=0.5,bv=5).RData")