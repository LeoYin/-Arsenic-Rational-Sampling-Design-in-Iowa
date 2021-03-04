library(Matrix)
library(MASS)
library(igraph)
library(deldir)
library(cccd)

Get_KNN<-function(coords,K,r){
  n=nrow(coords)
  g=nng(coords,k=K)
  g.MST.dataframe <- igraph::as_data_frame(g)
  from <- g.MST.dataframe$from
  to <- g.MST.dataframe$to
  for(i in 1:(n-1)){
    print(i)
    dst=sqrt((coords[,1]-coords[i,1])^2+(coords[,2]-coords[i,2])^2)
    candi=which(dst<=r)
    candi=candi[which(candi>i)]
    tmp=c(from[which(to==i)],to[which(from==i)])
    for(k in candi){
      if(!(k %in% tmp)){
        from=c(from,i)
        to=c(to,k)
      }
    }
  }
  g <- graph_from_edgelist(cbind(from,to))
  return(g)
}

Get_H <- function(g){
  np <- length(V(g))
  g.MST.dataframe <- igraph::as_data_frame(g)
  from <- g.MST.dataframe$from
  to <- g.MST.dataframe$to

  nv <- length(from)
  idi <- c(1:nv,1:nv)
  idj <- c(from,to)
  x <- c(rep(1,nv),rep(-1,nv))
  H <- sparseMatrix(i=idi, j=idj, x=x, dims=c(nv,np))
  
  return(H)
}


lasso_shrinkage<-function(a, kappa, w){
  
  sig=sign(a)
  y=abs(a)-kappa*w
  y[which(y<0)]=0
  y=sig*y
  return(y);
}

admm.Bern<-function(y, H, lambda.list, lr=0.01,
                       reltol=1e-04, abstol=1e-02/length(y), maxiter=1000, rho=1,
                       w=NULL, init=NULL, beta0=NULL){
  # 1. get parameters
  n=length(y); m=dim(H)[1]; B=m/n;
  if(is.null(w)){
    w=1
  }
  
  # 2. set ready
  if(is.null(init)){
    beta=runif(n)
  }else{
    beta=init
  }
  
  # 3. precompute static variables for x-update and factorization
  HtH = t(H)%*%H
  
  # 4. iteration
  sqrtn = sqrt(n);
  lambda.list.sort=sort(lambda.list,decreasing=TRUE)
  nlambda=length(lambda.list)
  output=list()
  output$lambda=lambda.list.sort
  output$x=matrix(0,n,nlambda)
  output$df=rep(0,nlambda)
  output$k=rep(0,nlambda)
  output$objval=rep(0,nlambda)
  output$r_norm=rep(0,nlambda)
  output$s_norm=rep(0,nlambda)
  output$eps_pri=rep(0,nlambda)
  output$eps_dual=rep(0,nlambda)
  
  for(j in 1:length(lambda.list.sort)){
    lambda=lambda.list.sort[j]
    beta=runif(n)
    z=H%*%beta
    u=rep(0,m)
    Hbeta = z
    
    h_objval=rep(0,maxiter)
    h_r_norm=rep(0,maxiter)
    h_s_norm=rep(0,maxiter)
    h_eps_pri=rep(0,maxiter)
    h_eps_dual=rep(0,maxiter)
    
    for(k in 1:maxiter){
      # 4-1. update 'x'
      b = exp(-beta)
      W = b/(1+b)^2
      D = Diagonal(x=W)+rho*HtH
      L = Cholesky(D, super=TRUE)
      q = -(y*b/(1+b)-(1-y)/(1+b))+rho*t(H)%*%(Hbeta-z+u)
      beta = beta-as.vector(solve(L,q))*lr
      
      # 4-2. update 'z'
      Hbeta = H%*%beta
      z_old = z;
      z = lasso_shrinkage(Hbeta + u, lambda/rho/B, w);
      
      # 4-3. update 'u'
      u = u + Hbeta - z;
      
      # 4-4. dianostics, reporting
      beta_norm = mean(abs(beta))
      z_norm = mean(abs(-z))
      # h_objval[k] = lasso_objective(X,y,v,Hbeta,lambda,beta);
      h_r_norm[k] = mean(abs(Hbeta-z));
      h_s_norm[k] = mean(abs(-rho*(z-z_old)));
      if (beta_norm>z_norm){
        h_eps_pri[k] = sqrtn*abstol + reltol*beta_norm;
      } else {
        h_eps_pri[k] = sqrtn*abstol + reltol*z_norm;
      }
      h_eps_dual[k] = sqrtn*abstol + reltol*mean(abs(rho*u));
      
      # 4-5. termination
      if ((h_r_norm[k] < h_eps_pri[k])&(h_s_norm[k]<h_eps_dual[k])){
        break;
      }
      
      print(k)
      if(is.null(beta0)==FALSE){
        print(mean((1/(1+exp(-beta))-as.vector(beta0))^2))
      }
    }
    
    output$x[,j] = as.vector(beta); 
    output$df[j]=length(unique(beta))
    output$k[j]=k
    output$objval[j] = h_objval[k]; 
    output$r_norm[j] = h_r_norm[k];
    output$s_norm[j] = h_s_norm[k];
    output$eps_pri[j] = h_eps_pri[k];
    output$eps_dual[j] = h_eps_dual[k];
  }
  
  # 5. report results
  return(output);  
}

admm.MN<-function(N, y, H, lambda.list, lr=0.01,
                    reltol=1e-04, abstol=1e-02/length(y), maxiter=1000, rho=1,
                    w=NULL, init=NULL, beta0=NULL){
  # 1. get parameters
  n=length(y); m=dim(H)[1]; B=m/n;
  if(is.null(w)){
    w=1
  }
  
  # 2. set ready
  if(is.null(init)){
    beta=runif(n)
  }else{
    beta=init
  }
  
  # 3. precompute static variables for x-update and factorization
  HtH = t(H)%*%H
  
  # 4. iteration
  sqrtn = sqrt(n);
  lambda.list.sort=sort(lambda.list,decreasing=TRUE)
  nlambda=length(lambda.list)
  output=list()
  output$lambda=lambda.list.sort
  output$x=matrix(0,n,nlambda)
  output$df=rep(0,nlambda)
  output$k=rep(0,nlambda)
  output$objval=rep(0,nlambda)
  output$r_norm=rep(0,nlambda)
  output$s_norm=rep(0,nlambda)
  output$eps_pri=rep(0,nlambda)
  output$eps_dual=rep(0,nlambda)
  
  for(j in 1:length(lambda.list.sort)){
    lambda=lambda.list.sort[j]
    beta=runif(n)
    z=H%*%beta
    u=rep(0,m)
    Hbeta = z
    
    h_objval=rep(0,maxiter)
    h_r_norm=rep(0,maxiter)
    h_s_norm=rep(0,maxiter)
    h_eps_pri=rep(0,maxiter)
    h_eps_dual=rep(0,maxiter)
    
    for(k in 1:maxiter){
      # 4-1. update 'x'
      b = exp(-beta)
      W = N*b/(1+b)^2
      D = Diagonal(x=W)+rho*HtH
      L = Cholesky(D, super=TRUE)
      q = -(y*b/(1+b)-(N-y)/(1+b))+rho*t(H)%*%(Hbeta-z+u)
      beta = beta-as.vector(solve(L,q))*lr
      
      # 4-2. update 'z'
      Hbeta = H%*%beta
      z_old = z;
      z = lasso_shrinkage(Hbeta + u, lambda/rho/B, w);
      
      # 4-3. update 'u'
      u = u + Hbeta - z;
      
      # 4-4. dianostics, reporting
      beta_norm = mean(abs(beta))
      z_norm = mean(abs(-z))
      # h_objval[k] = lasso_objective(X,y,v,Hbeta,lambda,beta);
      h_r_norm[k] = mean(abs(Hbeta-z));
      h_s_norm[k] = mean(abs(-rho*(z-z_old)));
      if (beta_norm>z_norm){
        h_eps_pri[k] = sqrtn*abstol + reltol*beta_norm;
      } else {
        h_eps_pri[k] = sqrtn*abstol + reltol*z_norm;
      }
      h_eps_dual[k] = sqrtn*abstol + reltol*mean(abs(rho*u));
      
      # 4-5. termination
      if ((h_r_norm[k] < h_eps_pri[k])&(h_s_norm[k]<h_eps_dual[k])){
        break;
      }
      
      print(k)
      if(is.null(beta0)==FALSE){
        print(mean((beta-as.vector(beta0))^2))
      }
    }
    
    output$x[,j] = as.vector(beta); 
    output$df[j]=length(unique(beta))
    output$k[j]=k
    output$objval[j] = h_objval[k]; 
    output$r_norm[j] = h_r_norm[k];
    output$s_norm[j] = h_s_norm[k];
    output$eps_pri[j] = h_eps_pri[k];
    output$eps_dual[j] = h_eps_dual[k];
  }
  
  # 5. report results
  return(output);  
}


Get_BIC <- function(fit, y, round.digits=3){
  n=length(y);nlambda=length(fit$lambda)
  output=list()
  output$BIC=output$k=output$dv=rep(0,nlambda)
  for(i in 1:nlambda){
    beta=round(fit$x[,i],digits=round.digits)
    beta=1/(1+exp(-beta))
    output$k[i]=length(unique(beta))
    output$dv[i]=-2*sum(y*log(beta)+(1-y)*log(1-beta))
    output$BIC[i]=output$dv[i]+output$k[i]*log(n)
  }
  return(output)
}

Delete_Duplicate<-function(data){
  coords=matrix(0,0,2);Y.new=c();
  for(i in 1:length(Y)){
    id=which((coords[,1]==coords[i,1])&(coords[,2]==coords[i,2]))
    if(!(min(id)<i)){
      coords.new=rbind(coords.new,coords[i,])
      tmp=mean(Y[id])
      if(tmp>=0.5){
        Y.new=c(Y.new,1)
      }else{
        Y.new=c(Y.new,0)
      }
    }
  }
  return(list(coords=coords.new,Y=Y.new))
}

# admm.Bern<-function(y, H, lambda.list, lr=0.01,
#                     reltol=1e-04, abstol=1e-02/length(y), maxiter=1000, rho=1,
#                     w=NULL, init=NULL, beta0=NULL){
#   # 1. get parameters
#   n=length(y); m=dim(H)[1]; B=m/n;
#   if(is.null(w)){
#     w=1
#   }
#   
#   # 2. set ready
#   if(is.null(init)){
#     beta=runif(n)
#   }else{
#     beta=init
#   }
#   
#   # 3. precompute static variables for x-update and factorization
#   HtH = t(H)%*%H
#   
#   # 4. iteration
#   sqrtn = sqrt(n);
#   lambda.list.sort=sort(lambda.list,decreasing=TRUE)
#   nlambda=length(lambda.list)
#   output=list()
#   output$lambda=lambda.list.sort
#   output$x=matrix(0,n,nlambda)
#   output$df=rep(0,nlambda)
#   output$k=rep(0,nlambda)
#   output$objval=rep(0,nlambda)
#   output$r_norm=rep(0,nlambda)
#   output$s_norm=rep(0,nlambda)
#   output$eps_pri=rep(0,nlambda)
#   output$eps_dual=rep(0,nlambda)
#   
#   for(j in 1:length(lambda.list.sort)){
#     lambda=lambda.list.sort[j]
#     #beta=runif(n)
#     z=H%*%beta
#     u=rep(0,m)
#     Hbeta = z
#     
#     h_objval=rep(0,maxiter)
#     h_r_norm=rep(0,maxiter)
#     h_s_norm=rep(0,maxiter)
#     h_eps_pri=rep(0,maxiter)
#     h_eps_dual=rep(0,maxiter)
#     
#     for(k in 1:maxiter){
#       # 4-1. update 'x'
#       W = y/beta^2+(1-y)/(1-beta)^2
#       D = Diagonal(x=W)+rho*HtH
#       L = Cholesky(D, super=TRUE)
#       q = (beta*W+y/beta-(1-y)/(1-beta))+rho*t(H)%*%(z-u)
#       beta = as.vector(solve(L,q))
#       
#       # 4-2. update 'z'
#       Hbeta = H%*%beta
#       z_old = z;
#       z = lasso_shrinkage(Hbeta + u, lambda/rho/B, w);
#       
#       # 4-3. update 'u'
#       u = u + Hbeta - z;
#       
#       # 4-4. dianostics, reporting
#       beta_norm = mean(abs(beta))
#       z_norm = mean(abs(-z))
#       # h_objval[k] = lasso_objective(X,y,v,Hbeta,lambda,beta);
#       h_r_norm[k] = mean(abs(Hbeta-z));
#       h_s_norm[k] = mean(abs(-rho*(z-z_old)));
#       if (beta_norm>z_norm){
#         h_eps_pri[k] = sqrtn*abstol + reltol*beta_norm;
#       } else {
#         h_eps_pri[k] = sqrtn*abstol + reltol*z_norm;
#       }
#       h_eps_dual[k] = sqrtn*abstol + reltol*mean(abs(rho*u));
#       
#       # 4-5. termination
#       if ((h_r_norm[k] < h_eps_pri[k])&(h_s_norm[k]<h_eps_dual[k])){
#         break;
#       }
#       
#       if(is.null(beta0)==FALSE){
#         print(k)
#         print(mean((beta-as.vector(beta0))^2))
#       }
#     }
#     
#     output$x[,j] = as.vector(beta); 
#     output$df[j]=length(unique(beta))
#     output$k[j]=k
#     output$objval[j] = h_objval[k]; 
#     output$r_norm[j] = h_r_norm[k];
#     output$s_norm[j] = h_s_norm[k];
#     output$eps_pri[j] = h_eps_pri[k];
#     output$eps_dual[j] = h_eps_dual[k];
#   }
#   
#   # 5. report results
#   return(output);  
# }

Wilson <- function(p,alpha,delta){
  q=qnorm(1-alpha/2)
  d=4*delta^2;a=2*p*(1-p)
  return((q^2*(-d+a+sqrt((d-a)^2-d*(d-1))))/d)
}

EL_Jeffrey <- function(n,p,alpha){
  sum(sapply(0:n, FUN=function(i){
    if(i==0){
      U=1-(alpha/2)^(1/n)
      L=0
    }else if(i==1){
      U=qbeta(1-alpha/2,2,n)
      L=0
    }else if(i==(n-1)){
      U=1
      L=qbeta(alpha/2,n,2)
    }else if(i==n){
      U=1
      L=(alpha/2)^(1/n)
    }else{
      U=qbeta(1-alpha/2,i+0.5,n-i+0.5)
      L=qbeta(alpha/2,i+0.5,n-i+0.5)
    }
    return(choose(n,i)*p^i*(1-p)^(n-i)*(U-L))
  }))/2
}

EL_Jeffrey_MC <- function(n,p,alpha){
  nsample=max(100*n,500000)
  a=rbinom(nsample,n,p)
  mean(sapply(a, FUN=function(i){
    if(i==0){
      U=1-(alpha/2)^(1/n)
      L=0
    }else if(i==1){
      U=qbeta(1-alpha/2,2,n)
      L=0
    }else if(i==(n-1)){
      U=1
      L=qbeta(alpha/2,n,2)
    }else if(i==n){
      U=1
      L=(alpha/2)^(1/n)
    }else{
      U=qbeta(1-alpha/2,i+0.5,n-i+0.5)
      L=qbeta(alpha/2,i+0.5,n-i+0.5)
    }
    return((U-L))
  }))/2
}





