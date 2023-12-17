########### CLASSICAL ANALYSIS ###############
# Install and load the necessary packages
install.packages("MASS")

library(MASS)

rm(list = ls(all=TRUE))

# Set parameters
theta0=0.5

phi_true <- 0.8     # True AR coefficient

sigma2_true <- 4    # True standard deviation of the AR process

n <- 100            # Number of observations per replication

num_replications <- 10  # Number of replications

df <- 8      # Degrees of freedom for t-distribution

################################################################################################
############## Create a matrix to store the simulated data

simulated_data <- matrix(NA, nrow = n, ncol = num_replications)


# Create vectors to store estimated parameters
est_theta0 <- numeric(num_replications)

est_phi <- numeric(num_replications)

est_df <- numeric(num_replications)

est_sigma2 <- numeric(num_replications)

################################################################
# Loop over replications
for (rep in 1:num_replications) {

    set.seed(rep)# Set seed for each replication
    

  errors <- rt(n, df)  # Generate t-distributed errors
  
  # Simulate AR(1) process
  ar_process <- numeric(n)
  
ar_process0=0

  ar_process[1] <-theta0+(phi_true*ar_process0) + errors[1]  # Initial value
  
  for (i in 2:n) {
  
    ar_process[i] <- theta0+(phi_true * ar_process[i-1]) + errors[i]
  }
 
  simulated_data[, rep] <- ar_process
  
  # Estimate AR(1) parameters using MLE
  
  #ar_fit <- arima(ar_process, order = c(1, 0, 0), method = "ML")
  
loglike=function(parm){

    theta0=parm[1];phi_true=parm[2];df=parm[3];sigma2_true=parm[4]
    
    
    eps=rep(0,n)
    
    for(i in 2:n){
    
    for(j in 1:rep){
    
      eps[i]=simulated_data[,j][i]-theta0-(phi_true*simulated_data[,j][i-1])
    }
}
    L=log((gamma(df+1)/2)/(sqrt(df)*sqrt(pi)*gamma(df/2)*sqrt(sigma2_true)))
    
    loglike=(n-1)*L-((df+1)/2)*sum(log(1+eps^2/(df*sigma2_true)))
    
    return(-loglike)
  }

  hess=nlm(loglike,c(0.5,0.8,8,0.0002),hessian=T)$hessian
  
  MLE=nlm(loglike,c(0.5,0.8,8,0.0002),hessian=T)$estimate
  
  theta0=MLE[1];phi_true=MLE[2];df=MLE[3];sigma2_true=MLE[4]
  
  cov=solve(hess)
  
  
est_theta0[rep]=theta0

  est_phi[rep] <- phi_true
  
est_df[rep]=df

  est_sigma2[rep] <- sigma2_true
}

######################################################### BAYESIAN ANALYSIS ###############

library(coda)

pmode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}
set.seed(1234)

n=80

theta0=0.5

theta1=0.8

df=8

sigma2=4

errors <- rt(n, df)  # Generate t-distributed errors


# Simulate AR(1) process
ar_process <- numeric(n)

ar_process0=0

ar_process[1] <-theta0+(theta1*ar_process0) + errors[1]  # Initial value


for (i in 2:n) {

  ar_process[i] <- theta0+(theta1 * ar_process[i-1]) + errors[i]
  
}

n=length(ar_process)
  
  # Estimate AR(1) parameters using MLE
  #ar_fit <- arima(ar_process, order = c(1, 0, 0), method = "ML")
  
  loglike=function(parm){
  
    theta0=parm[1];theta1=parm[2];df=parm[3];sigma2=parm[4]
    
    
    eps=rep(0,n)
    
    for(i in 2:n){
    
      eps[i]=ar_process[i]-theta0-theta1*ar_process[i-1]
    }
    L=log((gamma(df+1)/2)/(sqrt(df)*sqrt(pi)*gamma(df/2)*sqrt(sigma2)))
    
    loglike=(n-1)*L-((df+1)/2)*sum(log(1+eps^2/(df*sigma2)))
    
    return(-loglike)
  }
  
  hess=nlm(loglike,c(0.5,0.8,8,0.0002),hessian=T)$hessian
  
  MLE=nlm(loglike,c(0.5,0.8,8,0.0002),hessian=T)$estimate
  
  theta0=MLE[1];theta1=MLE[2];df=MLE[3];sigma2=MLE[4]
  
  cov=solve(hess)
  
  
c=0.7;s1=0.51;s2=0.1;s3=0.9;s4=0.5;Niter=20000;alpha=250;beta=2;

MCMC=matrix(0,nrow=Niter,ncol=4);MCMC[1,]=c(theta0,theta1,df,sigma2);

k1=0;k2=0;k3=0;k4=0

start.time <- Sys.time()

for(i in 2:Niter){

theta1=MCMC[i-1,2];df=MCMC[i-1,3];sigma2=MCMC[i-1,4];

fcond1=function(theta0)

{

eps=rep(0,n)
for(j in 2:n){

eps[j]=ar_process[j]-theta0-(theta1*ar_process[j-1])

}
fcond1=prod(1+(eps^2/(df*(sigma2))))^(-df-1)/2

return(fcond1)

}
repeat{

proptheta0=rnorm(1,mean=theta0,sd=c*s1)

if((proptheta0<100)&&(proptheta0>-100))break();

}
r1=fcond1(proptheta0)/fcond1(theta0)

if(runif(1)<r1){MCMC[i,1]=proptheta0;k1=k1+1} else{MCMC[i,1]=theta0}

theta0=MCMC[i,1]

#----------------------End of Full Conditional theta0---------------------###
fcond2=function(theta1)

{
eps=rep(0,n)

for(j in 2:n){

eps[j]=ar_process[j]-theta0-(theta1*ar_process[j-1])

}
fcond2=prod(1+(eps^2/(df*(sigma2))))^(-df-1)/2

return(fcond2)
}

repeat{
proptheta1=rnorm(1,mean=theta1,sd=c*s2)

if((proptheta1<100)&&(proptheta1>-100))break();

}
r2=fcond2(proptheta1)/fcond2(theta1)

if(runif(1)<r2){MCMC[i,2]=proptheta1;k2=k2+1} else{MCMC[i,2]=theta1}

theta1=MCMC[i,2]

#--------------------End of Full Conditional theta1------------------------

fcond3=function(df){

eps=rep(0,n)

for(j in 2:n){

eps[j]=ar_process[j]-theta0-(theta1*ar_process[j-1])
}

fcond3=exp(-beta*df)*(df^(alpha-1-n/2))*((gamma((df+1)/2)/(gamma(df/2)*sqrt(pi)))^n)*(prod(1+(eps^2/(df*(sigma2))))^(-df-1)/2)

return(fcond3)
}

repeat{
propnu=rnorm(1,mean=df,sd=c*s3)

if((propnu<8.1)&&(propnu>3))break();
}

r3=fcond3(propnu)/fcond3(df)

if(runif(1)<r3){MCMC[i,3]=propnu;k3=k3+1} else{MCMC[i,3]=df}

df=MCMC[i,3]


#--------------------End of Full Conditional nu------------------------

fcond4=function(sigma2)
{

eps=rep(0,n)

for(j in 2:n){

eps[j]=ar_process[j]-theta0-(theta1*ar_process[j-1])
}
fcond4=(sqrt(sigma2)^(-n+1))*(prod(1+(eps^2/(df*(sigma2))))^(-df-1)/2)

return(fcond4)
}

propsig=rnorm(1,mean=sigma2,sd=c*s4)

r4=fcond4(propsig)/fcond4(sigma2)

if(runif(1)<r4){MCMC[i,4]=propsig;k4=k4+1} else{MCMC[i,4]=sigma2}

sigma2=MCMC[i,4];

MC=matrix(0, i, 2); 

if(i==1000) {MC=MCMC[1:i, 3:4]; s3=sqrt(var(MC[,1])); s4=sqrt(var(MC[,2]))}
}

end.time <- Sys.time()

time.taken <- end.time - start.time

time.taken

#####################################################
