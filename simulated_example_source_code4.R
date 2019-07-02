generate_weights<-function(Smod,Amod, data){
  S1data<-subset(data, S==1)
  w_reg<-glm(Smod, family="binomial", data=data)
  ps<- predict(w_reg,newdata=data, type="response") 
  w_reg2<-glm(Amod, family="binomial", data=S1data)
  pa<- predict(w_reg2,newdata=data, type="response") 
  w= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
  data$w<-w
  list<-list(dat=data, Smod=w_reg, Amod=w_reg2)
  return(list)
}


#functions that do not account for estimation of the working models 
#most appropriate for bootstrap estimation of SEs

OM_est<-function(data){
  S1data_A1<-subset(data, S==1 & A==1)
  OM1mod<-glm(formula=Y~X1+X2+X3, data=S1data_A1)
  p1<- predict(OM1mod,newdata=data, type="response") 
  data$p1<-p1
  S1data_A0<-subset(data, S==1 & A==0)
  OM0mod<-glm(formula=Y~X1+X2+X3, data=S1data_A0)
  p0<- predict(OM0mod,newdata=data, type="response") 
  data$p0<-p0
  S0sub<-subset(data, S==0)
  OM_1<-mean(S0sub$p1)
  OM_0<-mean(S0sub$p0)
  OM<-mean(S0sub$p1)-mean(S0sub$p0)
  list<-list(OM_1=OM_1, OM_0=OM_0, OM=OM, p1=p1, p0=p0, OM1mod=OM1mod, OM0mod=OM0mod)
  return(list)
}


IOW1_est<-function(data){ 
  A<-data$A
  S<-data$S
  w<-data$w
  Y<-data$Y
  IOW1_1 <-(sum((1-S))^-1)* sum(A*S*w*Y)
  IOW1_0 <-(sum((1-S))^-1)* sum((1-A)*S*w*Y)
  IOW1 = IOW1_1 - IOW1_0
  return(list(IOW1_1=IOW1_1,IOW1_0=IOW1_0, IOW1=IOW1))
}


IOW2_est<-function(data){
  S0data<-subset(data, S==0)  
  S1data_A1<-subset(data, S==1 & A==1)
  IOW1mod<-glm(formula=Y~1, data=S1data_A1, weights=w)
  p1<- predict(IOW1mod,newdata=S0data, type="response") 
  S1data_A0<-subset(data, S==1 & A==0)
  IOW0mod<-glm(formula=Y~1, data=S1data_A0, weights=w)
  p0<- predict(IOW0mod,newdata=S0data, type="response") 
  IOW2_1<-mean(p1)
  IOW2_0<-mean(p0)
  IOW2<-mean(p1)-mean(p0)
  list<-list(IOW2_1=IOW2_1,IOW2_0=IOW2_0, IOW2=IOW2,IOW1mod=IOW1mod,IOW0mod=IOW0mod)
  return(list)
}


DR1_est<-function(data){  
  A<-data$A
  S<-data$S
  Y<-data$Y
  p1<-data$p1
  p0<-data$p0
  w<-data$w
  DR1_1<-(sum((1-S))^-1)* sum(S*A*w*(Y-p1) + (1-S)*p1)
  DR1_0<-(sum((1-S))^-1)* sum(S*(1-A)*w*(Y-p0) + (1-S)*p0)
  DR1<-DR1_1-DR1_0
  list<-list(DR1_1=DR1_1,DR1_0=DR1_0, DR1=DR1)
  return(list)
}


DR2_est<-function(data){
  A<-data$A
  S<-data$S
  Y<-data$Y
  p1<-data$p1
  p0<-data$p0
  w<-data$w
  sum1_DR2<-sum(S*A*w*(Y-p1)) 
  sum0_DR2<-sum(S*(1-A)*w*(Y-p0)) 
  norm1<-(sum(S*A*w))^-1
  norm0<-(sum(S*(1-A)*w))^-1
  DR2_1<-norm1*sum1_DR2 + (sum(1-S)^-1)*sum((1-S)*p1)
  DR2_0<-norm0*sum0_DR2 + (sum(1-S)^-1)*sum((1-S)*p0)
  DR2<-DR2_1-DR2_0
  list<-list(DR2_1=DR2_1,DR2_0=DR2_0, DR2=DR2)
  return(list)
}


DR3_est<-function(data){
  S0data<-subset(data, S==0)
  S1data_A1<-subset(data, S==1 & A==1)
  DR1mod<-glm(formula=Y~X1+X2+X3, data=S1data_A1, weights=w)
  p1<- predict(DR1mod,newdata=S0data, type="response") 
  S1data_A0<-subset(data, S==1 & A==0)
  DR0mod<-glm(formula=Y~X1+X2+X3, data=S1data_A0, weights=w)
  p0<- predict(DR0mod,newdata=S0data, type="response") 
  DR3_1<- mean(p1)
  DR3_0<- mean(p0)
  DR3<-mean(p1)-mean(p0) 
  list<-list(DR3_1=DR3_1,DR3_0=DR3_0, DR3=DR3,DR1mod=DR1mod, DR0mod=DR0mod)
  return(list)
}



#functions for M-estimation (geex)

OM_EE <- function(data){
  A<-data$A
  S<- data$S
  Y <- data$Y
  X <- cbind(1, data$X1, data$X2, data$X3) 
  matA <- cbind(1, data$A) 
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  function(theta){ 
    #outcome model 
    beta<-theta[1:4]
    alpha<-theta[5:8]
    mu1<-theta[9]
    mu0<-theta[10]
    muate<-theta[11]
    m_A1 <-X %*% beta
    m_A0<-X %*% alpha
    ols_A1 <-crossprod(X, (S*A)*(Y - m_A1))
    ols_A0 <-crossprod(X, (S*(1-A))*(Y - m_A0))
    #estimates
    mean1<-(1-S)*(m_A1-mu1) 
    mean0 <- (1-S)*(m_A0-mu0) 
    ate<-(1-S)*(m_A1-m_A0-muate) 
    c(ols_A1,ols_A0,mean1, mean0,ate)
  }
}


IOW1_EE <- function(data){ 
  A<-data$A
  S<- data$S
  Y <- data$Y  
  X <- cbind(1, data$X1, data$X2, data$X3) 
  matA <- cbind(1, data$A) 
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  function(theta){
    mu1<-theta[9]
    mu0<-theta[10]
    mu_ate<-theta[11]
    #participation model
    lp  <- X %*% theta[1:4]
    ps <- plogis(lp)
    score_eqns<-crossprod(X, S-ps)
    #treatment model
    lp2  <- X %*% theta[5:8] 
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,S*(A - pa) )
    w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
    #estimates
    summand1<-(w*A*S*Y)  
    summand0<-(w*(1-A)*S*Y)  
    mean1<-(summand1)- (1-S)*mu1
    mean0<-(summand0)- (1-S)*mu0
    ate<-(summand1)-(summand0)- (1-S)*mu_ate
    c(score_eqns,score_eqns2,mean1,mean0, ate)
  }
}


IOW2_EE <- function(data){
  A<-data$A
  S<- data$S
  Y <- data$Y
  X <- cbind(1, data$X1, data$X2, data$X3)
  matA <- cbind(1, data$A) 
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  function(theta){
    #participation model
    lp  <- X %*% theta[1:4]
    ps <- plogis(lp)
    score_eqns<-crossprod(X, S-ps)
    #treatment model
    lp2  <- X %*% theta[5:8] 
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,S*(A - pa) )
    w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
    #outcome model
    m_A1<-1 %*% theta[9]
    m_A0<-1 %*% theta[10]
    linear_eqns1<-crossprod(1, (S*A*w)*(Y -  m_A1) )
    linear_eqns0<-crossprod(1, (S*(1-A)*w)*(Y - 1 %*% theta[10]) )
    mu1<-theta[11]
    mu0<-theta[12]
    muate<-theta[13]
    #estimates
    mean1 <- (1-S)*( m_A1 -mu1) 
    mean0 <- (1-S)*( m_A0 -mu0)
    ate <- (1-S)*(m_A1- m_A0 - muate)
    c(score_eqns,score_eqns2, linear_eqns1,linear_eqns0,  mean1,  mean0,ate)
  }
}


DR1_EE <- function(data){
  A<-data$A
  S<- data$S
  Y <- data$Y
  X <- cbind(1, data$X1, data$X2, data$X3) 
  matA <- cbind(1, data$A)  
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  function(theta){
    #participation model
    lp  <- X %*% theta[1:4]
    ps <- plogis(lp)
    score_eqns<-crossprod(X, S-ps)
    #treatment model
    lp2  <- X %*% theta[5:8] 
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,S*(A - pa) )
    w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
    #outcome model
    beta<-theta[9:12]
    alpha<-theta[13:16]
    mu1<-theta[17]
    mu0<-theta[18]
    mu<-theta[19]
    m_A1 <-X %*% beta
    m_A0<-X %*% alpha
    ols_A1 <-crossprod(X, (S*A)*(Y - m_A1))
    ols_A0 <-crossprod(X, (S*(1-A))*(Y - m_A0))
    ey1<-w*S*A*(Y-m_A1) + (1-S)*m_A1
    ey0<-w*S*(1-A)*(Y-m_A0) + (1-S)*m_A0
    #estimates
    mean1<-ey1-(1-S)*mu1
    mean0<-ey0-(1-S)*mu0
    ate<-ey1-ey0- (1-S)*mu
    c(score_eqns,score_eqns2,ols_A1,ols_A0,mean1,mean0, ate)   
  }
}


DR2_EE <- function(data){
  A<-data$A
  S<- data$S
  Y <- data$Y
  X <- cbind(1, data$X1, data$X2, data$X3) 
  matA <- cbind(1, data$A) 
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0
  function(theta){
    #participation model
    lp  <- X %*% theta[1:4]
    ps <- plogis(lp)
    score_eqns<-crossprod(X, S-ps)
    #treatment model
    lp2  <- X %*% theta[5:8] 
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,S*(A - pa) ) 
    w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa)) 
    #outcome model 
    beta<-theta[9:12]
    alpha<-theta[13:16]
    mu1<-theta[20]
    mu0<-theta[21]
    mu<-theta[22]
    m_A1 <-X %*% beta
    m_A0<-X %*% alpha
    ols_A1 <-crossprod(X, (S*A)*(Y - m_A1))
    ols_A0 <-crossprod(X, (S*(1-A))*(Y - m_A0))
    #normalizing term
    mu_S<-theta[17]
    propS1<-S-mu_S
    one_over<-(1/(1-mu_S))
    mu_norm1<-theta[18]
    norm1eq<-(A*S*w)-mu_norm1
    norm1<-1/mu_norm1
    mu_norm0<-theta[19]
    norm0eq<-((1-A)*S*w)-mu_norm0
    norm0<-1/mu_norm0
    ey1<-norm1*((w*S*A*(Y-m_A1))) + one_over*((1-S)*m_A1)
    ey0<-norm0*((w*S*(1-A)*(Y-m_A0))) + one_over*((1-S)*m_A0)
    #estimates
    mean1<-ey1-mu1
    mean0<-ey0-mu0
    ate<-ey1-ey0-mu
    c(score_eqns,score_eqns2,ols_A1,ols_A0,propS1, norm1eq, norm0eq,mean1,mean0,ate)
  }
}


DR3_EE <- function(data){
  A<-data$A
  S<- data$S
  Y <- data$Y
  X <- cbind(1, data$X1, data$X2, data$X3) 
  matA <- cbind(1, data$A) 
  A[is.na(A)] <- 0 
  Y[is.na(Y)] <- 0 
  function(theta){
    #participation model
    lp  <- X %*% theta[1:4]
    ps <- plogis(lp)
    score_eqns<-crossprod(X, S-ps)
    #treatment model
    lp2  <- X %*% theta[5:8] 
    pa<- plogis(lp2)
    score_eqns2<-crossprod(X,S*(A - pa) )
    w = (A * S*(1-ps))/(ps*pa) + ((1 - A)*S*(1-ps))/(ps*(1-pa))
    #outcome model
    beta<-theta[9:12]
    alpha<-theta[13:16]
    mu1<-theta[17]
    mu0<-theta[18]
    mu<-theta[19]
    m_A1 <-X %*% beta
    m_A0<-X %*% alpha
    ols_A1 <-crossprod(X, (S*A*w)*(Y - m_A1))
    ols_A0 <-crossprod(X, (S*(1-A)*w)*(Y - m_A0))
    #estimates   
    mean1 <- (1-S)*(m_A1-mu1) 
    mean0 <- (1-S)*(m_A0-mu0) 
    mean <- (1-S)*(m_A1-m_A0-mu) 
    c(score_eqns,score_eqns2,ols_A1,ols_A0,mean1,mean0, mean)
  }
}


#Function to extract point estimate and SE from geex output
extractEST<-function(geex_output=OM_mest, est_name="m1",param_start=param_start_OM){
  param_num_EST<-match(est_name,names(param_start))
  EST<-geex_output@estimates[param_num_EST]
  sandwich_se <- diag(geex_output@vcov)^0.5 
  SE<-sandwich_se[param_num_EST]
  return(c(EST, SE=SE))
}
