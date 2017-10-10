FIS <- function(Y,Z,R,TT,Q,S){
  
  library(corpcor)
  
  # %______________________________________________________________________
  # % Fixed intervall smoother (see Harvey, 1989, p. 154)
  # % FIS returns the smoothed state vector AmT and its covar matrix PmT             
  # % Use this in conjunction with function SKF
  # %______________________________________________________________________
  # % INPUT  
  # %        Y         Data                                 (nobs x n)  
  # %        S Estimates from Kalman filter SKF                                                          
  # %          S.Am   : Estimates     a_t|t-1                  (nobs x m) 
  # %          S.Pm   : P_t|t-1 = Cov(a_t|t-1)             (nobs x m x m)
  # %          S.AmU  : Estimates     a_t|t                    (nobs x m) 
  # %          S.PmU  : P_t|t   = Cov(a_t|t)               (nobs x m x m)       
  # % OUTPUT 
  # %        S Smoothed estimates added to above
  # %          S.AmT  : Estimates     a_t|T                    (nobs x m) 
  # %          S.PmT :  P_t|T   = Cov(a_t|T)               (nobs x m x m)
  # %          S.PmT_1 : Cov(a_ta_t-1|T)
  # %        where m is the dim of state vector and t = 1 ...T is time
  
  m<-dim(S$Am)[1]
  nobs<-dim(S$Am)[2]
  
  S$AmT           = zeros(m,nobs+1)
  S$PmT           = array(0,c(m,m,nobs+1))
  S$AmT[,nobs+1] <- S$AmU[,nobs+1]
  S$PmT[,,nobs+1] <- S$PmU[,,nobs+1]
  S$PmT_1[,,nobs] <- (eye(m)-S$KZ)%*%TT%*%S$PmU[,,nobs]
  
  pinv<-corpcor::pseudoinverse(S$Pm[,,nobs])
  
  J_2 <- S$PmU[,,nobs]%*%t(TT)%*%pinv
  
  for (t in nobs:1){ 
    PmU <- S$PmU[,,t]
    Pm1 <- S$Pm[,,t]
    P_T <- S$PmT[,,t+1]
    P_T1 <- S$PmT_1[,,t]
    
    J_1 <- J_2
    
    S$AmT[,t] <- S$AmU[,t] + J_1%*%(S$AmT[,t+1] - TT%*%S$AmU[,t])
    S$PmT[,,t] <- PmU + J_1%*%(P_T - Pm1)%*%t(J_1) 
    
    if(t>1){
      pinv<-corpcor::pseudoinverse(S$Pm[,,t-1])
      J_2 <- S$PmU[,,t-1]%*%t(TT)%*%pinv
      S$PmT_1[,,t-1] = PmU%*%t(J_2)+J_1%*%(P_T1-TT%*%PmU)%*%t(J_2)
    }
  }
  
  return(S)
  
}

