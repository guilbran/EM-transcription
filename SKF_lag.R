SKF_lag <-function(Y,Z,R,TT,Q,A_0,P_0){
  
  
  # %______________________________________________________________________
  # % Kalman filter for stationary systems with time-varying system matrices
  # % and missing data.
  # %
  # % The model is        y_t   = Z * a_t + eps_t       
  # %                     a_t+1 = TT * a_t + u_t       
  # %
  # %______________________________________________________________________
  # % INPUT  
  # %        Y         Data                                 (nobs x n)  
  # % OUTPUT 
  # %        S.Am       Predicted state vector  A_t|t-1      (nobs x m)  
  # %        S.AmU      Filtered  state vector  A_t|t        (nobs+1 x m)  
  # %        S.Pm       Predicted covariance of A_t|t-1      (nobs x m x m)  
  # %        S.PmU      Filtered  covariance of A_t|t        (nobs+1 x m x m)  
  # %        S.loglik   Value of likelihood function
  # 
  # % Output structure & dimensions
  
  n <- dim(Z)[1]
  m <- dim(Z)[2]
  nobs  <- size(Y,2)
  
  S<-list()
  S$Am <- array(NA,c(m,nobs))
  S$Pm <- array(NA,c(m,m,nobs))
  S$AmU <- array(NA,c(m,nobs+1))
  S$PmU <- array(NA,c(m,m,nobs+1))
  S$loglik <- 0
  
  
  # ______________________________________________________________________
  Au <- A_0;  # A_0|0;
  Pu <- P_0;  # P_0|0
  
  S$AmU[,1]    = Au;
  S$PmU[,,1]  = Pu;
  
  
  
  for(t in 1:nobs){
    #       t
    # A = A_t|t-1   & P = P_t|t-1
    
    A <- TT%*%Au;
    P <- TT%*%Pu%*%t(TT) + Q;
    P <-  0.5*(P+t(P))
    
    # handling the missing data
    res_MissData <- MissData(Y[,t],Z,R)
    
    y_t <- res_MissData$y
    Z_t <- res_MissData$C
    R_t <- res_MissData$R
    L_t <- res_MissData$L
    
    # if(is.null(y_t)){
    if(sum(is.na(y_t))==length(y_t)){
      Au <- A
      Pu <- P
    } else {
      
      
      if(!is.matrix(Z_t)){
        Z_t<-t(as.matrix(Z_t))
      }
      
      
      PZ  <- P%*%t(Z_t)
      iF  <- ginv(Z_t%*%PZ + R_t)
      PZF <- PZ%*%iF
      
      V <- y_t - Z_t%*%A
      Au <- A  + PZF%*%V
      Pu <- P  - PZF%*%t(PZ)
      Pu <-  0.5*(Pu+t(Pu))
      
    }
    
    S$Am[,t] <- A
    S$Pm[,,t] <- P
    
    # Au = A_t|t   & Pu = P_t|t
    
    S$AmU[,t+1] <- Au
    S$PmU[,,t+1] <- Pu
    S$loglik <- S$loglik + 0.5*(log(det(iF))  - t(V)%*%iF%*%V)
  } # t
  
  if(sum(is.na(y_t))==length(y_t)){
    S$KZ <- zeros(m,m)
  }else{
    S$KZ <- PZF%*%Z_t
  }
  
  return(S)
  
}