# %%%  Replication files for:
#   %%%  ""Nowcasting", 2010, (by Marta Banbura, Domenico Giannone and Lucrezia Reichlin), 
# %%% in Michael P. Clements and David F. Hendry, editors, Oxford Handbook on Economic Forecasting.
# %%%
# %%% The software can be freely used in applications. 
# %%% Users are kindly requested to add acknowledgements to published work and 
# %%% to cite the above reference in any resulting publications
# %--------------------------------------------------------------------------
# % KALMAN FILTER
# %--------------------------------------------------------------------------

runKF <- function(y, A, C, Q, R, x_0, Sig_0){
  S <- SKF(y,C,R,A,Q, x_0, Sig_0);
  S <- FIS(y,C,R,A,Q,S);
  
  return(S)  

}

SKF <-function(Y,Z,R,TT,Q,A_0,P_0){
# %______________________________________________________________________
# % Kalman filter for stationary systems with time-varying system matrices
# % and missing data.
# %
# % The model is        y_t   = Z * a_t + eps_t       
# %                     a_t+1 = T * a_t + u_t       
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

  S$Am <- array(NA,c(m,nobs))
  S$Pm <- array(NA,c(m,m,nobs))
  S$AmU <- array(NA,c(m,nobs+1))
  S$PmU <- array(NA,c(m,m,nobs+1))
  S$loglik <- 0
  

# %______________________________________________________________________
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

  if(is.null(y_t)){
    Au <- A
    Pu <- P
  } else {
  
    PZ  <- P%*%t(Z_t)
    iF  <- solve(Z_t%*%PZ + R_t)
    PZF <- PZ%*%iF
  
    V <- y_t - Z_t%*%A
    Au <- A  + PZF%*%V
    Pu <- P  - PZF%*%t(PZ)
    Pu <-  0.5%*%(Pu+t(Pu))
    S$loglik <- S$loglik + 0.5*(log(det(iF))  - t(V)%*%iF%*%V)
  }

  S$Am[,t] <- A
  S$Pm[,,t] <- P

  # Au = A_t|t   & Pu = P_t|t

  S$AmU[,t+1] <- Au
  S$PmU[,,t+1] <- Pu
} # t

  if(is.null(y_t)){
  S$KZ <- zeros(m,m)
  }else{
  S$KZ <- PZF%*%Z_t
  }

  return(S)

}




function S = FIS(Y,Z,R,TT,Q,S){

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
  
  S$AmT           = zeros(m,nobs+1))
  S$PmT           = array(0,c(m,m,nobs+1))
  S$AmT[,nobs+1] <- S$AmU[,nobs+1]
  S$PmT[,,nobs+1] <- S$PmU[:,:,nobs+1]
  S$PmT_1[,,nobs] <- (eye(m)-S$KZ)%*%TT%*%S%PmU[,,nobs]

  J_2 <- S$PmU[,,nobs]%*%t(TT)%*%corpcor::pseudoinverse(S$Pm[,,nobs])

  for (t in nobs:1){ 
  PmU <- S$PmU[,,t]
  Pm1 <- S$Pm[,,t]
  P_T <- S$PmT[,,t+1]
  P_T1 <- S$PmT_1[,,t]
  
  J_1 <- J_2
  
  S$AmT[,t] <- S$AmU[,t] + J_1%*%(S.AmT[,t+1] - TT%*%S$AmU[,t]) ; 
  S$PmT[,,t] <- PmU +  J_1%*%(P_T - Pm1)%*%t(J_1) 
  
  if(t>1){
  J_2 <- S$PmU[,,t-1]%*%t(TT)%*%pseudoinverse(S$Pm[,,t-1])
  S$PmT_1[,,t-1] = PmU%*%t(J_2)+J_1%*%(P_T1-TT%*%PmU)%*%t(J_2)
  }
}

%______________________________________________________________________
function [y,C,R,L]  = MissData(y,C,R);
%______________________________________________________________________
% PROC missdata                                                        
% PURPOSE: eliminates the rows in y & matrices Z, G that correspond to     
%          missing data (NaN) in y                                                                                  
% INPUT    y             vector of observations at time t  (n x 1 )    
%          S             KF system matrices             (structure)
%                        must contain Z & G
% OUTPUT   y             vector of observations (reduced)   (# x 1)     
%          Z G           KF system matrices     (reduced)   (# x ?)     
%          L             To restore standard dimensions     (n x #)     
%                        where # is the nr of available data in y
%______________________________________________________________________
ix = ~isnan(y);
e  = eye(size(y,1));
L  = e(:,ix);

y  =    y(ix);
C  =  C(ix,:);  
R  =  R(ix,ix);
