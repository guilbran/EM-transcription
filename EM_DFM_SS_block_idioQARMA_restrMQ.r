EM_DFM_SS_block_idioQARMA_restrMQ<-function(X,Par){
  
  library(matlab)
  
  thresh = 1e-4
  r = Par$r
  p = Par$p
  max_iter = Par$max_iter
  i_idio = Par$i_idio
  R_mat = Par$Rconstr
  q = Par$q
  nQ = Par$nQ
  blocks = Par$blocks
  
  ### Prepara??o dos dados
  
  TT <- dim(X)[1]
  N <- dim(X)[2]
  
  ### Standardise X
  Mx = colMeans(X,na.rm=T)
  Wx = sapply(1:N,function(x) sd(X[,x],na.rm = T))
  xNaN = (X-repmat(Mx,TT,1))/repmat(Wx,TT,1)
  
  ### Initial conditions
  
  # Removing missing values (for initial estimator)
  optNaN<-list()
  optNaN$method = 2; # Remove leading and closing zeros
  optNaN$k = 3;
  
  res_InitCond<-InitCond(xNaN,r,p,blocks,optNaN,R_mat,q,nQ,i_idio);
  
  A<-res_InitCond$A
  C<-res_InitCond$C
  Q<-res_InitCond$Q
  R<-res_InitCond$R
  Z_0<-res_InitCond$initZ
  V_0<-res_InitCond$initV
  
  # some auxiliary variables for the iterations
  previous_loglik = -Inf
  num_iter = 0
  LL = -Inf
  converged = F
  
  # y for the estimation is WITH missing data
  y = t(xNaN)
  
  #--------------------------------------------------------------------------
  #THE EM LOOP
  #--------------------------------------------------------------------------
  
  #The model can be written as
  #y = C*Z + e;
  #Z = A*Z(-1) + v
  #where y is NxT, Z is (pr)xT, etc
  
  #remove the leading and ending nans for the estimation
  optNaN$method = 3
  y_est <- remNaNs_spline(xNaN,optNaN)
  y_est_indNaN<-t(y_est$indNaN)
  y_est<-t(y_est$X)
  
  
  while ((num_iter < max_iter) & !converged){
    
    res_EMstep = EMstep(y_est, A, C, Q, R, Z_0, V_0, r,p,R_mat,q,nQ,i_idio,blocks)
    # res_EMstep <- list(C_new, R_new, A_new, Q_new, Z_0, V_0, loglik)
    
    C = res_EMstep$C_new;
    R = res_EMstep$R_new;
    A = res_EMstep$A_new;
    Q = res_EMstep$Q_new;
    
    # Checking convergence
    if (num_iter>2){
      res_em_converged = em_converged(loglik, previous_loglik, thresh,1)
      # res_em_converged<-list(converged,decrease[num_iter+1])
      
      converged<-res_em_converged$converged
      decreasse<-res_em_converged$decrease
    }
    
    LL <- cbind(LL, loglik)
    previous_loglik <- loglik
    num_iter <-  num_iter + 1
  }
  
  # final run of the Kalman filter
  Zsmooth = runKF(y, A, C, Q, R, Z_0, V_0)
  Zsmooth<-t(Zsmooth)
  x_sm <- Zsmooth[2:dim(Zsmooth)[1],]%*%t(C)
  
  Res<-list()
  
  Res$X_sm <- repmat(Wx,TT,1)*x_sm+repmat(Mx,TT,1)
  Res$FF <- Zsmooth[2:dim(Zsmooth)[1],]
  
  #--------------------------------------------------------------------------
  #  Loading the structure with the results
  #--------------------------------------------------------------------------
  Res$C <- C;
  Res$R <- R;
  Res$A <- A;
  Res$Q <- Q;
  Res$Mx <- Mx;
  Res$Wx <- Wx;
  Res$Z_0 <- Z_0;
  Res$V_0 <- V_0;
  Res$r <- r;
  Res$p <- p;
  
  return(Res)
  
  
}
