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
  
  xsmooth <- S$AmT;
  Vsmooth <- S$PmT;
  VVsmooth <- S$PmT_1;
  loglik <- S$loglik;
  
  return(list(xsmooth = xsmooth,Vsmooth = Vsmooth,VVsmooth = VVsmooth,loglik = loglik))  

}



