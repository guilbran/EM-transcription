em_converged <- function(loglik = NULL, previous_loglik = NULL, threshold = NULL, check_increased = NULL){
  # EM_CONVERGED Has EM converged?
  # [converged, decrease] = em_converged(loglik, previous_loglik, threshold)
  #
  # We have converged if the slope of the log-likelihood function falls below 'threshold',
  # i.e., |f(t) - f(t-1)| / avg < threshold,
  # where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
  # 'threshold' defaults to 1e-4.
  #
  # This stopping criterion is from Numerical Recipes in C p423
  #
  # If we are doing MAP estimation (using priors), the likelihood can decrase,
  # even though the mode of the posterior is increasing.
  
  nargin <- 4 - sum(is.null(loglik, previous_loglik, threshold, check_increased))
  
  if(nargin < 3){threshold <- 1e-4}
  if(nargin < 4){check_increased <- 1}
  
  converged <- 0
  decrease <- 0
  
  if(!is.null(check_increased)){
    if(loglik - previous_loglik < -1e-3){ # allow for a little imprecision
      print(paste(1, '******likelihood decreased from %6.4f to %6.4f!\n', previous_loglik, loglik))
      decrease <- 1
    }
  }
  
  delta_loglik <- abs(loglik - previous_loglik)
  avg_loglik <- (abs(loglik) + abs(previous_loglik) + 2.2204e-16)/2
  
  if((delta_loglik / avg_loglik) < threshold){converged <- 1}
  
  # output
  list(converged = converged, decrease = decrease)
  
}


