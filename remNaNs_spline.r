


remNaNs_spline <-function(X,options){
  
  TT <- dim(X)[1]
  N <- dim(X)[2]
  k <- options$k
  indNaN <- is.na(X)
  
  if(options$method == 1){ # replace all the missing values (método Giannone et al. 2008)
    for (i in 1:N){  
    x = X[,i]
    x[indNaN[,i]] = median(x,na.rm = T);
    x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
    x_MA=x_MA[(2*k+1):length(x_MA)]
    x[indNaN[,i]]=x_MA[indNaN[,i]]
    X[,i]=x;
    }
    
  }else if(options$method == 2){ # replace missing values after removing leading and closing zeros
    
    rem1 <- (rowSums(indNaN)>N*0.8)
    nanLead <- which(rem1)
    # nanEnd <- which(rem1[length(rem1):1])
    # nanLE <- c(nanEnd,nanLead)
    nanLE<-nanLead
    X[nanLE,] <- NA
    indNaN=is.na(X)
    for (i in 1:N){  
      x = X[,i]
      isnanx = is.na(x)
      t1 = min(which(!isnanx))
      t2 = max(which(!isnanx))
      
      x1<-spline(x[t1:t2],xout = 1:(t2-t1+1))
      xx<-x1$y
      x[t1:t2]<-x1$y
      isnanx<-is.na(x)
      x[isnanx] <- median(x,na.rm = T)
  
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }
    
  }else if(options$method == 3){ # only remove rows with leading and closing zeros
    
    rem1 <- (rowSums(indNaN)==N)
    nanLead <- which(rem1)
    # nanEnd <- which(rem1[length(rem1):1])
    # nanLE <- c(nanEnd,nanLead)
    nanLE<-nanLead
    X[nanLE,] <- NA
    indNaN=is.na(X)
    
  }else if(options$method == 4){ # remove rows with leading and closing zeros & replace missing values
    
    rem1 <- (rowSums(indNaN)==N)
    nanLead <- which(rem1)
    # nanEnd <- which(rem1[length(rem1):1])
    # nanLE <- c(nanEnd,nanLead)
    nanLE<-nanLead
    X[nanLE,] <- NA
    indNaN=is.na(X)
    for (i in 1:N){  
      x = X[,i]
      isnanx = is.na(x)
      t1 = min(which(!isnanx))
      t2 = max(which(!isnanx))
      
      x1<-spline(x[t1:t2],xout = 1:(t2-t1+1))
      xx<-x1$y
      x[t1:t2]<-x1$y
      isnanx<-is.na(x)
      x[isnanx] <- median(x,na.rm = T)
      
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }

  }else if(options$method == 5){
    indNaN=is.na(X)
    
    for (i in 1:N){  
      x = X[,i]
      isnanx = is.na(x)
      t1 = min(which(!isnanx))
      t2 = max(which(!isnanx))
      
      x1<-spline(x[t1:t2],xout = 1:(t2-t1+1))
      xx<-x1$y
      x[t1:t2]<-x1$y
      isnanx<-is.na(x)
      x[isnanx] <- median(x,na.rm = T)
      
      x_MA<-filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
      x_MA=x_MA[(2*k+1):length(x_MA)]
      x[indNaN[,i]]=x_MA[indNaN[,i]]
      X[,i]=x;
    }
  }

  return(list(X = X,indNaN=indNaN))
  
}
