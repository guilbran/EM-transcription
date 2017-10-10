library(matlab)
library(R.matlab)
library(Matrix)

y <- as.matrix(read.csv("arquivos pra fç EMstep/y.csv", header = F))
A <- as.matrix(read.csv("arquivos pra fç EMstep/A.csv", header = F))
C <- as.matrix(read.csv("arquivos pra fç EMstep/C.csv", header = F))
Q <- as.matrix(read.csv("arquivos pra fç EMstep/Q.csv", header = F))
R <- as.matrix(read.csv("arquivos pra fç EMstep/R.csv", header = F))
Z_0 <- as.matrix(read.csv("arquivos pra fç EMstep/Z_0.csv", header = F))
V_0 <- as.matrix(read.csv("arquivos pra fç EMstep/V_0.csv", header = F))
r <- as.matrix(read.csv("arquivos pra fç EMstep/r0.csv", header = F))
p <- as.matrix(read.csv("arquivos pra fç EMstep/p.csv", header = F))
R_mat <- as.matrix(read.csv("arquivos pra fç EMstep/R_mat.csv", header = F))
q <- as.matrix(read.csv("arquivos pra fç EMstep/q0.csv", header = F))
nQ <- as.matrix(read.csv("arquivos pra fç EMstep/nQ.csv", header = F))
i_idio <- as.matrix(read.csv("arquivos pra fç EMstep/i_idio.csv", header = F))
blocks <- as.matrix(read.csv("arquivos pra fç EMstep/blocks.csv", header = F))
Zsmooth <- as.matrix(read.csv("arquivos pra fç EMstep/Zsmooth.csv", header = F))
Vsmooth <- readMat("arquivos pra fç EMstep/Vsmooth.mat")$Vsmooth
VVsmooth <- readMat("arquivos pra fç EMstep/VVsmooth.mat")$VVsmooth
loglik <- as.matrix(read.csv("arquivos pra fç EMstep/loglik.csv", header = F))


EMstep <- function(y = NULL, A = NULL, C = NULL, Q = NULL, R = NULL, Z_0 = NULL, V_0 = NULL, 
                   r = NULL, p = NULL, R_mat = NULL , q  = NULL, nQ  = NULL, i_idio  = NULL, blocks  = NULL){
  
  n <- size(y,1)
  TT <- size(y,2)
  nM <- n - nQ
  pC <- size(R_mat,2)
  ppC <- max(p,pC)
  n_b <- size(blocks,2)
  
  # Compute the (expected) sufficient statistics for a single Kalman filter sequence.
  
  #Running the Kalman filter with the current estimates of the parameters
  ################ [Zsmooth, Vsmooth, VVsmooth, loglik] = runKF(y, A, C, Q, R, Z_0, V_0);
  
  A_new <- A
  Q_new <- Q
  V_0_new <- V_0
  
  for(i in 1:n_b){
    r_i <- r[i]
    rp <- r_i*p
    if(i == 1){
      rp1 <- 0*ppC
    }else{
      rp1 <- sum(r[1:(i-1)])*ppC
    }
    A_i <- A[(rp1+1):(rp1+r_i*ppC), (rp1+1):(rp1+r_i*ppC)]
    Q_i <- Q[(rp1+1):(rp1+r_i*ppC), (rp1+1):(rp1+r_i*ppC)]
    
    EZZ <- t(Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)] + sum(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),2:dim(Vsmooth)[3]])  # E(Z'Z)
    EZZ_BB <- t(Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)]) %*% Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)] + sum(Vsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),1:(dim(Vsmooth)[3]-1)]) #E(Z(-1)'Z_(-1))
    EZZ_FB <- t(Zsmooth[(rp1+1):(rp1+rp),2:ncol(Zsmooth)]) %*% Zsmooth[(rp1+1):(rp1+rp),1:(ncol(Zsmooth)-1)] + sum(VVsmooth[(rp1+1):(rp1+rp),(rp1+1):(rp1+rp),]) #E(Z'Z_(-1))
    
    A_i[1:r_i,1:rp] <- EZZ_FB[1:r_i,1:rp] %*% solve(EZZ_BB[1:rp,1:rp])
    Q_i[1:r_i,1:r_i] <- (EZZ[1:r_i,1:r_i] - A_i[1:r_i,1:rp] %*% t(EZZ_FB[1:r_i,1:rp])) / TT
    
    A_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] <- A_i 
    Q_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] <- Q_i;
    V_0_new[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC)] <- Vsmooth[(rp1+1):(rp1+r_i*ppC),(rp1+1):(rp1+r_i*ppC),1]
  }
  
  rp1 <- rowSums(r)*ppC
  niM <- sum(i_idio[1:nM])
  
  # idiosyncratic
  EZZ <- diag(diag(Zsmooth[(rp1+1):nrow(Zsmooth),2:ncol(Zsmooth)] %*% t(Zsmooth[(rp1+1):nrow(Zsmooth),2:ncol(Zsmooth)]))) + diag(diag(apply(Vsmooth[(rp1+1):dim(Vsmooth)[1],(rp1+1):dim(Vsmooth)[2],2:dim(Vsmooth)[3]], MARGIN = 1:2, FUN = sum))) #E(Z'Z)
  EZZ_BB <- diag(diag(Zsmooth[(rp1+1):nrow(Zsmooth),1:(ncol(Zsmooth)-1)] %*% t(Zsmooth[(rp1+1):nrow(Zsmooth),1:(ncol(Zsmooth)-1)]))) + diag(diag(apply(Vsmooth[(rp1+1):dim(Vsmooth)[1],(rp1+1):dim(Vsmooth)[2],1:(dim(Vsmooth)[3]-1)], MARGIN = 1:2, FUN = sum)))  #E(Z(-1)'Z_(-1))
  EZZ_FB <- diag(diag(Zsmooth[(rp1+1):nrow(Zsmooth),2:ncol(Zsmooth)] %*% t(Zsmooth[(rp1+1):nrow(Zsmooth),1:(ncol(Zsmooth)-1)]))) + diag(diag(apply(VVsmooth[(rp1+1):dim(VVsmooth)[1],(rp1+1):dim(VVsmooth)[2],], MARGIN = 1:2, FUN = sum))) #E(Z'Z_(-1)) 
  
  A_i <- EZZ_FB %*% diag(1/diag(EZZ_BB))
  Q_i <- (EZZ - A_i %*% t(EZZ_FB)) / TT
  
  A_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] <- A_i[1:niM,1:niM] 
  Q_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] <- Q_i[1:niM,1:niM]
  V_0_new[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM)] <- diag(diag(Vsmooth[(rp1+1):(rp1+niM),(rp1+1):(rp1+niM),1]))
  
  Z_0 <- Zsmooth[,1] #zeros(size(Zsmooth,1),1); #
  
  nanY <- is.nan(y)
  y[nanY] <- 0
  
  # LOADINGS
  C_new <- C
  
  # Blocks
  bl <- unique(blocks)
  n_bl <- size(bl,1)
  bl_idxM <- NULL
  bl_idxQ <- NULL
  R_con <- NULL
  q_con <- NULL
  
  for(i in 1:n_b){
    bl_idxQ <- cbind(bl_idxQ, repmat(bl[,i],1,r[i]*ppC))
    bl_idxM <- cbind(bl_idxM, repmat(bl[,i],1,r[i]), zeros(n_bl,r[i]*(ppC-1)))
    if(i == 1){
      R_con <- kronecker(R_mat,eye(r[i]))
    }else{
      R_con <- as.matrix(bdiag(R_con, kronecker(R_mat,eye(r[i]))))
    }
    q_con <- rbind(q_con, zeros(r[i]*size(R_mat,1),1))
  }
  
  bl_idxM <- bl_idxM == 1
  bl_idxQ <- bl_idxQ == 1
  
  #idio
  i_idio_M <- i_idio[1:nM]
  n_idio_M <- length(find(i_idio_M))
  c_i_idio <- cumsum(i_idio)
  
  for(i in 1:n_bl){
    bl_i <- bl[i,]
    rs <- sum(r[bl_i == 1])
    
    idx_i <- NULL
    for(k in 1:nrow(blocks)){
      idx_i[k] <- sum(blocks[k,] == bl_i) == 3
    }
    idx_i <- find(idx_i)
    
    # MONTHLY
    idx_iM <- idx_i[idx_i < (c(nM) + 1)]
    n_i <- length(idx_iM)
    
    denom <- zeros(n_i*rs,n_i*rs)
    nom <- zeros(n_i,rs)
    
    i_idio_i <- i_idio_M[idx_iM] == 1
    i_idio_ii <- c_i_idio[idx_iM]
    i_idio_ii <- i_idio_ii[i_idio_i]
    
    for(t in 1:TT){
      nanYt <- diag(!nanY[idx_iM,t])
      nn <- sum(bl_idxM[i,])
      denom <- denom + kronecker(Zsmooth[bl_idxM[i,],t+1][1:nn] %*% t(Zsmooth[bl_idxM[i,],t+1][1:nn]) + Vsmooth[bl_idxM[i,],bl_idxM[i,],t+1][1:nn,1:nn], nanYt)
      nom <- nom + y[idx_iM,t] %*% t(Zsmooth[bl_idxM[i,],t+1][1:nn]) - nanYt[,i_idio_i] %*% (Zsmooth[rp1+i_idio_ii,t+1] %*% t(Zsmooth[bl_idxM[i,],t+1][1:nn]) + Vsmooth[rp1+i_idio_ii,bl_idxM[i,],t+1][,1:nn])
    }
    
    vec_C <- solve(denom) %*% c(nom)
    C_new[idx_iM,bl_idxM[i,]][,1:nn] <- reshape(vec_C,n_i,rs)
    
    # QUARTERLY
    idx_iQ <- idx_i[idx_i>c(nM)]
    rps <- rs*ppC
    
    R_con_i <- R_con[,bl_idxQ[i,]]
    q_con_i <- q_con
    no_c <- !(rowSums(R_con_i == 0) == ncol(R_con_i))
    R_con_i <- R_con_i[no_c,]
    q_con_i <- q_con_i[no_c,] 
    
    if(i != 1){
      for(j in idx_iQ){
        denom <- zeros(rps,rps)
        nom <- zeros(1,rps)
        idx_jQ <- j-c(nM)
        i_idio_jQ <- (rp1+n_idio_M+5*(idx_jQ-1)+1):(rp1+n_idio_M+5*idx_jQ)
        V_0_new[i_idio_jQ,i_idio_jQ] <- Vsmooth[i_idio_jQ,i_idio_jQ,1]
        A_new[i_idio_jQ[1],i_idio_jQ[1]] <- A_i[i_idio_jQ[1]-rp1,i_idio_jQ[1]-rp1]
        Q_new[i_idio_jQ[1],i_idio_jQ[1]] <- Q_i[i_idio_jQ[1]-rp1,i_idio_jQ[1]-rp1]
        
        for(t in 1:TT){
          nanYt <- as.vector(!nanY[j,t])*1
          nn2 <- sum(bl_idxQ[i,])
          denom <- denom + kronecker(Zsmooth[bl_idxQ[i,],t+1][1:nn2] %*% t(Zsmooth[bl_idxQ[i,],t+1][1:nn2]) + Vsmooth[bl_idxQ[i,],bl_idxQ[i,],t+1][1:nn2,1:nn2],nanYt)
          nom <- nom + y[j,t] %*% t(Zsmooth[bl_idxQ[i,],t+1][1:nn2])
          nom <- nom - nanYt %*% (matrix(c(1,2,3,2,1), nrow = 1) %*% Zsmooth[i_idio_jQ,t+1] %*% t(Zsmooth[bl_idxQ[i,],t+1][1:nn2]) +
                                    matrix(c(1,2,3,2,1), nrow = 1) %*% Vsmooth[i_idio_jQ,bl_idxQ[i,],t+1][,1:nn2])
        }
        C_i <- solve(denom) %*% t(nom)
        C_i_constr <- C_i - solve(denom) %*% t(R_con_i) %*% solve(R_con_i %*% solve(denom) %*% t(R_con_i)) %*% (R_con_i %*% C_i - q_con_i)
        C_new[j,bl_idxQ[i,]] <- C_i_constr
      }
    }
  }
  
  R_new <- zeros(n,n)
  for(t in 1:TT){
    nanYt <- diag(!nanY[,t])
    R_new <- R_new + (y[,t] - nanYt %*% C_new %*% Zsmooth[,t+1]) %*% t(y[,t] - nanYt %*% C_new %*% Zsmooth[,t+1]) + nanYt %*% C_new %*% Vsmooth[,,t+1] %*% t(C_new) %*% nanYt + (eye(n)-nanYt) %*% R %*% (eye(n)-nanYt)
  }
  
  R_new <- R_new/TT
  RR <- diag(R_new) #RR(RR<1e-2) = 1e-2;
  RR[i_idio_M] <- 1e-04
  RR[(nM+1):length(RR)] <- 1e-04
  R_new <- diag(RR)
  
  # output
  k <- list(C_new = C_new, R_new = R_new, A_new = A_new, Q_new = Q_new, Z_0 = Z_0, V_0 = V_0, loglik = loglik)
  
}
