setwd('C:/Users/b32750/Downloads/EM/R')

source('remNaNs_spline.r')

# Dados e parâmetros ------------------------------------------------------

# importo a base de dados
library(xlsx)
dados<-read.xlsx2('dados.xlsx',sheetIndex = 1,
                  header = F,colIndex = 1:26,
                  colClasses = rep('numeric',26))
nas<-data.frame(data.frame(matrix(NA,13,26)))
X<-rbind(dados,nas)

for(i in 1:26){
X[is.nan(X[,i]),i]<-rep(NA,sum(is.nan(X[,i])))
}


# Defino os parâmetros

Par<-list(r=c(1,1,1),p=1,max_iter=500,i_idio=c(rep(T,25),F),
          Rconstr = matrix(c(
            c(2,3,2,1),
            c(-1,0,0,0),
            c(0,-1,0,0),
            c(0,0,-1,0),
            c(0,0,0,-1))
            ,4,5),
          q = matrix(rep(0,4),4,1),nQ = 1,
          blocks = matrix(c(1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    0,
                    0,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    1,
                    0,
                    1,
                    1,
                    0,
                    1,
                    1,
                    1,
                    0,
                    0,
                    0,
                    0,
                    1,
                    0,
                    0,
                    1,
                    1,
                    1,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    1,
                    0,
                    0,
                    1,
                    0,
                    0,
                    0,
                    1,
                    1,
                    1,
                    1,
                    0,
                    1,
                    1,
                    0),26,3))



# EM_DFM_SS_block_idioQARMA_restrMQ ---------------------------------------

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
  
  ### Preparação dos dados
  
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
  
  Resinitcond<-InitCond(xNaN,r,p,blocks,optNaN,R_mat,q,nQ,i_idio);
  
  A<-Resinitcond$A
  C<-Resinitcond$C
  Q<-Resinitcond$Q
  R<-Resinitcond$R
  initZ<-Resinitcond$initZ
  initV<-Resinitcond$initV
  
  
  
  
}



# InitCond ----------------------------------------------------------------

InitCond<-function(x,r,p,blocks,optNaN,Rcon,q,nQ,i_idio){
  
  library(magic)
  
  pC = size(Rcon,2)
  ppC = max(p,pC)
  n_b = size(blocks,2)
  
  OPTS<-list()
  OPTS$disp=0
  
  res_remNaNs_spline <- remNaNs_spline(x,optNaN)
  xBal<- res_remNaNs_spline$X
  indNaN <- res_remNaNs_spline$indNaN
  TT <- dim(xBal)[1]
  N <- dim(xBal)[2]
  NM <- N-nQ
  
  xNaN = xBal
  
  for(i in 1:N){
    xNaN[indNaN[,i],i] <- NA
  }  
  
  
  C = {}
  A = {}
  Q = {}
  initV = {}
  
  res = xBal
  resNaN = xNaN
  indNaN[1:pC-1,] <- T
  
  
  for(i in 1:n_b){    # roda o loop em cada um dos blocos (geral, real e nominal)
    r_i<-r[i]
    
    ########################
    # Observation equation #
    ########################
    
    C_i = zeros(N,r_i*ppC)
    idx_i = find(blocks[,i])
    
    # Atenção aqui funciona pois a(s) última(S) variável(is) da base de dados é(são) trimestral(is)!
    idx_iM = idx_i[idx_i<NM+1];   # índice representando variável mesal
    idx_iQ = idx_i[idx_i>NM];     # índice representando variável trimestral
    
    eig<-eigen(cov(res[,idx_iM]))
    v<-eig$vectors[,1:r_i]
    d<-eig$values[1:r_i]
    
    C_i[idx_iM,1:r_i] = v
    f = as.matrix(res[,idx_iM])%*%as.matrix(v)
    for(kk in 0:(max(p+1,pC)-1)){
      if(kk == 0){
        FF<-f[(pC-kk):(dim(f)[1]-kk),]
      }else{
        FF <- cbind(FF,f[(pC-kk):(dim(f)[1]-kk),])
      }
    }
    
    Rcon_i = kronecker(Rcon,eye(r_i))
    q_i = kronecker(q,zeros(r_i,1));
    ff = FF[,1:(r_i*pC)]
    
    for(j in idx_iQ){     # Coeficiente "loadings" de Variáveis trimestrais
      xx_j = resNaN[pC:dim(resNaN)[1],j]
      if(sum(!is.na(xx_j)) < size(ff,2)+2){
        xx_j = res[pC:dim(res)[1],j]
      }
      ff_j = ff[!is.na(xx_j),]
      xx_j = xx_j[!is.na(xx_j)]
      iff_j = solve(t(ff_j)%*%ff_j)
      Cc = iff_j%*%t(ff_j)%*%xx_j
      Cc = Cc - iff_j%*%t(Rcon_i)%*%solve(Rcon_i%*%iff_j%*%t(Rcon_i))%*%(Rcon_i%*%Cc-q_i);
      C_i[j,1:(pC*r_i)] <- t(Cc)
    }
    
    ff = rbind(zeros(pC-1,pC*r_i),ff)
    res = res - ff%*%t(C_i)
    resNaN = res
    for(i_aux in 1:dim(indNaN)[2]){
      resNaN[indNaN[,i_aux],i_aux] <- NA
    }
    C <- cbind(C,C_i)
    
    #######################    
    # Transition Equation #
    #######################  
    z <- FF[,1:r_i]
    Z <- FF[,(r_i+1):(r_i*(p+1))]
    A_i = t(zeros(r_i*ppC,r_i*ppC))
    A_temp = solve(t(Z)%*%Z)%*%t(Z)%*%z
    A_i[1:r_i,1:r_i*p] <- t(A_temp)
    A_i[(r_i+1):dim(A_i)[1],1:(r_i*(ppC-1))] <- eye(r_i*(ppC-1))
    
    ##########################
    
    Q_i = zeros(ppC*r_i,ppC*r_i)
    e = z  - Z%*%A_temp         # VAR residuals
    Q_i[1:r_i,1:r_i] = cov(e);  # VAR covariance matrix
    
    initV_i = reshape(solve(eye((r_i*ppC)^2)-kronecker(A_i,A_i))%*%c(Q_i),r_i*ppC,r_i*ppC);
    
    if(is.null(A)){
      A<-A_i
    }else{
      A <- adiag(A,A_i)  
    }
    
    if(is.null(Q)){
      Q<-Q_i
    }else{
      Q <- adiag(Q,Q_i)  
    }
    
    if(is.null(initV)){
      initV<-initV_i
    }else{
      initV <- adiag(initV,initV_i)  
    }
    
    # linha 401 do código em matlab
    
  }
  
  
  R = diag(diag(var(resNaN,na.rm = T)))
  
  eyeN = eye(N)
  eyeN<-eyeN[,i_idio]
  
  # Initial conditions
  C=cbind(C,eyeN)
  
  ii_idio = find(i_idio)
  n_idio = length(ii_idio)
  B = zeros(n_idio)
  S = zeros(n_idio)
  
  BM<-zeros(n_idio)
  SM<-zeros(n_idio)
  
  for (i in 1:n_idio){
    R[ii_idio[i],ii_idio[i]] <- 1e-04
    
    res_i = resNaN[,ii_idio[i]]
    # number of leading zeros
    leadZero = max( find( t(1:TT) == cumsum(is.na(res_i)) ) )
    endZero = max( find( t(1:TT) == cumsum(is.na(res_i[length(res_i):1])) ) );
    
    res_i<-res_i[(leadZero+1):(length(res_i)-endZero)]
    
    BM[i,i] = solve(t(res_i[1:(length(res_i)-1)])%*%res_i[1:(length(res_i)-1)])%*%t(res_i[1:(length(res_i)-1)])%*%res_i[2:length(res_i)] 
    SM[i,i] = var(res_i[2:length(res_i)]-res_i[1:(length(res_i)-1)]*BM[i,i])
    # SM[i,i] = var(res_i[2:length(res_i)]-res_i[1:(length(res_i)-1)]*B[i,i])
    # ATENÇÃO: Aqui os autores usam B[i,i], porém esse valor é 0. Então eu uso BM[i,i]
  }
  
  initViM = diag(1/diag(eye(size(BM,1))-BM^2))%*%SM;
  
  
  C<-cbind(C,rbind(zeros(NM,5*nQ),t(kronecker(eye(nQ),c(1,2,3,2,1)))))
  Rdiag<-diag(R)
  sig_e <- Rdiag[(NM+1):N]/19
  Rdiag[(NM+1):N] <- 1e-04
  R = diag(Rdiag)
  
  rho0<-0.1
  
  BQ <- kronecker(eye(nQ),rbind(cbind(rho0,zeros(1,4)),cbind(eye(4),zeros(4,1))))
  temp = zeros(5)
  temp[1,1] = 1
  if(is.matrix(sig_e)){
    SQ = kronecker(diag((1-rho0^2)*sig_e),temp)
  }else{
    SQ = kronecker((1-rho0^2)*sig_e,temp)
  }
  
  initViQ = reshape(solve(eye((5*nQ)^2)-kronecker(BQ,BQ))%*%c(SQ),5*nQ,5*nQ)
  
  # BQ = kronecker(eye(nQ),rbind(zeros(1,5),cbind(eye(4),zeros(4,1))))
  # temp = zeros(5)
  # temp[1,1] = 1
  # if(is.matrix(sig_e)){
  #   SQ = kronecker(diag(sig_e),temp)
  # }else{
  #   SQ = kronecker(diag(as.matrix(sig_e)),temp)  
  # }
  # temp = matrix(c(19, 16, 10, 4, 1, 16, 19, 16, 10, 4, 10, 16, 19, 16, 10, 4, 10, 16, 19, 16, 1, 4, 10, 16, 19),
  #               5,5)
  # if(is.matrix(sig_e)){
  #   initViQ = kronecker(diag(sig_e),temp);
  # }else{
  #   initViQ = kronecker(sig_e,temp);
  # }
  
  A1<-magic::adiag(A,BM,BQ)
  Q1<-magic::adiag(Q, SM, SQ)
  
  A<-A1
  Q<-Q1
  
  # Initial conditions
  initZ = zeros(size(A,1),1); ##[randn(1,r*(nlag+1))]';
  initV = magic::adiag(initV, initViM, initViQ)
  
  return(list(A = A, C = C, Q = Q, R = R, initZ = initZ, initV = initV))
  
}
