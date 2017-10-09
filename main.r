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
  
  
  
}



# InitCond ----------------------------------------------------------------

InitCond<-function(x,r,p,blocks,optNaN,Rcon,q,nQ,i_idio){
  
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
  
}









