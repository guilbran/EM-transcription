
rm(list=ls())

source('remNaNs_spline.r')
source('InitCond.r')
source('EMstep.R')
source('FIS.r')
source('MissData.r')
source('runKF.r')
source('SKF.r')
source('EM_DFM_SS_block_idioQARMA_restrMQ.r')


# Dados e parametros ------------------------------------------------------

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


# Defino os par?metros

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


# function ----------------------------------------------------------------



EM_DFM_SS_block_idioQARMA_restrMQ(X,Par)


