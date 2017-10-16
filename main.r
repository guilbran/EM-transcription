
rm(list=ls())

source('remNaNs_spline.r')
source('InitCond.r')
source('EMstep.R')
source('FIS.r')
source('MissData.r')
source('runKF.r')
source('SKF.r')
source('EM_DFM_SS_block_idioQARMA_restrMQ.r')
source('em_converged.r')


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

Res<-EM_DFM_SS_block_idioQARMA_restrMQ(X,Par)


# testes ------------------------------------------------------------------

X1<-as.data.frame(BRGDP)
adc<-X1[1:5,]
adc[1:5,]<-NA
X1<-rbind(X1,adc)
X2<-Bpanel(X1,rep(4,dim(X1)[2]))
X2<-cbind(X2,c(rep(NA,27),diff(diff(X1[,8],12),3)))

P<-Par
P$i_idio<-c(rep(T,dim(X2)[2]-1),F)
P$blocks<-matrix(rep(1,dim(X2)[2]*3),dim(X2)[2],3)

Res<-EM_DFM_SS_block_idioQARMA_restrMQ(X2,P)


cbind(X2[,7],Res$X_sm[,7])

ts.plot(Res$X_sm[,7])

ts.plot(cbind(Res$FF[,1],Res$FF[,6],Res$FF[,11]),col=1:3)

i<-26
plot(X[,i],main=i)
lines(Res$X_sm[,i],col=2)
lines(Res$FF[,41],col=4)
ts.plot(Res$FF[,41])
cor(Res$X_sm[,26],Res$FF[,41])

ts.plot(Res$FF[,1:10],col=1:10)
ts.plot(Res$FF[,11:15],col=11:15)
ts.plot(Res$FF[,16:17],col=16:17)

for(i in 1:45){
  Sys.sleep(1)
  ts.plot(Res$FF[,i],main=i)
}

ts.plot(Res$FF[,c(1,6,11)],col=1:3)


1:5  
6:10
11:15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
c(42,43,44,45)

Res$X_sm[,1]

i<-24
ts.plot(cbind(Res$X_sm[,i],X[,i]),col=2:1)


