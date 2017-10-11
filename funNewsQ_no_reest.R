library(R.matlab)
library(matlab)
library(readxl)
library(zoo)
P <- readMat("arquivos pra fç EMstep/P.mat")
P <- P$P[,,1]
P$DataFile <- "C:\\Users\\daiane.mattos\\Desktop\\NowCastingOxfordReplicationWeb\\Data\\Data 2010-01-21 Back oil15.xls"
P$BlockFile <- "C:\\Users\\daiane.mattos\\Desktop\\NowCastingOxfordReplicationWeb\\Data\\RNBen.xls"

funNewsQ_no_reest <- function(P){
  
  
  StartEst <- P$StartEst
  Qnews <- P$Qnews
  SerNews <- P$SerNews
  StartEv <- P$StartEv
  EndEv <- P$EndEv
  P.max_iter <- 500
  fcstH <- 1
  
  # %--------------------------------------------------------------------------
  # % Loading monthly data
  # %--------------------------------------------------------------------------
  
  DataFile <-  P$DataFile
  BlockFile <- P$BlockFile
  
  a <- data.frame(read_excel(BlockFile, sheet = 1))
  ListM <- a[,4]
  ListM <- find(ListM)
  Blocks <- a[,4:ncol(a)]
  Blocks <- Blocks[ListM,]
  
  aa <- data.frame(read_excel(DataFile, sheet = 3))
  a <- aa
  a[,2:3] <- NaN
  b <- aa[,-1]
  b[,-c(1:2)] <- NA
  b <- b[2:nrow(b),]
  GroupM <- b[ListM,2]
  SeriesM <- b[ListM,3]
  # Transformation
  TransfM <- a[ListM,4:5]
  # unbalancedeness patterns
  UnbM <- a[ListM,6:11] 
  
  a <- data.frame(read_excel(DataFile, sheet = 1, skip = 3, col_names = F)[,-1])
  b <- data.frame(read_excel(DataFile, sheet = 1, col_names = T))
  b[3:nrow(b),2:ncol(b)] <- NA
  
  DataM <- a[,ListM]
  # if strcmp(version('-release'),'2006b')
  
  DatesM <- as.Date(data.frame(b[3:nrow(b),1])[,1])
  # else
  #     DatesM <- datenum(b(4:end,1));
  # end
  DatesMV <- data.frame(ano = as.numeric(substr(DatesM,1,4)),
                        mes = as.numeric(substr(DatesM,6,7)))
  
  TT <- length(DatesM)
  
  # MoM transformations
  DataMM <- DataM
  DataMM[,c(TransfM[,1] == 1)] <- 100*log(DataMM[,c(TransfM[,1] == 1)])
  DataMM[2:nrow(DataMM),c(TransfM[,2] == 1)] <- DataMM[2:nrow(DataMM),c(TransfM[,2] == 1)] - DataMM[1:(nrow(DataMM)-1),c(TransfM[,2] == 1)]
  DataMM[1,c(TransfM[,2] == 1)] <- NaN
  
  GroupSurveys <- c('ECSurv','ECSurvNom','PMI','PMInom')
  
  if(P$SL == 1){
    DataMM[,GroupM %in% GroupSurveys] <- DataM[,GroupM %in% GroupSurveys];
  }
  
  
  DataMTrf <- DataMM
  
  tM <- nrow(DataMTrf)
  nM <- ncol(DataMTrf)
  
  x <- matrix(NaN, ncol = nM, nrow = TT-tM)
  colnames(x) <- colnames(DataMTrf)
  DataMTrf <- rbind(DataMTrf,x)
  x <- matrix(NaN, ncol = nM, nrow = TT-tM)
  colnames(x) <- colnames(DataM)
  DataMM <- rbind(DataMM, x)
  
  # %--------------------------------------------------------------------------
  # % Loading quarterly data
  # %--------------------------------------------------------------------------
  a <- data.frame(read_excel(BlockFile, sheet = 2, col_names = T))
  ListQ <- a[,4]
  ListQ <- find(ListQ)
  BlocksQ <- a[,4:ncol(a)]
  BlocksQ <- BlocksQ[ListQ,]
  
  aa <- data.frame(read_excel(DataFile, sheet = 4))
  a <- aa
  a[,2:3] <- NaN
  b <- aa[,-1]
  b[,-c(1:2)] <- NA
  b <- b[2:nrow(b),]
  GroupQ <- b[ListQ,2]
  SeriesQ <- b[ListQ,3]
  # Transformation
  Transf <- a[ListQ,4:5]
  # unbalancedeness patterns
  UnbQ <- a[ListQ,6:11]
  
  a <- data.frame(read_excel(DataFile, sheet = 2, skip = 3, col_names = F)[,-1])
  b <- data.frame(read_excel(DataFile, sheet = 2, col_names = T))
  b[3:nrow(b),2:ncol(b)] <- NA
  DataQ <- a[,ListQ]
  DataQTrf <- data.frame(DataQ)
  DataQTrf[,Transf[,1] == 1] <- log(DataQTrf[,Transf[,1] == 1])
  DataQTrf[2:nrow(DataQTrf),Transf[,2] == 1] <- 100*(DataQTrf[2:nrow(DataQTrf),Transf[,2] == 1] - DataQTrf[1:(nrow(DataQTrf)-1),Transf[,2] == 1])
  DataQTrf[1,Transf[,2] == 1] <- NaN
  
  # quarterly at monthly frequency
  DataQMTrf <- kronecker(as.matrix(DataQTrf),c(NaN,NaN,1))
  
  tQ <- nrow(DataQMTrf)
  nQ <- ncol(DataQMTrf)
  
  x <- matrix(NaN, ncol = nQ, nrow = TT-tQ)
  colnames(x) <- colnames(DataQMTrf)
  DataQMTrf <- rbind(DataQMTrf,x)
  
  # %--------------------------------------------------------------------------
  # % complete dataset
  # %--------------------------------------------------------------------------
  
  Data <- cbind(DataMTrf,DataQMTrf)
  Series <- rbind(SeriesM,SeriesQ)
  Group <- rbind(GroupM,GroupQ)
  UnbPatt <- rbind(UnbM,UnbQ)
  
  P$blocks <- rbind(Blocks,BlocksQ)
  
  iEst <- find(DatesMV[,1] == StartEst[1] & DatesMV[,2] == StartEst[2])
  
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Data <- Data[iEst:nrow(Data),]     
  Dates <- DatesM[iEst:length(DatesM)]     
  DatesV <- DatesMV[iEst:nrow(DatesMV),]   
  idxM <- t(1:nM)                  
  idxQ <- t((nM+1):(nM+nQ))
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  DataMM <- DataMM[iEst:nrow(DataMM),]
  nVar <- nM+nQ
  
  # %--------------------------------------------------------------------------
  # % unbalancedness patterns
  # %--------------------------------------------------------------------------
  nn <- min(UnbPatt)
  nn <- min(nn,0)
  
  UnbPattM1_1 <- zeros(12-nn,nVar)
  UnbPattM1_2 <- zeros(12-nn,nVar)
  UnbPattM2_1 <- zeros(12-nn,nVar)
  UnbPattM2_2 <- zeros(12-nn,nVar)
  UnbPattM3_1 <- zeros(12-nn,nVar)
  UnbPattM3_2 <- zeros(12-nn,nVar)
  
  nUnb <- 12-nn
  
  for(i in 1:(nVar-2)){
    UnbPattM1_1[(nrow(UnbPattM1_1)-UnbPatt[i,1]+1+nn):nrow(UnbPattM1_1),i] <- NaN
    UnbPattM2_1[(nrow(UnbPattM2_1)-UnbPatt[i,2]+1+nn):nrow(UnbPattM2_1),i] <- NaN
    UnbPattM3_1[(nrow(UnbPattM3_1)-UnbPatt[i,3]+1+nn):nrow(UnbPattM3_1),i] <- NaN
    if(i < (nUnb-1)){
      UnbPattM1_2[(nrow(UnbPattM1_2)-UnbPatt[i,4]+1+nn):nrow(UnbPattM1_2),i] <- NaN
      UnbPattM2_2[(nrow(UnbPattM2_2)-UnbPatt[i,5]+1+nn):nrow(UnbPattM2_2),i] <- NaN
      UnbPattM3_2[(nrow(UnbPattM3_2)-UnbPatt[i,6]+1+nn):nrow(UnbPattM3_2),i] <- NaN
    }
  }
  
  # %--------------------------------------------------------------------------
  # % restrictions
  # %--------------------------------------------------------------------------
  P$nQ <- nQ
  P$Rconstr <- matrix(c(2, -1, 0, 0, 0,
                        3, 0, -1, 0, 0,
                        2, 0, 0, -1, 0,
                        1, 0, 0, 0, -1), byrow = T, ncol = 5, nrow = 4)
  P$q <- zeros(4,1)
  P$restr <- '_restrMQ'
  
  # %--------------------------------------------------------------------------
  # % out-of-sample evaluation
  # %--------------------------------------------------------------------------
  
  iS <- find(DatesV[,1] == StartEv[1] & DatesV[,2] == StartEv[2])
  iE <- find(DatesV[,1] == EndEv[1] & DatesV[,2] == EndEv[2])
  iQ <- find(DatesV[,1] == Qnews[1] & DatesV[,2] == Qnews[2])
  iSer <- find(Series %in% SerNews)
  
  
  Month <- mod(DatesV[,2],3)
  Month[Month == 0] <- 3
  
  P$i_idio <- rbind(ones(nM,1),zeros(nQ,1)) == 1
  
  Month_i <- Month[iS-1]
  # second unbalancedeness pattern
  
  eval(parse(text = paste0('UnbP = UnbPattM',Month_i,'_2')))
  
  X <- Data[1:(iS-1-nn),]
  temp <- X[(nrow(X)-nUnb+1):nrow(X),]
  temp[is.nan(UnbP)] <- NaN
  X[(nrow(X)-nUnb+1):nrow(X),] <- temp
  x <- matrix(NaN, nrow = max(0,(fcstH+1)*3-Month_i+nn), ncol = nM+nQ)
  colnames(x) <- colnames(X)
  X_old <- rbind(X,x)
  
  OldFcst <- zeros(2*(iE-iS+1),1);
  NewFcst <- zeros(2*(iE-iS+1),1);
  GroupNews <- zeros(2*(iE-iS+1),length(unique(Group)));
  SerNews <- zeros(2*(iE-iS+1),nM+nQ);
  Gain <- zeros(2*(iE-iS+1),nM+nQ);
  
  for(i in iS:iE){
    
    Date_i <- DatesV[i,]
    Month_i <- Month[i];
    message(paste0('Computing the news for the vintages: y', DatesV[i,1],' m', DatesV[i,2]))
    
    # first unbalancedeness pattern
    eval(paste(text = paste0('UnbP = UnbPattM',Month_i,'_1;')))
    X <- Data[1:(i-nn),]
    temp <- X[(nrow(X)-nUnb+1):nrow(X),]
    temp[is.nan(UnbP)] <- NaN
    X[(nrow(X)-nUnb+1):nrow(X),] <- temp
    x <- matrix(NaN, nrow = max(0,(fcstH+1)*3-Month_i+nn), ncol = nM+nQ)
    colnames(x) <- colnames(X)
    X_new <- rbind(X,x)
    
    T_o <- nrow(X_old)
    T_n <- nrow(X_new)
    x <- matrix(NaN, nrow = T_n-T_o, ncol = nM+nQ)
    colnames(x) <- colnames(X_old)
    X_old <- rbind(X_old,x)
    
    if(i == iS){
      eval(parse(text = paste0('R_new <- EM_DFM_SS',P$method,P$idio,P$restr,'(X_new,P)')))
      R_new$Groups <- Group
      R_new$Series <- Series
    }
    
    # ATENÇÃO AQUI BICHAUM
    out <- News_DFM_ML(X_old,X_new,R_new,iQ,iSer)
    
    OldFcst[2*(i-iS)+1,1] <- out$OldFcst
    NewFcst[2*(i-iS)+1,1] <- out$NewFcst
    GroupNews[2*(i-iS)+1,] <- out$GroupNews
    SerNews[2*(i-iS)+1,] <- out$SerNews
    gainT <- out$gainT
    serGainT <- out$serGainT
    Actual[,2*(i-iS)+1] <- out$Actual
    Fcst[,2*(i-iS)+1] <- out$Fcst
    Filt[,2*(i-iS)+1] <- out$Filt  
    Gain[2*(i-iS)+1,Series %in% serGainT] <- gainT
    X_old <- X_new
    
    # second unbalancedeness pattern
    eval(parse(text = paste0('UnbP <- UnbPattM',Month_i,'_2')))
    X <- Data[1:(i-nn),]
    temp <- X[(nrow(X)-nUnb+1):nrow(X),]
    temp[is.nan(UnbP)] <- NaN
    X[(nrow(X)-nUnb+1):nrow(X),] <- temp
    x <- matrix(NaN, nrow = max(0,(fcstH+1)*3-Month_i+nn), ncol = nM+nQ)
    colnames(x) <- colnames(X)
    X_new <- rbind(X,x)
    
    # ATENÇÃO AQUI DE NOVO BICHAUM
    out2 <- News_DFM_ML(X_old,X_new,R_new,iQ,iSer)
    OldFcst[2*(i-iS)+2,1] <- out2$OldFcst
    NewFcst[2*(i-iS)+2,1] <- out2$NewFcst
    GroupNews[2*(i-iS)+2,] <- out2$GroupNews
    SerNews[2*(i-iS)+2,] <- out2$SerNews
    gainT <- out2$gainT
    serGainT <- out2$serGainT
    Actual[,2*(i-iS)+2] <- out2$Actual
    Fcst[,2*(i-iS)+2] <- out2$Fcst
    Filt[,2*(i-iS)+2] <- out2$Filt  
    Gain[2*(i-iS)+2,Series %in% serGainT] <- gainT
    X_old <- X_new
    
  }
  
  DatesNews <- matrix(NaN, nrow = length(OldFcst), ncol = 2)
  DatesNews[seq(1,nrow(DatesNews), by = 2),] <- DatesV[iS:iE,]
  DatesNews[seq(1,nrow(DatesNews), by = 2),] <- DatesV[iS:iE,]
  
  GroupNames <- t(unique(Group))
  TrueSer <- Data[iQ,iSer]
  
  # check whether the new forecats is equal to the old forecast plus the news 
  check <- NewFcst-OldFcst-matrix(rowSums(GroupNews))
  
  datafile <- paste0('news',P$DF,P$method,P$idio,paste0(P$r, collapse =""),P$p)
  # datafile <- strrep(datafile,' ','_');
  
  # output
  list(OldFcst = OldFcst,NewFcst = NewFcst, TrueSer = TrueSer, DatesNews = DatesNews, GroupNews = GroupNews, 
       SerNews = SerNews, GroupNames = GroupNames, Gain = Gain, Fcst = Fcst, Actual = Actual, Filt = Filt,
       Series = Series, Group = Group, P = P)
  
}