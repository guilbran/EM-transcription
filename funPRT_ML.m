%%%  Replication files for:
%%%  ""Nowcasting", 2010, (by Marta Banbura, Domenico Giannone and Lucrezia Reichlin), 
%%% in Michael P. Clements and David F. Hendry, editors, Oxford Handbook on Economic Forecasting.
%%%
%%% The software can be freely used in applications. 
%%% Users are kindly requested to add acknowledgements to published work and 
%%% to cite the above reference in any resulting publications
function funPRT_ML(P)

StartEst = P.StartEst;


P.max_iter = 500;

fcstH_M = P.fcstH_M;
fcstH_Q = P.fcstH_Q;

%--------------------------------------------------------------------------
% Loading monthly data
%--------------------------------------------------------------------------
DataFile =  P.DataFile;
BlockFile = P.BlockFile;

a = xlsread(BlockFile,'BlocksM');
ListM = a(:,4);
ListM = find(ListM);
Blocks = a(:,4:end);
Blocks = Blocks(ListM,:);

[a,b] = xlsread(DataFile,'LegendMonthly');

b = b(2:end,:);
GroupM = b(ListM,2);
SeriesM = b(ListM,3);
%Transformation
TransfM = a(ListM,4:5);
% unbalancedeness patterns
UnbM = a(ListM,6:11); 

[a,b] = xlsread(DataFile,'DataMonthly');

DataM = a(:,ListM);
% if strcmp(version('-release'),'2006b')
    DatesM = datenum(b(4:end,1),'dd/m/yy');
% else
%     DatesM = datenum(b(4:end,1));
% end
DatesMV = datevec(DatesM);
DatesMV = DatesMV(:,1:2);

T = length(DatesM);


% MoM transformations
DataMM = DataM;
DataMM(:,TransfM(:,1) == 1) = log(DataMM(:,TransfM(:,1) == 1));
DataMM(2:end,TransfM(:,2) == 1) = 100*(DataMM(2:end,TransfM(:,2) == 1) ...
    - DataMM(1:end-1,TransfM(:,2) == 1));
DataMM(1,TransfM(:,2) == 1) = nan;

GroupSurveys = {'ECSurv','ECSurvNom','PMI','PMInom'};
if P.SL
    DataMM(:,ismember(GroupM,GroupSurveys)) = DataM(:,ismember(GroupM,GroupSurveys));
end

DataMTrf = DataMM;

[tM,nM] = size(DataMTrf);
DataMTrf = [DataMTrf;nan(T-tM,nM)];
DataMM = [DataMM;nan(T-tM,nM)];
%--------------------------------------------------------------------------
% Loading quarterly data
%--------------------------------------------------------------------------
a = xlsread(BlockFile,'BlocksQ');
ListQ = a(:,4);
ListQ = find(ListQ);
BlocksQ = a(:,4:end);
BlocksQ = BlocksQ(ListQ,:);

[a,b] = xlsread(DataFile,'LegendQuarterly');

b = b(2:end,:);
GroupQ = b(ListQ,2);
SeriesQ = b(ListQ,3);
%Transformation
Transf = a(ListQ,4:5);
% unbalancedeness patterns
UnbQ = a(ListQ,6:11); 

[a,b] = xlsread(DataFile,'DataQuarterly');

DataQ = a(:,ListQ);

DataQTrf = DataQ;
DataQTrf(:,Transf(:,1) == 1) = log(DataQTrf(:,Transf(:,1) == 1));
DataQTrf(2:end,Transf(:,2) == 1) = 100*(DataQTrf(2:end,Transf(:,2) == 1) ...
    - DataQTrf(1:end-1,Transf(:,2) == 1));
DataQTrf(1,Transf(:,2) == 1) = nan;

% quarterly at monthly frequency
DataQMTrf = kron(DataQTrf,[nan;nan;1]);

[tQ,nQ] = size(DataQMTrf);
DataQMTrf = [DataQMTrf;nan(T-tQ,nQ)];
%--------------------------------------------------------------------------
% complete dataset
%--------------------------------------------------------------------------

Data = [DataMTrf DataQMTrf];
Series = [SeriesM;SeriesQ];
Group = [GroupM;GroupQ];
UnbPatt = [UnbM;UnbQ];

P.blocks = [Blocks;BlocksQ];

iEst = find(DatesMV(:,1) == StartEst(1) & DatesMV(:,2) == StartEst(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data = Data(iEst:end,:);         %
Dates = DatesM(iEst:end,:);      %
DatesV = DatesMV(iEst:end,:);     %
idxM = (1:nM)';                  %
idxQ = (nM+1:nM+nQ)';            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DataMM = DataMM(iEst:end,:);

nVar = nM+nQ;

%--------------------------------------------------------------------------
% unbalancedness patterns
%--------------------------------------------------------------------------
nn = min(min(UnbPatt));
nn = min(nn,0);

UnbPattM1_1 = zeros(12-nn,nVar);
UnbPattM1_2 = zeros(12-nn,nVar);
UnbPattM2_1 = zeros(12-nn,nVar);
UnbPattM2_2 = zeros(12-nn,nVar);
UnbPattM3_1 = zeros(12-nn,nVar);
UnbPattM3_2 = zeros(12-nn,nVar);

nUnb = 12-nn;

for i = 1:nVar
    UnbPattM1_1(end-UnbPatt(i,1)+1+nn:end,i) = nan;
    UnbPattM2_1(end-UnbPatt(i,2)+1+nn:end,i) = nan;
    UnbPattM3_1(end-UnbPatt(i,3)+1+nn:end,i) = nan;
    UnbPattM1_2(end-UnbPatt(i,4)+1+nn:end,i) = nan;
    UnbPattM2_2(end-UnbPatt(i,5)+1+nn:end,i) = nan;
    UnbPattM3_2(end-UnbPatt(i,6)+1+nn:end,i) = nan;
end


%--------------------------------------------------------------------------
% restrictions
%--------------------------------------------------------------------------
P.nQ = nQ;
P.Rconstr = [2 -1 0 0 0;...
        3 0 -1 0 0;...
        2 0 0 -1 0;...
        1 0 0 0 -1];
P.q = zeros(4,1);
P.restr = '_restrMQ';

%--------------------------------------------------------------------------
% out-of-sample evaluation
%--------------------------------------------------------------------------

iS = find(DatesV(:,1) == P.StartFcst(1) & DatesV(:,2) == P.StartFcst(2));
iE = find(DatesV(:,1) == P.EndFcst(1) & DatesV(:,2) == P.EndFcst(2));

FcstQ = nan(iE,(fcstH_Q+2)*2,nQ);
FcstM1_1 = nan(iE,(fcstH_M+2),nM);
FcstM1_2 = nan(iE,(fcstH_M+2),nM);
FcstM12_1 = nan(iE,(fcstH_M+2),nM);
FcstM12_2 = nan(iE,(fcstH_M+2),nM);

Month = mod(DatesV(:,2),3);
Month(Month == 0) = 3;

P.i_idio = logical([ones(nM,1);zeros(nQ,1)]);
P.Series = Series;

% for i = iS:iE
    
    i = iE
    
    Date_i = DatesV(i,:);
    Month_i = Month(i);
    disp(['Computing the predictions for the vintages: y', num2str(DatesV(i,1)),' m',num2str(DatesV(i,2))])
    
    
    % first unbalancedeness pattern
    eval(['UnbP = UnbPattM',int2str(Month_i),'_1;']);
    
    X = Data(1:i-nn,:);
    temp = X(end-nUnb+1:end,:);
    temp(isnan(UnbP)) = nan;
    X(end-nUnb+1:end,:) = temp;
    n_nan = max([0,fcstH_M+nn,(fcstH_Q+1)*3-Month_i+nn]);
    X = [X;nan(n_nan,nM+nQ)];
    
    xlswrite('dados.xlsx',X)

    eval(['Res = EM_DFM_SS',P.method,P.idio,P.restr,'(X,P);'])
    
    XM = DataMM(1:i-nn,:);
    temp = XM(end-nUnb+1:end,:);
    temp(isnan(UnbP(:,1:nM))) = nan;
    XM(end-nUnb+1:end,:) = temp;
    
    tempFcstM = cnvrt2Mgr(Res.X_sm(:,1:nM),X(:,1:nM),XM,TransfM(:,2),P.transf);
    FcstM1_1(i,:,:) = tempFcstM([i-1:i+fcstH_M],:);
    tfm = filter(ones(12,1),1,tempFcstM);
    FcstM12_1(i,:,:) = tfm([i-1:i+fcstH_M],:);
   
    FcstQ(i,1:2:end,:) = Res.X_sm([i-Month_i:3:i-Month_i+3*(fcstH_Q+1)],nM+1:end);
    
    % second unbalancedeness pattern
    eval(['UnbP = UnbPattM',int2str(Month_i),'_2;']);
    
    X = Data(1:i-nn,:);
    temp = X(end-nUnb+1:end,:);
    temp(isnan(UnbP)) = nan;
    X(end-nUnb+1:end,:) = temp;
    X = [X;nan(n_nan,nM+nQ)];
    
    eval(['Res = EM_DFM_SS',P.method,P.idio,P.restr,'(X,P);'])
    
    XM = DataMM(1:i-nn,:);
    temp = XM(end-nUnb+1:end,:);
    temp(isnan(UnbP(:,1:nM))) = nan;
    XM(end-nUnb+1:end,:) = temp;
    
    tempFcstM = cnvrt2Mgr(Res.X_sm(:,1:nM),X(:,1:nM),XM,TransfM(:,2),P.transf);
    FcstM1_2(i,:,:) = tempFcstM([i-1:i+fcstH_M],:);
    tfm = filter(ones(12,1),1,tempFcstM);
    FcstM12_2(i,:,:) = tfm([i-1:i+fcstH_M],:);
   
    FcstQ(i,2:2:end,:) = Res.X_sm([i-Month_i:3:i-Month_i+3*(fcstH_Q+1)],nM+1:end);
    
   
% end



% forecasts for the quarterly variables
FcstQQ = [];
TrueQQ = [];
for k = iS:3:iE-8
    FcstQQ = cat(1,FcstQQ,cat(2,FcstQ(k,5:6,:),FcstQ(k+1,5:6,:),FcstQ(k+2,5:6,:),...
        FcstQ(k+3,3:4,:),FcstQ(k+4,3:4,:),FcstQ(k+5,3:4,:),...
        FcstQ(k+6,1:2,:),FcstQ(k+7,1:2,:),FcstQ(k+8,1:2,:)));
    TrueQQ = [TrueQQ;Data(k+5,nM+1:end)];
end
DateQQ = Dates(iS+5:3:iE-3);
DateQQ_V = DatesV(iS+5:3:iE-3,:);

% Forecasts for the annual growth rates of month;y variables
TrueMM12 = filter(ones(12,1),1,DataMM);
TrueMM12 = TrueMM12(iS+12:iE-1,:);
FcstMM12 = nan(size(TrueMM12,1),28,nM);
for j = 1:14
    FcstMM12(:,2*j-1,:) = FcstM12_1(iS-1+j:iE-14+j,15-j,:);
    FcstMM12(:,2*j,:) = FcstM12_2(iS-1+j:iE-14+j,15-j,:);
end
DateMM = Dates(iS+12:iE-1);
DateMM_V = DatesV(iS+12:iE-1,:);



% RMSFE for GDP
iVar = find(ismember(SeriesQ,'GDP'));
iS = find(DateQQ_V(:,1) == P.StartEvQ(1) & DateQQ_V(:,2) == P.StartEvQ(2));
iE = find(DateQQ_V(:,1) == P.EndEvQ(1) & DateQQ_V(:,2) == P.EndEvQ(2));
rmsfeGDP  = sqrt(mean((FcstQQ(iS:iE,1:14,1) - repmat(TrueQQ(iS:iE,1),1,14)).^2)');

% RMSFE for HICP
iVar = find(ismember(SeriesQ,'HICP total'));
iS = find(DateMM_V(:,1) == P.StartEvM(1) & DateMM_V(:,2) == P.StartEvM(2));
iE = find(DateMM_V(:,1) == P.EndEvM(1) & DateMM_V(:,2) == P.EndEvM(2));
rmsfeHICP = sqrt(mean((FcstMM12(iS:iE,:,1) - repmat(TrueMM12(iS:iE,1),1,28)).^2)');


% save the results
datafile = ['fcst',P.DF,P.method,P.idio,strrep(int2str(P.r),' ',''),int2str(P.p)];
eval(['save ',datafile,...
    ' FcstMM12 TrueMM12 DateMM DateMM_V FcstQQ TrueQQ DateQQ DateQQ_V rmsfe* Data Series* Group Dates DatesV P'])
    



