% Save as csv for using in R

clear all
close all

%dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ/';

load([dpath 'Oneloc_hist_phenol_all.mat'])
load([dpath 'Oneloc_hist_all.mat'])

%% Recruitment
t=1:length(sp);
yr=t((end-(30*365)+1):end);
PL=Mlp(yr,:);
FL=Mmf(yr,:);
DL=Mld(yr,:);
PA=lp(yr,:);
FA=mf(yr,:);
DA=ld(yr,:);

pPL=pMlp(yr,:);
pFL=pMmf(yr,:);
pDL=pMld(yr,:);
pPA=plp(yr,:);
pFA=pmf(yr,:);
pDA=pld(yr,:);

%%
st=1:365:length(yr);
en=365:365:length(yr);
PLy = NaN*ones(30,6);
FLy = PLy;
DLy = PLy;
PAy = PLy;
FAy = PLy;
DAy = PLy;
pPLy = PLy;
pFLy = PLy;
pDLy = PLy;
pPAy = PLy;
pFAy = PLy;
pDAy = PLy;
for n=1:30
    PLy(n,:) = nansum(PL(st(n):en(n),:));
    FLy(n,:) = nansum(FL(st(n):en(n),:));
    DLy(n,:) = nansum(DL(st(n):en(n),:));
    PAy(n,:) = nansum(PA(st(n):en(n),:));
    FAy(n,:) = nansum(FA(st(n):en(n),:));
    DAy(n,:) = nansum(DA(st(n):en(n),:));
    pPLy(n,:) = nansum(pPL(st(n):en(n),:));
    pFLy(n,:) = nansum(pFL(st(n):en(n),:));
    pDLy(n,:) = nansum(pDL(st(n):en(n),:));
    pPAy(n,:) = nansum(pPA(st(n):en(n),:));
    pFAy(n,:) = nansum(pFA(st(n):en(n),:));
    pDAy(n,:) = nansum(pDA(st(n):en(n),:));
end

%%
fL = table(FLy(:,1),pFLy(:,1),FLy(:,2),pFLy(:,2),FLy(:,3),pFLy(:,3),...
    FLy(:,4),pFLy(:,4),FLy(:,5),pFLy(:,5),FLy(:,6),pFLy(:,6),...
    'VariableNames',{'GBc','GBp','EBSc','EBSp','OSPc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

fA = table(FAy(:,1),pFAy(:,1),FAy(:,2),pFAy(:,2),FAy(:,3),pFAy(:,3),...
    FAy(:,4),pFAy(:,4),FAy(:,5),pFAy(:,5),FAy(:,6),pFAy(:,6),...
    'VariableNames',{'GBc','GBp','EBSc','EBSp','OSPc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

pL = table(PLy(:,1),pPLy(:,1),PLy(:,2),pPLy(:,2),PLy(:,3),pPLy(:,3),...
    PLy(:,4),pPLy(:,4),PLy(:,5),pPLy(:,5),PLy(:,6),pPLy(:,6),...
    'VariableNames',{'GBc','GBp','EBSc','EBSp','OSPc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

pA = table(PAy(:,1),pPAy(:,1),PAy(:,2),pPAy(:,2),PAy(:,3),pPAy(:,3),...
    PAy(:,4),pPAy(:,4),PAy(:,5),pPAy(:,5),PAy(:,6),pPAy(:,6),...
    'VariableNames',{'GBc','GBp','EBSc','EBSp','OSPc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

dL = table(DLy(:,1),pDLy(:,1),DLy(:,2),pDLy(:,2),DLy(:,3),pDLy(:,3),...
    DLy(:,4),pDLy(:,4),DLy(:,5),pDLy(:,5),DLy(:,6),pDLy(:,6),...
    'VariableNames',{'GBc','GBp','EBSc','EBSp','OSDc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

dA = table(DAy(:,1),pDAy(:,1),DAy(:,2),pDAy(:,2),DAy(:,3),pDAy(:,3),...
    DAy(:,4),pDAy(:,4),DAy(:,5),pDAy(:,5),DAy(:,6),pDAy(:,6),...
    'VariableNames',{'GBc','GBp','EBSc','EBSp','OSPc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

%%
writetable(fL,[dpath 'onelocs_hist_forage_larv.csv'])
writetable(fA,[dpath 'onelocs_hist_forage_adult.csv'])
writetable(pL,[dpath 'onelocs_hist_pisc_larv.csv'])
writetable(pA,[dpath 'onelocs_hist_pisc_adult.csv'])
writetable(dL,[dpath 'onelocs_hist_dem_larv.csv'])
writetable(dA,[dpath 'onelocs_hist_dem_adult.csv'])

