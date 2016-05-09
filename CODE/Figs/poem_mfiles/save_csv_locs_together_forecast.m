% Save as csv for using in R

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';

load('Oneloc_fore_pris_phenol_all.mat')
load('Oneloc_fore_pris_all.mat')

%% Recruitment
t=1:length(sp);
yr=t(1:end);
PL=sp(yr,:);
FL=sf(yr,:);
PJ=mp(yr,:);
PA=lp(yr,:);
FA=mf(yr,:);
pPL=psp(yr,:);
pFL=psf(yr,:);
pPJ=pmp(yr,:);
pPA=plp(yr,:);
pFA=pmf(yr,:);

st=1:365:length(yr);
en=365:365:length(yr);
PLy = NaN*ones(95,6);
FLy = PLy;
PJy = PLy;
PAy = PLy;
FAy = PLy;
pPLy = PLy;
pFLy = PLy;
pPJy = PLy;
pPAy = PLy;
pFAy = PLy;
for n=1:95
    PLy(n,:) = nansum(PL(st(n):en(n),:));
    FLy(n,:) = nansum(FL(st(n):en(n),:));
    PJy(n,:) = nansum(PJ(st(n):en(n),:));
    PAy(n,:) = nansum(PA(st(n):en(n),:));
    FAy(n,:) = nansum(FA(st(n):en(n),:));
    pPLy(n,:) = nansum(pPL(st(n):en(n),:));
    pFLy(n,:) = nansum(pFL(st(n):en(n),:));
    pPJy(n,:) = nansum(pPJ(st(n):en(n),:));
    pPAy(n,:) = nansum(pPA(st(n):en(n),:));
    pFAy(n,:) = nansum(pFA(st(n):en(n),:));
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

pJ = table(PJy(:,1),pPJy(:,1),PJy(:,2),pPJy(:,2),PJy(:,3),pPJy(:,3),...
    PJy(:,4),pPJy(:,4),PJy(:,5),pPJy(:,5),PJy(:,6),pPJy(:,6),...
    'VariableNames',{'GBc','GBp','EBSc','EBSp','OSPc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

pA = table(PAy(:,1),pPAy(:,1),PAy(:,2),pPAy(:,2),PAy(:,3),pPAy(:,3),...
    PAy(:,4),pPAy(:,4),PAy(:,5),pPAy(:,5),PAy(:,6),pPAy(:,6),...
    'VariableNames',{'GBc','GBp','EBSc','EBSp','OSPc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

%%
writetable(fL,[dpath 'onelocs_fore_forage_larv.csv'])
writetable(fA,[dpath 'onelocs_fore_forage_adult.csv'])
writetable(pL,[dpath 'onelocs_fore_pisc_larv.csv'])
writetable(pJ,[dpath 'onelocs_fore_pisc_juve.csv'])
writetable(pA,[dpath 'onelocs_fore_pisc_adult.csv'])

