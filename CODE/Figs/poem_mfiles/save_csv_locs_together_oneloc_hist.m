% Save as csv for using in R

clear all
close all

%dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
%dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ_TrefOrig/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ_TrefOrig/';

load([dpath 'Oneloc_hist_phenol_all.mat'])
load([dpath 'Oneloc_hist_all.mat'])

%% Recruitment
x=1:length(sp); 
yr=x((end-(30*365)+1):end);
SPL=Mlp(yr,:);
SFL=Mmf(yr,:);
SDL=Mld(yr,:);
PA=lp(yr,:);
FA=mf(yr,:);
DA=ld(yr,:);

pSPL=pMlp(yr,:);
pSFL=pMmf(yr,:);
pSDL=pMld(yr,:);
pPA=plp(yr,:);
pFA=pmf(yr,:);
pDA=pld(yr,:);

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
    PLy(n,:) = nansum(SPL(st(n):en(n),:));
    FLy(n,:) = nansum(SFL(st(n):en(n),:));
    DLy(n,:) = nansum(SDL(st(n):en(n),:));
    PAy(n,:) = nansum(PA(st(n):en(n),:));
    FAy(n,:) = nansum(FA(st(n):en(n),:));
    DAy(n,:) = nansum(DA(st(n):en(n),:));
    pPLy(n,:) = nansum(pSPL(st(n):en(n),:));
    pFLy(n,:) = nansum(pSFL(st(n):en(n),:));
    pDLy(n,:) = nansum(pSDL(st(n):en(n),:));
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

%% FIGURES

close all
for n=1:length(spots)
    loc=spots{n};
    
    f10=figure(10);
    subplot(2,3,n)
    plot(1976:2005,FLy(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pFLy(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    ylabel('Forage fish recruits (g m^-^2)')
    if (n==5)
        legend('C','S')
        legend('location','northeast')
    end
    title(loc)
    
    f11=figure(11);
    subplot(2,3,n)
    plot(1976:2005,PLy(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pPLy(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    title(loc)
    if (n==5)
        legend('C','S')
        legend('location','northeast')
    end
    ylabel('Pelagic piscivore recruits (g m^-^2)')
    
    f12=figure(12);
    subplot(2,3,n)
    plot(1976:2005,DLy(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pDLy(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    title(loc)
    if (n==5)
        legend('C','S')
        legend('location','northeast')
    end
    ylabel('Demersal piscivore recruits (g m^-^2)')
end
print(f10,'-dpng',[fpath 'Forage_oneloc_compare_recruitment.png'])
print(f11,'-dpng',[fpath 'Pisc_oneloc_compare_recruitment.png'])
print(f12,'-dpng',[fpath 'Dem_oneloc_compare_recruitment.png'])

%% Convert from biomass to numbers
cid=csvread('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/grid_360x200_id_locs_area.csv',1,1);
A = cid(:,4);

PLn = NaN*ones(30,6);
FLn = PLn;
DLn = PLn;
pPLn = PLn;
pFLn = PLn;
pDLn = PLn;
for s=1:6
    PLn(:,s) = round(PLy(:,s) * A(s) / 2529.822128134704) +1;
    FLn(:,s) = round(FLy(:,s) * A(s) / 2.5298221281347053) +1;
    DLn(:,s) = round(DLy(:,s) * A(s) / 2529.822128134704) +1;
    pPLn(:,s) = round(pPLy(:,s) * A(s) / 2529.822128134704) +1;
    pFLn(:,s) = round(pFLy(:,s) * A(s) / 2.5298221281347053) +1;
    pDLn(:,s) = round(pDLy(:,s) * A(s) / 2529.822128134704) +1;
end
%%
for n=1:length(spots)
    loc=spots{n};
    
    f1=figure(1);
    subplot(2,3,n)
    plot(1976:2005,log(FLn(:,n)),'Linewidth',2); hold on;
    plot(1976:2005,log(pFLn(:,n)),'Linewidth',2); hold on;
    xlim([1976 2005])
    %ylim([0 22])
    ylabel('log Forage fish recruits')
%     if (n==5)
%         legend('C','S')
%         legend('location','northeast')
%     end
    title(loc)
    
    f2=figure(2);
    subplot(2,3,n)
    plot(1976:2005,log(PLn(:,n)),'Linewidth',2); hold on;
    plot(1976:2005,log(pPLn(:,n)),'Linewidth',2); hold on;
    xlim([1976 2005])
    %ylim([0 22])
    title(loc)
%     if (n==5)
%         legend('C','S')
%         legend('location','northeast')
%     end
    ylabel('log Pelagic piscivore recruits')
    
    f3=figure(3);
    subplot(2,3,n)
    plot(1976:2005,log(DLn(:,n)),'Linewidth',2); hold on;
    plot(1976:2005,log(pDLn(:,n)),'Linewidth',2); hold on;
    xlim([1976 2005])
    %ylim([0 22])
    title(loc)
%     if (n==5)
%         legend('C','S')
%         legend('location','northeast')
%     end
    ylabel('log Demersal piscivore recruits')
end
print(f1,'-dpng',[fpath 'Forage_oneloc_compare_recruit_num.png'])
print(f2,'-dpng',[fpath 'Pisc_oneloc_compare_recruit_num.png'])
print(f3,'-dpng',[fpath 'Dem_oneloc_compare_recruit_num.png'])

%% save numbers
year=[1976:2005]';
fL = table(year,FLn(:,1),pFLn(:,1),FLn(:,2),pFLn(:,2),FLn(:,3),pFLn(:,3),...
    FLn(:,4),pFLn(:,4),FLn(:,5),pFLn(:,5),FLn(:,6),pFLn(:,6),...
    'VariableNames',{'Year','GBc','GBp','EBSc','EBSp','OSPc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

pL = table(year,PLn(:,1),pPLn(:,1),PLn(:,2),pPLn(:,2),PLn(:,3),pPLn(:,3),...
    PLn(:,4),pPLn(:,4),PLn(:,5),pPLn(:,5),PLn(:,6),pPLn(:,6),...
    'VariableNames',{'Year','GBc','GBp','EBSc','EBSp','OSPc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

dL = table(year,DLn(:,1),pDLn(:,1),DLn(:,2),pDLn(:,2),DLn(:,3),pDLn(:,3),...
    DLn(:,4),pDLn(:,4),DLn(:,5),pDLn(:,5),DLn(:,6),pDLn(:,6),...
    'VariableNames',{'Year','GBc','GBp','EBSc','EBSp','OSDc','OSPp','HOTc',...
    'HOTp','BATSc','BATSp','NSc','NSp'});

writetable(fL,[dpath 'onelocs_hist_forage_larv_num.csv'])
writetable(pL,[dpath 'onelocs_hist_pisc_larv_num.csv'])
writetable(dL,[dpath 'onelocs_hist_dem_larv_num.csv'])
