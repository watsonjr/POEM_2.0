%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/'];
% npath8 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/'];
% npath9 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/'];
% npath10 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath11 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/'];
% npath12 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
% npath13 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath14 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/'];
% npath15 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/'];
% npath16 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/'];
% npath17 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath18 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/'];
% npath19 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
% npath20 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath21 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/'];
npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE025/'];
npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/'];
npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP25/'];
npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP50/'];
npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP75/'];
npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE075/'];
npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10/'];
npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE15_BP25/'];
npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE15_BP50/'];
npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE15_BP75/'];
npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20/'];
npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE25_BP25/'];
npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE25_BP50/'];
npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE25_BP75/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30/'];
% npath1 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/'];
% npath2 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/'];
% npath3 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath4 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20/'];
% npath5 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25/'];
% npath6 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30/'];


fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

% dp = {dpath1; dpath2; dpath3; dpath4; dpath5; dpath6};
% sims = {'f05','f10','f15','f20','f30','f40'};
% cfile = 'PDc_MFeqMP_fcrit_comp';
% dp = {dpath3;dpath4;dpath2;dpath5;dpath6;dpath7;dpath8};
% sims = {'MP0.50','MP0.75','Meq','MF1.5','MF2.0','MF2.5','MF3.0'};
% cfile = 'PDc_MFbetter_enc_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5;dpath6;dpath7};
% sims = {'con','T','L','T&L','simpQ','compQ1','compQ1'};
% cfile = 'PDc_MFeqMP_mort_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5};
% sims = {'full','redSM','T-redSM','redS','T-redS'};
% cfile = 'PDc_MFeqMP_cannibal_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5;dpath6;dpath7;dpath8};
% sims = {'0.0+N','0.01+N','0.025+N','0.05+N','0.10+N','0.20+N','0.30+N','0.40+N'};
% cfile = 'PDc_MFeqMP_fishing_nmort_comp';
% dp = {dpath1;dpath9;dpath10;dpath11;dpath12;dpath13;dpath14;dpath15};
% sims = {'0.0','0.01','0.025','0.05','0.10','0.20','0.30','0.40'};
% cfile = 'PDc_MFeqMP_fishing_comp';
% dp = {dpath2;dpath3;dpath4;dpath1;dpath5;dpath6;dpath7};
% sims = {'0.0','0.1','0.5','1.0','MP0.0','MP0.1','MP0.5'};
% cfile = 'PDc_MFeqMP_MZ_comp';
% dp = {dpath2;dpath3;dpath1;dpath4;dpath5};
% sims = {'0.5','0.75','1.0','1.25','1.5'};
% cfile = 'PDc_MFeqMP_Mcann_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5};
% sims = {'none','const','T','L','T&L'};
% cfile = 'PDc_MFeqMP_mort2_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;...
%     npath10;npath11};
% sims = {'0.0MZ01','0.025MZ01','0.05MZ01','0.10MZ01','0.20MZ01','0.30MZ01',...
%     '0.40MZ01','0.45MZ01','0.50MZ01','0.55MZ01','0.60MZ01','0.70MZ01'};
% cfile = 'PDc_MFeqMP_MZ01_fishing_comp';
% dp = {npath0;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9};
% sims = {'0.0','0.05','0.10','0.20','0.30','0.40','0.50','0.60','0.70'};
% cfile = 'PDc_MFeqMP_MZ01_halfM_fishing_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5};
% sims = {'0.0','0.1','0.25','0.5','0.75','1.0'};
% cfile = 'Dc_MFeqMP_NOnmort_MZpref_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;...
%     npath9};
% sims = {'0.025','0.05','0.075','0.10','0.20','0.30','0.40','0.50',...
%     '0.60','0.70'};
% cfile = 'Dc_MFeqMP_MZ01_fishing_comp';
% dp = {npath1;npath2;npath3;npath4;npath5};
% sims = {'none','const','T','L','T&L'};
% cfile = 'Dc_MFeqMP_MZ01_mort_comp';
% dp = {npath1;npath2;npath3;npath4;npath5};
% sims = {'0.1','0.2','0.3','0.4','0.5'};
% cfile = 'Dc_MFeqMP_MZ01_fcrit_comp';
% dp = {npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11};
% sims = {'0.10','0.20','0.30','0.40','0.50','0.60','0.70','0.80','0.90'};
% cfile = 'Dc_MFeqMP_MZ01_fishing_half_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% sims = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
% cfile = 'Dc_MFeqMP_MZ01_reproeff_comp';
% dp = {npath1;npath2;npath3};
% sims = {'noPDc','PDc','Dc'};
% cfile = 'K&Hparams_all_MFeqMP_MZ01_PDc_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11;npath12};
% sims = {'no25','D25','PD25','no50','D50','PD50','no75','D75','PD75','no100','D100','PD100'};
% cfile = 'K&Hparams_all_MFeqMP_MZ01_PDc_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7};
% sims = {'0.05','0.075','0.1','0.2','0.3','0.4','0.5'};
% cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit_comp';
% dp = {npath8;npath9;npath10;npath11;npath12;npath13;npath14};
% sims = {'0.05','0.075','0.1','0.2','0.3','0.4','0.5'};
% cfile = 'PDc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit_comp';
% dp = {npath15;npath16;npath17;npath18;npath19;npath20;npath21};
% sims = {'0.05','0.075','0.1','0.2','0.3','0.4','0.5'};
% cfile = 'NoPDc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit_comp';
% dp = {npath17;npath3;npath10};
% sims = {'noPDc','Dc','PDc'};
% cfile = 'Hartvig_cmax-metab_MFeqMP_MZ01_fcrit10_PDc_comp';
% dp = {npath21;npath7;npath14};
% sims = {'noPDc','Dc','PDc'};
% cfile = 'Hartvig_cmax-metab_MFeqMP_MZ01_fcrit50_PDc_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7};
% sims = {'.025','.05','.075','.1','.15','.2','.3'};
% cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit30_BentEff_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6};
% sims = {'.05','.10','.15','.20','.25','.30'};
% cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit40_BentEff_comp';
dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11;npath12;...
    npath13;npath14;npath15};
sims = {'.025-1','.05-1/4','.05-1/2','.05-3/4','.05-1','.075','.1-1','.15-1/4','.15-1/2',...
    '.15-3/4','.15-1','.2-1','.25-1/4','.25-1/2','.25-3/4'};
cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit30_BentEff_BentPref_comp';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

%%

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met'};
cols=cols';

ndp = length(dp);

load('cmap_ppt_angles.mat')

%% Predicted biomass
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/grid_360x200_id_locs_area_dep.mat','depth')
% Using linear equations based on depth from Wei et al. 2010
% log10 biomass in mg C /m2
%non-log10
meio = 10.^(2.18-(2.39e-4)*depth);
macro = 10.^(3.05-(5.15e-4)*depth);
mega = 10.^(1.81-(3.07e-4)*depth);
%convert from mg C to g wet weight
meio = meio * 9e-3;
macro = macro * 9e-3;
mega = mega * 9e-3;

allb = meio+macro+mega;

%%
for i=1:length(spots)
    loc = spots{i};
    lname = [sname2 loc '_'];
    %%
    for s=1:ndp
        dpath = char(dp(s));
        load([dpath sname sname2 'lastyr_bent_biom.mat']);
        load([dpath sname sname2 'lastyr_sum_mean_biom.mat']);
        
        %Demersal fishes (M&L)
        dem = nansum(squeeze(all_mean(2:3,3,:)));
        BF = Bmean + dem';
        
        %% Bent biom by location
        f1=figure(1);
        plot(s,log10(allb(i)),'.r','MarkerSize',25); hold on;
        plot(s,log10(Bmean(i)),'.k','MarkerSize',25); hold on;
        plot(s,log10(BF(i)),'.b','MarkerSize',25); hold on;
        ylim([-2 2])
        xlim([0 1+ndp])
        set(gca,'XTick',1:ndp,'XTickLabel',[]);
        for t=1:ndp
            text(t,-2.1,sims{t},'HorizontalAlignment','right','Rotation',45)
        end
        if (s==ndp)
            ylabel('log10 Mean Benthic Biom (g m^-^2) in final year')
            title([loc ' Benthic biomass'])
            legend('predict','modelB','modelBF')
            stamp(cfile)
        end
       
        f2=figure(2);
        subplot(3,3,i)
        plot(s,log10(allb(i)),'.r','MarkerSize',25); hold on;
        plot(s,log10(Bmean(i)),'.k','MarkerSize',25); hold on;
        plot(s,log10(BF(i)),'.b','MarkerSize',25); hold on;
        ylim([-2 2])
        xlim([0 1+ndp])
        set(gca,'XTick',1:ndp,'XTickLabel',[]);
        for t=1:ndp
            text(t,-2.1,sims{t},'HorizontalAlignment','right','Rotation',45)
        end
        title(loc)
        if (s==ndp)
            stamp(cfile)
        end
        if (i==4)
            ylabel('log10 Mean Benthic Biom (g m^-^2) in final year')
        end
        if (i==length(spots))
            legend('p','B','BF')
        end
        
        
        
    end
    print(f1,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'mean_bent_biom.png'])
    close(f1)
    
end
print(f2,'-dpng',[fpath sname sname2 cfile '_mean_bent_biom.png'])
