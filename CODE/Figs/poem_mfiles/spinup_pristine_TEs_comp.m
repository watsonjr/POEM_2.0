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
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE025/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP25/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP50/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP75/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE075/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10_BP25/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10_BP50/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10_BP75/'];
% npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10/'];
% npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE15_BP25/'];
% npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE15_BP50/'];
% npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE15_BP75/'];
% npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
% npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_BP25/'];
% npath16 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_BP50/'];
% npath17 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_BP75/'];
% npath18 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20/'];
% npath19 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE25_BP25/'];
% npath20 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE25_BP50/'];
% npath21 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE25_BP75/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_BP25/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_BP50/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_BP75/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10_BP25/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10_BP50/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10_BP75/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE15_BP25/'];
% npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE15_BP50/'];
% npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE15_BP75/'];
% npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_BP25/'];
% npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_BP50/'];
% npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_BP75/'];
% npath16 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20/'];
% npath17 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25_BP25/'];
% npath18 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25_BP50/'];
% npath19 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25_BP75/'];
% npath20 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25/'];
% npath21 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30/'];
npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE05_BP25/'];
npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE05_BP50/'];
npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE05_BP75/'];
npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE05_BP100/'];
npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE10_BP25/'];
npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE10_BP50/'];
npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE10_BP75/'];
npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE10_BP100/'];
npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE15_BP25/'];
npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE15_BP50/'];
npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE15_BP75/'];
npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE15_BP100/'];
npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE20_BP25/'];
npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE20_BP50/'];
npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE20_BP75/'];
npath16 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE20_BP100/'];


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
% sims = {'.05','.075','.1','.2','.3','.4','.5'};
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
% dp = {npath1;npath2;npath3;npath4;npath5;npath6};
% sims = {'.025','.05','.075','.1','.15','.2'};
% cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit30_BentEff_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6};
% sims = {'.05','.1','.15','.2','.25','.3'};
% cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit40_BentEff_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11;npath12;...
%     npath13;npath14;npath15;npath16;npath17;npath18;npath19;npath20;npath21};
% sims = {'.025-100','.05-25','.05-50','.05-75','.05-100','.075-100','.1-25','.1-50','.1-75',...
%     '.1-100','.15-25','.15-50','.15-75','.15-100','.2-25','.2-50','.2-75','.2-100','.25-25',...
%     '.25-50','.25-75'};
% cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit30_BentEff_BentPref_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11;npath12;...
%     npath13;npath14;npath15;npath16;npath17;npath18;npath19;npath20};
% sims = {'.05-25','.05-50','.05-75','.05-100','.1-25','.1-50','.1-75','.1-100','.15-25','.15-50',...
%     '.15-75','.15-100','.2-25','.2-50','.2-75','.2-100','.25-25','.25-50','.25-75','.25-100',};
% cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit40_BentEff_BentPref_comp';
dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11;npath12;...
    npath13;npath14;npath15;npath16};
sims = {'.05-25','.05-50','.05-75','.05-100','.1-25','.1-50','.1-75',...
    '.1-100','.15-25','.15-50','.15-75','.15-100','.2-25','.2-50','.2-75','.2-100'};
cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit35_BentEff_BentPref_comp';

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

%%
for i=1:length(spots)
    close all
    loc = spots{i};
    lname = [sname2 loc '_'];
    %%
    for s=1:ndp
        dpath = char(dp(s));
        load([dpath sname sname2 'lastyr_TEs.mat']);
        
        %% TEs by location
        %         f1=figure(1);
        %         subplot(2,1,1)
        %         plot(s,TEcon(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Medium/Small')
        %         title([loc ' mean con TE in final year'])
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %
        %         subplot(2,1,2)
        %         plot(s,TEcon(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Large/Medium')
        %         xlabel('sim')
        %
        %         f11=figure(11);
        %         subplot(2,1,1)
        %         plot(s,TEcon(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([0 1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Medium/Small')
        %         title([loc ' mean con TE in final year'])
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %
        %         subplot(2,1,2)
        %         plot(s,TEcon(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([0 1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Large/Medium')
        %         xlabel('sim')
        %
        %         %% Production
        %         f2=figure(2);
        %         subplot(2,1,1)
        %         plot(s,TEprod(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Medium/Small')
        %         title([loc ' mean prod TE in final year'])
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %
        %         subplot(2,1,2)
        %         plot(s,TEprod(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Large/Medium')
        %         xlabel('sim')
        %
        %         f12=figure(12);
        %         subplot(2,1,1)
        %         plot(s,TEprod(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([0 1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Medium/Small')
        %         title([loc ' mean prod TE in final year'])
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %
        %         subplot(2,1,2)
        %         plot(s,TEprod(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([0 1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Large/Medium')
        %         xlabel('sim')
        %
        %         %% effective TE
        %         f3=figure(3);
        %         subplot(2,1,1)
        %         plot(s,TEeff(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Medium')
        %         title([loc ' mean eff TE in final year'])
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %
        %         subplot(2,1,2)
        %         plot(s,TEeff(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         ylabel('Large')
        %         xlabel('sim')
        
        %% TL prey
        f20=figure(20);
        subplot(2,1,1)
        plot(s-0.25,TLprey(4,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,TLprey(5,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,TLprey(6,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',[]);
        for t=1:ndp
            text(t,2.95,sims{t},'HorizontalAlignment','right','Rotation',90)
        end
        title('Mean TL in final year')
        ylabel('M')
        if (s==ndp)
            stamp(cfile)
        end
        
        subplot(2,1,2)
        plot(s,TLprey(7,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,TLprey(8,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        set(gca,'XTick',1:ndp,'XTickLabel',[]);
        for t=1:ndp
            text(t,2.95,sims{t},'HorizontalAlignment','right','Rotation',90)
        end
        ylabel('L')
        
        
    end
    %     print(f1,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'TEs_con.png'])
    %     print(f2,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'TEs_prod.png'])
    %     print(f3,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'TEs_eff.png'])
    %     print(f11,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'TEs_con_axes.png'])
    %     print(f12,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'TEs_prod_axes.png'])
    print(f20,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'TLs.png'])
    
end

%%
for i=1:length(spots)
    loc = spots{i};
    lname = [sname2 loc '_'];
    %%
    for s=1:ndp
        dpath = char(dp(s));
        load([dpath sname sname2 'lastyr_TEs.mat']);
        
        %% TEs subplot
        %         f4=figure(4);
        %         subplot(3,3,i)
        %         plot(s,TEcon(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Medium/Small consumption')
        %         end
        %
        %         f5=figure(5);
        %         subplot(3,3,i)
        %         plot(s,TEcon(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Large/Medium consumption')
        %         end
        %
        %         f13=figure(13);
        %         subplot(3,3,i)
        %         plot(s,TEcon(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([0 1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Medium/Small consumption')
        %         end
        %
        %         f14=figure(14);
        %         subplot(3,3,i)
        %         plot(s,TEcon(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([0 1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Large/Medium consumption')
        %         end
        %
        %         %% production
        %         f6=figure(6);
        %         subplot(3,3,i)
        %         plot(s,TEprod(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Medium/Small production')
        %         end
        %
        %         f7=figure(7);
        %         subplot(3,3,i)
        %         plot(s,TEprod(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Large/Medium production')
        %         end
        %
        %         f15=figure(15);
        %         subplot(3,3,i)
        %         plot(s,TEprod(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([0 1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Medium/Small production')
        %         end
        %
        %         f16=figure(16);
        %         subplot(3,3,i)
        %         plot(s,TEprod(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([0 1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Large/Medium production')
        %         end
        %
        %         %% effective
        %         f8=figure(8);
        %         subplot(3,3,i)
        %         plot(s,TEeff(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Medium eff TE')
        %         end
        %
        %         f9=figure(9);
        %         subplot(3,3,i)
        %         plot(s,TEeff(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Large eff TE')
        %         end
        %
        %         f18=figure(18);
        %         subplot(3,3,i)
        %         plot(s,TEeff(1,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([-25e-3 10e-2])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Medium eff TE')
        %         end
        %
        %         f19=figure(19);
        %         subplot(3,3,i)
        %         plot(s,TEeff(2,i),'.k','MarkerSize',25); hold on;
        %         xlim([0 ndp+1])
        %         ylim([-1e-3 15e-3])
        %         set(gca,'XTick',1:ndp,'XTickLabel',sims)
        %         title(loc)
        %         if (s==ndp)
        %             stamp(cfile)
        %         end
        %         if (i==4)
        %             ylabel('Large eff TE')
        %         end
        %
        %% TL prey
        f21=figure(21);
        subplot(3,3,i)
        plot(s-0.25,TLprey(4,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,TLprey(5,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,TLprey(6,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([3 4])
        set(gca,'XTick',1:ndp,'XTickLabel',[]);
        for t=1:ndp
            text(t,2.95,sims{t},'HorizontalAlignment','right','Rotation',45)
        end
        title(loc)
        if (i==4)
            ylabel('Mean M TL in final year')
        end
        if (s==ndp)
            stamp(cfile)
        end
        
        f22=figure(22);
        subplot(3,3,i)
        plot(s,TLprey(7,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,TLprey(8,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([3 5])
        set(gca,'XTick',1:ndp,'XTickLabel',[]);
        for t=1:ndp
            text(t,2.95,sims{t},'HorizontalAlignment','right','Rotation',45)
        end
        title(loc)
        if (i==4)
            ylabel('Mean L TL in final year')
        end
        if (s==ndp)
            stamp(cfile)
        end
        
    end
end
% print(f4,'-dpng',[fpath sname sname2 cfile '_TEs_con_M.png'])
% print(f5,'-dpng',[fpath sname sname2 cfile '_TEs_con_L.png'])
% print(f6,'-dpng',[fpath sname sname2 cfile '_TEs_prod_M.png'])
% print(f7,'-dpng',[fpath sname sname2 cfile '_TEs_prod_L.png'])
% print(f8,'-dpng',[fpath sname sname2 cfile '_TEs_eff_M.png'])
% print(f9,'-dpng',[fpath sname sname2 cfile '_TEs_eff_L.png'])
% print(f13,'-dpng',[fpath sname sname2 cfile '_TEs_con_M_axes.png'])
% print(f14,'-dpng',[fpath sname sname2 cfile '_TEs_con_L_axes.png'])
% print(f15,'-dpng',[fpath sname sname2 cfile '_TEs_prod_M_axes.png'])
% print(f16,'-dpng',[fpath sname sname2 cfile '_TEs_prod_L_axes.png'])
% print(f18,'-dpng',[fpath sname sname2 cfile '_TEs_eff_M_axes.png'])
% print(f19,'-dpng',[fpath sname sname2 cfile '_TEs_eff_L_axes.png'])
print(f21,'-dpng',[fpath sname sname2 cfile '_TLs_M.png'])
print(f22,'-dpng',[fpath sname sname2 cfile '_TLs_L.png'])


