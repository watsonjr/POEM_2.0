%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_flev1e4/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_flev1e4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_flev4e4/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_flev4e4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_flev8e4/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_flev8e4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e4_NoWgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e4_NoWgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoMetab_TrefO_1e4/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoMetab_TrefO_1e4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoMetab_TrefO_1e4_HalfC/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoMetab_TrefO_1e4_HalfC/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e4_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e4_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_4e4_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_4e4_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_8e4_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_8e4_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e5_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e5_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e6_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e6_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_Cmax_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_Cmax_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoMet_TrefO_Cmax_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoMet_TrefO_Cmax_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e6_noC&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e6_noC&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_TrefO_KHparams_cmax-metab/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_TrefO_KHparams_cmax-metab/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_TrefO_KHparams_cmax-metab_MFeatS/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_TrefO_KHparams_cmax-metab_MFeatS/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MprefLZoverMZ/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MprefLZoverMZ/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFprefZ/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFprefZ/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit01/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit01/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFmet/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFmet/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFpred/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFpred/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LTmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LTmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_simpQmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_simpQmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_compQmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_compQmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_bioQmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_bioQmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF50/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF50/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF75/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF75/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_sameA/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_sameA/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1_Tmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2_Tmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Tmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Tmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_Fenc2x_Tmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_Fenc2x_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Lmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Lmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Lmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Lmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit20/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit20/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit30/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit30/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit40/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit40/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05_Tmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc15/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc15/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc20/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc20/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc25/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc25/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc30/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc30/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc075/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc075/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc050/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc050/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc025/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc025/';

cfile = 'fcrit10_MPenc050';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

load([dpath sname sname2 'consump.mat'],'mclev','Zcon');

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod'};
cols=cols';

load('cmap_ppt_angles.mat')

%% Zoop con

figure(12)
subplot(2,1,1)
bar(Zcon(:,1),'k')
ylim([0 1])
set(gca,'XTickLabel',spots);
title('Med zoo')
ylabel('Fraction of times overconsumed')

subplot(2,1,2)
bar(Zcon(:,2),'k')
ylim([0 1])
set(gca,'XTickLabel',spots);
title('Large zoo')
ylabel('Fraction of times overconsumed')
stamp(cfile)
print('-dpng',[fpath sname sname2 'All_oneloc_zoo_con.png'])

%%
Psum = NaN*ones(3,length(spots));
Fsum = NaN*ones(2,length(spots));
Dsum = Psum;
Pmean = Psum;
Fmean = Fsum;
Dmean = Psum;
Pmgr = Psum;
Fmgr = Fsum;
Dmgr = Psum;
Pprod = Psum;
Fprod = Fsum;
Dprod = Psum;
Pcon = Psum;
Fcon = Fsum;
Dcon = Psum;
Prep = Fsum;
Frep = Fsum;
Drep = Fsum;
all_mean=NaN*ones(3,3,length(spots));
z = NaN*ones(length(spots),3);

%%
for s=1:length(spots)
    %%
    loc = spots{s};
    lname = [sname2 loc '_'];
    SP = csvread([dpath sname lname 'Sml_p.csv']);
    SF = csvread([dpath sname lname 'Sml_f.csv']);
    SD = csvread([dpath sname lname 'Sml_d.csv']);
    MP = csvread([dpath sname lname 'Med_p.csv']);
    MF = csvread([dpath sname lname 'Med_f.csv']);
    MD = csvread([dpath sname lname 'Med_d.csv']);
    LP = csvread([dpath sname lname 'Lrg_p.csv']);
    LD = csvread([dpath sname lname 'Lrg_d.csv']);
    C = csvread([dpath sname lname 'Cobalt.csv']);
    
    %% Plots over time
    x=1:length(SP);
    y=x/365;
    lstd=length(SP);
    id1 = 0:365:(lstd-1);
    id2 = 365:365:(lstd);
    ID  = [id1 id2];
    
    
    %% Final mean biomass in each size
    t=1:length(SP);
    lyr=t((end-365+1):end);
    
    SP_sum=sum(SP(lyr,1));
    SF_sum=sum(SF(lyr,1));
    SD_sum=sum(SD(lyr,1));
    MP_sum=sum(MP(lyr,1));
    MF_sum=sum(MF(lyr,1));
    MD_sum=sum(MD(lyr,1));
    LP_sum=sum(LP(lyr,1));
    LD_sum=sum(LD(lyr,1));
    
    SP_mean=mean(SP(lyr,1));
    SF_mean=mean(SF(lyr,1));
    SD_mean=mean(SD(lyr,1));
    MP_mean=mean(MP(lyr,1));
    MF_mean=mean(MF(lyr,1));
    MD_mean=mean(MD(lyr,1));
    LP_mean=mean(LP(lyr,1));
    LD_mean=mean(LD(lyr,1));
    
    P_sum=[SP_sum;MP_sum;LP_sum];
    F_sum=[SF_sum;MF_sum];
    D_sum=[SD_sum;MD_sum;LD_sum];
    P_mean=[SP_mean;MP_mean;LP_mean];
    F_mean=[SF_mean;MF_mean];
    D_mean=[SD_mean;MD_mean;LD_mean];
    
    Psum(:,s) = P_sum;
    Fsum(:,s) = F_sum;
    Dsum(:,s) = D_sum;
    Pmean(:,s) = P_mean;
    Fmean(:,s) = F_mean;
    Dmean(:,s) = D_mean;
    
    
    all_mean(1:2,1,s) = F_mean;
    all_mean(:,2,s) = P_mean;
    all_mean(:,3,s) = D_mean;
    
    f1 = figure(1);
    subplot(3,3,s)
    plot(0.5:2:5.5,log10(squeeze(all_mean(:,1,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,log10(squeeze(all_mean(:,2,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.5:2:6.5,log10(squeeze(all_mean(:,3,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 6])
    %ylim([-20 10])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    if (s==4)
        %legend('F','P','D')
        %legend('location','southeast')
        ylabel('log10 Mean Biom (g m^-^2) in final year')
    end
    title(loc)
    xlabel('Stage')
    if (s==3)
        stamp(cfile)
    end
    
    f21 = figure(21);
    subplot(3,3,s)
    plot(0.5:2:5.5,log10(squeeze(all_mean(:,1,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,log10(squeeze(all_mean(:,2,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.5:2:6.5,log10(squeeze(all_mean(:,3,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 6])
    ylim([-150 2])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    if (s==4)
        ylabel('log10 Mean Biom (g m^-^2) in final year')
    end
    title(loc)
    xlabel('Stage')
    if (s==3)
        stamp(cfile)
    end
    
    %% Feeding level
    f2=figure(2);
    subplot(3,3,s)
    bar(mclev(s,:),'k')
    ylim([0 1])
    set(gca,'XTickLabel',[]);
    for n=1:8
        text(n-0.5,-0.2,stage{n},'Rotation',45)
    end
    title(spots{s})
    if (s==4)
        ylabel('Feeding level')
    end
    stamp(cfile)
    
    %% Growth rate (nu - energy for biomass production)
    SP_mgr=nanmean(SP(lyr,15));
    SF_mgr=nanmean(SF(lyr,15));
    SD_mgr=nanmean(SD(lyr,15));
    MP_mgr=nanmean(MP(lyr,15));
    MF_mgr=nanmean(MF(lyr,15));
    MD_mgr=nanmean(MD(lyr,15));
    LP_mgr=nanmean(LP(lyr,15));
    LD_mgr=nanmean(LD(lyr,15));
    
    P_mgr=[SP_mgr;MP_mgr;LP_mgr];
    F_mgr=[SF_mgr;MF_mgr];
    D_mgr=[SD_mgr;MD_mgr;LD_mgr];
    
    Pmgr(:,s) = P_mgr;
    Fmgr(:,s) = F_mgr;
    Dmgr(:,s) = D_mgr;
    
    f3 = figure(3);
    subplot(1,3,1)
    plot(s-0.25,F_mgr(1),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,P_mgr(1),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,D_mgr(1),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-0.1 0.1])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha1=gca;
        for n=1:7
            text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean growth rate (g g^-^1 d^-^1) in final year')
    title('S')
    
    subplot(1,3,2)
    plot(s-0.25,(F_mgr(2)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(P_mgr(2)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(D_mgr(2)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-2 5])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha2=gca;
        for n=1:7
            text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean growth/repro rate (g g^-^1 d^-^1) in final year')
    title('M')
    
    subplot(1,3,3)
    plot(s,(P_mgr(3)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(D_mgr(3)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-2 7])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha3=gca;
        for n=1:7
            text(n-0.5,ha3.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean repro rate (g g^-^1 d^-^1) in final year')
    title('L')
    if (s==3)
        stamp(cfile)
    end
    
    %% Consump per biomass (I)
    SP_con=nanmean(SP(lyr,14));
    SF_con=nanmean(SF(lyr,14));
    SD_con=nanmean(SD(lyr,14));
    MP_con=nanmean(MP(lyr,14));
    MF_con=nanmean(MF(lyr,14));
    MD_con=nanmean(MD(lyr,14));
    LP_con=nanmean(LP(lyr,14));
    LD_con=nanmean(LD(lyr,14));
    
    P_con=[SP_con;MP_con;LP_con];
    F_con=[SF_con;MF_con];
    D_con=[SD_con;MD_con;LD_con];
    
    Pcon(:,s) = P_con;
    Fcon(:,s) = F_con;
    Dcon(:,s) = D_con;
    
    f4 = figure(4);
    subplot(3,3,s)
    plot(0.5:2:3.5,(F_con),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,(P_con),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.5:2:6.5,(D_con),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 6])
    %ylim([-25 15])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    if (s==4)
        %legend('F','P','D')
        %legend('location','southeast')
        ylabel('Mean consumption (g g^-^1 d^-^1) in final year')
    end
    title(loc)
    xlabel('Stage')
    if (s==3)
        stamp(cfile)
    end
    
    %% Fraction zoop losses consumed
    z(s,1) = nanmean(C(lyr,2));
    z(s,2) = nanmean(C(lyr,3));
    z(s,3) = nanmean(C(lyr,4));
    
    f5 = figure(5);
    subplot(3,3,s)
    bar(z(s,:)); hold on;
    xlim([0 4])
    ylim([0 1])
    set(gca,'XTick',1:3,'XTickLabel',{'MZ','LZ','Det'})
    if (s==4)
        ylabel('Fraction flux consumed')
    end
    title(loc)
    if (s==3)
        stamp(cfile)
    end
    
    %% Size spectrum (sum stages)
    spec = nansum(all_mean(:,:,s),2);
    
    f6 = figure(6);
    subplot(3,3,s)
    plot(1:2:6,log10(spec),'sk',...
        'MarkerFaceColor','k',...
        'MarkerSize',15); hold on;
    xlim([0 6])
    %ylim([-4 4])clo
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    title(loc)
    if (s==4)
        ylabel('log Mean Biom (g m^-^2) in final year')
    end
    xlabel('Size')
    if (s==3)
        stamp(cfile)
    end
    
    f7 = figure(7);
    stamp(cfile)
    %     plot(1:2:6,log10(spec),'color',cmap_ppt(s,:),...
    %         'LineWidth',2); hold on;
    plot(1:2:6,log10(spec),'LineWidth',2); hold on;
    xlim([0 6])
    %ylim([-25 15])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    if (s==7)
        legend(spots)
        legend('location','northwest')
    end
    ylabel('log Mean Biom (g m^-^2) in final year')
    xlabel('Size class')
    if (s==1)
        stamp(cfile)
    end
    
    %% Production (= nu * biom)
    SP_prod=nanmean(SP(lyr,24));
    SF_prod=nanmean(SF(lyr,24));
    SD_prod=nanmean(SD(lyr,24));
    MP_prod=nanmean(MP(lyr,24));
    MF_prod=nanmean(MF(lyr,24));
    MD_prod=nanmean(MD(lyr,24));
    LP_prod=nanmean(LP(lyr,24));
    LD_prod=nanmean(LD(lyr,24));
    
    P_prod=[SP_prod;MP_prod;LP_prod];
    F_prod=[SF_prod;MF_prod];
    D_prod=[SD_prod;MD_prod;LD_prod];
    
    Pprod(:,s) = P_prod;
    Fprod(:,s) = F_prod;
    Dprod(:,s) = D_prod;
    
    f8 = figure(8);
    subplot(1,3,1)
    plot(s-0.25,F_prod(1),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,P_prod(1),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,D_prod(1),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-0.1 0.1])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha1=gca;
        for n=1:7
            text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
    title('S')
    
    subplot(1,3,2)
    plot(s-0.25,(F_prod(2)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(P_prod(2)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(D_prod(2)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-2 5])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha2=gca;
        for n=1:7
            text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
    title('M')
    
    subplot(1,3,3)
    plot(s,(P_prod(3)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(D_prod(3)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-2 7])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha3=gca;
        for n=1:7
            text(n-0.5,ha3.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
    title('L')
    if (s==1)
        stamp(cfile)
    end
    
    %% Reproduction
    F_rep(1,1)=nanmean(MF(lyr,18));
    D_rep(1,1)=nanmean(LD(lyr,18));
    P_rep(1,1)=nanmean(LP(lyr,18));
    F_rep(2,1)=nanmean(MF(lyr,1).*MF(lyr,18));
    D_rep(2,1)=nanmean(LD(lyr,1).*LD(lyr,18));
    P_rep(2,1)=nanmean(LP(lyr,1).*LP(lyr,18));
    
    Prep(:,s) = P_rep;
    Frep(:,s) = F_rep;
    Drep(:,s) = D_rep;
    
    f9 = figure(9);
    subplot(1,2,1)
    plot(s-0.25,F_rep(1),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,P_rep(1),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,D_rep(1),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-0.1 0.1])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha1=gca;
        for n=1:7
            text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean repro rate (g g^-^1 d^-^1) in final year')
    
    subplot(1,2,2)
    plot(s-0.25,(F_rep(2)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(P_rep(2)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(D_rep(2)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-2 5])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha2=gca;
        for n=1:7
            text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean biom reproduced (g d^-^1) in final year')
    if (s==1)
        stamp(cfile)
    end
    
    
end
print(f1,'-dpng',[fpath sname sname2 'All_oneloc_Logmean_biomass.png'])
print(f21,'-dpng',[fpath sname sname2 'All_oneloc_Logmean_biomass_axes.png'])
print(f2,'-dpng',[fpath sname sname2 'All_oneloc_con_level.png'])
print(f3,'-dpng',[fpath sname sname2 'All_oneloc_nu.png'])
print(f4,'-dpng',[fpath sname sname2 'All_oneloc_consump.png'])
print(f5,'-dpng',[fpath sname sname2 'All_oneloc_frac_zoop_loss.png'])
print(f6,'-dpng',[fpath sname sname2 'All_oneloc_size_spec_sub.png'])
print(f7,'-dpng',[fpath sname sname2 'All_oneloc_size_spec.png'])
print(f8,'-dpng',[fpath sname sname2 'All_oneloc_prod.png'])
print(f9,'-dpng',[fpath sname sname2 'All_oneloc_rep.png'])

save([dpath sname sname2 'lastyr_sum_mean_biom'],'Psum','Fsum',...
    'Dsum','Pmean','Fmean','Dmean','all_mean',...
    'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
    'Prep','Frep','Drep');

%% Sum mean biom over stages
sumspec = squeeze(nansum(nansum(all_mean)));

figure(10);
subplot(2,1,1)
plot(1:7,log10(sumspec),'k.','MarkerSize',25); hold on;
xlim([0 8])
%ylim([-2 1])
set(gca,'XTick',1:7,'XTickLabel',[])
for n=1:7
    text(n,-0.6,spots{n},'HorizontalAlignment','center')
end
ylabel('log10 Mean Biom (g m^-^2) in final year')
title('All fishes and stages')

subplot(2,1,2)
plot(1:7,(sumspec),'k.','MarkerSize',25); hold on;
xlim([0 8])
set(gca,'XTick',1:7,'XTickLabel',[])
for n=1:7
    text(n,-1,spots{n},'HorizontalAlignment','center')
end
ylabel('Mean Biom (g m^-^2) in final year')
stamp(cfile)
print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_spec.png'])

%% All on one
figure(11)
subplot(2,2,1)
bar(log(Fmean))
xlim([0 4])
ylim([-10 3])
xlabel('Stage')
title('Forage')
ylabel('log Mean Biomass (g m^-^2) in final year',...
    'HorizontalAlignment','right')

subplot(2,2,2)
bar(log(Pmean))
xlim([0 4])
ylim([-10 3])
xlabel('Stage')
title('Pel Pisc')

subplot(2,2,3)
bar(log(Dmean))
xlim([0 4])
ylim([-10 3])
xlabel('Stage')
title('Dem Pisc')

subplot(2,2,4)
bar(log(Fsum))
xlim([0 0.5])
ylim([-0.5 0.5])
legend(spots)
legend('location','west')
stamp(cfile)
print('-dpng',[fpath sname sname2 'All_oneloc_biomass_spec.png'])






