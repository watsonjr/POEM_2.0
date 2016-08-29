%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

% dpath = [datap 'NoPDc_NoAct_TrefO_flev1e4/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_flev4e4/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_flev8e4/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e4_NoWgt/'];
% dpath = [datap 'NoPDc_NoMetab_TrefO_1e4/'];
% dpath = [datap 'NoPDc_NoMetab_TrefO_1e4_HalfC/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e4_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_4e4_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_8e4_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e5_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e6_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_Cmax_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoMet_TrefO_Cmax_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e6_noC&Mwgt/'];
% dpath = [datap 'NoPDc_TrefO_KHparams_cmax-metab/'];
% dpath = [datap 'NoPDc_TrefO_KHparams_cmax-metab_MFeatS/'];
% dpath = [datap 'NoPDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MprefLZoverMZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFprefZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit01/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFmet/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFpred/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/'];
dpath1 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LTmort/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_simpQmort/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_compQmort/'];
% dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_bioQmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF50/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF75/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_sameA/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_Fenc2x_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Lmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Lmort/'];
% dpath1 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit15/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit20/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit30/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit40/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05_Tmort/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc15/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc20/'];
% dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc25/'];
% dpath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc30/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc075/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc050/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann_Tmort/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann_Tmort/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20/'];
% dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30/'];
% dpath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40/'];
% dpath9 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010_NOnmort/'];
% dpath10 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025_NOnmort/'];
% dpath11 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05_NOnmort/'];
% dpath12 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10_NOnmort/'];
% dpath13 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20_NOnmort/'];
% dpath14 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30_NOnmort/'];
% dpath15 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40_NOnmort/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ00/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ05/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ00/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ01/'];
% dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ05/'];
dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_05/'];
dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_075/'];
dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_125/'];
dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_15/'];

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

% dp = {dpath1; dpath2; dpath3; dpath4; dpath5; dpath6};
% sims = {'f05','f10','f15','f20','f30','f40'};
% cfile = 'MFeqMP_fcrit_comp';
% dp = {dpath3;dpath4;dpath2;dpath5;dpath6;dpath7;dpath8};
% sims = {'MP0.50','MP0.75','Meq','MF1.5','MF2.0','MF2.5','MF3.0'};
% cfile = 'MFbetter_enc_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5;dpath6;dpath7};
% sims = {'con','T','L','T&L','simpQ','compQ1','compQ1'};
% cfile = 'MFeqMP_mort_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5};
% sims = {'full','redSM','T-redSM','redS','T-redS'};
% cfile = 'MFeqMP_cannibal_comp';
% dp = {dpath9;dpath10;dpath11;dpath12;dpath13;dpath14;dpath15;dpath1;dpath4;...
%     dpath5;dpath6;dpath7;dpath8};
% sims = {'0.01','0.025','0.05','0.10','0.20','0.30','0.40','0.0+N','0.01+N',...
%     '0.025+N','0.05+N','0.10+N','0.20+N','0.30+N'};
% cfile = 'MFeqMP_fishing_comp';
% dp = {dpath2;dpath3;dpath4;dpath1;dpath5;dpath6;dpath7};
% sims = {'0.0','0.1','0.5','1.0','MP0.0','MP0.1','MP0.5'};
% cfile = 'MFeqMP_MZ_comp';
dp = {dpath2;dpath3;dpath1;dpath4;dpath5};
sims = {'0.5','0.75','1.0','1.25','1.5'};
cfile = 'MFeqMP_Mcann_comp';

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

load('cmap_ppt_angles.mat')

%%
for i=1:length(spots)
    close all
        
    loc = spots{i};
    lname = [sname2 loc '_'];
    
    %%
    for s=1:length(dp)
        
        dpath = char(dp(s));
        load([dpath sname sname2 'consump.mat'],'mclev','Zcon');
        load([dpath sname sname2 'lastyr_sum_mean_biom']);
        
        %%
        f1 = figure(1);
        subplot(3,1,1)
        plot(s-0.25,log10(squeeze(all_mean(1,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(1,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(1,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,log10(squeeze(all_mean(2,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(2,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(2,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title('M')
        
        subplot(3,1,3)
        plot(s,log10(squeeze(all_mean(3,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(3,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title('L')
        xlabel('Sim')
        
        
        %% Feeding level
        f2=figure(2);
        subplot(3,3,1)
        plot(s,mclev(i,1),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title(stage(1))
        
        subplot(3,3,2)
        plot(s,mclev(i,2),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title([loc ' ' stage(2)])
        
        subplot(3,3,3)
        plot(s,mclev(i,3),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            stamp(cfile)
        end
        title(stage(3))
        
        subplot(3,3,4)
        plot(s,mclev(i,4),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('Feeding level')
        end
        title(stage(4))
        
        subplot(3,3,5)
        plot(s,mclev(i,5),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title(stage(5))
        
        subplot(3,3,6)
        plot(s,mclev(i,6),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title(stage(6))
        
        subplot(3,3,7)
        plot(s,mclev(i,7),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title(stage(7))
        
        subplot(3,3,8)
        plot(s,mclev(i,8),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title(stage(8))
        
        %% Growth rate (nu - energy for biomass production)
        f3 = figure(3);
        subplot(3,1,1)
        plot(s-0.25,Fmgr(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Pmgr(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dmgr(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,(Fmgr(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Pmgr(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmgr(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        ylabel('Mean growth/repro rate (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(3,1,3)
        plot(s,(Pmgr(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmgr(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title('L')
        
        %% Consump per biomass (I)
        f4 = figure(4);
        subplot(3,1,1)
        plot(s-0.25,Fcon(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Pcon(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dcon(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,(Fcon(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Pcon(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dcon(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        ylabel('Mean consumption rate (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(3,1,3)
        plot(s,(Pcon(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dcon(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title('L')
        
        %% Fraction zoop losses consumed
        f5 = figure(5);
        subplot(3,1,1)
        plot(s,z(i,1),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' Med zoo'])
        
        subplot(3,1,2)
        plot(s,z(i,2),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title('Large zoo')
        ylabel('Fraction of flux consumed')
        
        subplot(3,1,3)
        plot(s,z(i,3),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title('Detritus')
        
        
        %% Zoop con
        f12 = figure(12);
        subplot(2,1,1)
        plot(s,Zcon(i,1),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' Med zoo'])
        ylabel('Fraction of times overconsumed')
        
        subplot(2,1,2)
        plot(s,Zcon(i,2),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title([loc ' Large zoo'])
        ylabel('Fraction of times overconsumed')
        
        
        %% Size spectrum (sum stages)
        spec = nansum(all_mean(:,:,i),2);
        
        f7 = figure(7);
        stamp(cfile)
        plot(1:2:6,log10(spec),'LineWidth',2); hold on;
        xlim([0 6])
        set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
        if (s==length(dp))
            legend(sims)
            legend('location','northwest')
            stamp(cfile)
        end
        ylabel('log Mean Biom (g m^-^2) in final year')
        xlabel('Size class')
        title(loc)
        
        %% Production (= nu * biom)
        f8 = figure(8);
        subplot(3,1,1)
        plot(s-0.25,Fprod(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Pprod(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dprod(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,(Fprod(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Pprod(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dprod(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(3,1,3)
        plot(s,(Pprod(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dprod(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        title('L')
        
        %% Reproduction
        f9 = figure(9);
        subplot(2,1,1)
        plot(s-0.25,Frep(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Prep(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Drep(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            stamp(cfile)
        end
        ylabel('Mean repro rate (g g^-^1 d^-^1)')
        title(loc)
        
        subplot(2,1,2)
        plot(s-0.25,(Frep(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Prep(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Drep(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        ylabel('Mean biom reproduced (g d^-^1)')
        
        %% Metabolism
%         f10 = figure(10);
%         subplot(3,1,1)
%         plot(s-0.25,Fmet(1,i),'sk',...
%             'MarkerFaceColor',cmap_ppt(3,:),...
%             'MarkerSize',15); hold on;
%         plot(s,Pmet(1,i),'sk',...
%             'MarkerFaceColor',cmap_ppt(1,:),...
%             'MarkerSize',15); hold on;
%         plot(s+0.25,Dmet(1,i),'sk',...
%             'MarkerFaceColor',cmap_ppt(2,:),...
%             'MarkerSize',15); hold on;
%         xlim([0 length(dp)+1])
%         if (s==length(dp))
%             set(gca,'XTick',1:length(dp),'XTickLabel',sims);
%             stamp(cfile)
%         end
%         title([loc ' S'])
%         
%         subplot(3,1,2)
%         plot(s-0.25,(Fmet(2,i)),'sk',...
%             'MarkerFaceColor',cmap_ppt(3,:),...
%             'MarkerSize',15); hold on;
%         plot(s,(Pmet(2,i)),'sk',...
%             'MarkerFaceColor',cmap_ppt(1,:),...
%             'MarkerSize',15); hold on;
%         plot(s+0.25,(Dmet(2,i)),'sk',...
%             'MarkerFaceColor',cmap_ppt(2,:),...
%             'MarkerSize',15); hold on;
%         xlim([0 length(dp)+1])
%         if (s==length(dp))
%             set(gca,'XTick',1:length(dp),'XTickLabel',sims);
%         end
%         ylabel('Mean metabolic rate (g g^-^1 d^-^1) in final year')
%         title('M')
%         
%         subplot(3,1,3)
%         plot(s,(Pmet(3,i)),'sk',...
%             'MarkerFaceColor',cmap_ppt(1,:),...
%             'MarkerSize',15); hold on;
%         plot(s+0.25,(Dmet(3,i)),'sk',...
%             'MarkerFaceColor',cmap_ppt(2,:),...
%             'MarkerSize',15); hold on;
%         xlim([0 length(dp)+1])
%         if (s==length(dp))
%             set(gca,'XTick',1:length(dp),'XTickLabel',sims);
%         end
%         title('L')
        
        %% Predation
%         f11 = figure(11);
%         subplot(2,1,1)
%         plot(s-0.25,Fpred(1,i),'sk',...
%             'MarkerFaceColor',cmap_ppt(3,:),...
%             'MarkerSize',15); hold on;
%         plot(s,Ppred(1,i),'sk',...
%             'MarkerFaceColor',cmap_ppt(1,:),...
%             'MarkerSize',15); hold on;
%         plot(s+0.25,Dpred(1,i),'sk',...
%             'MarkerFaceColor',cmap_ppt(2,:),...
%             'MarkerSize',15); hold on;
%         xlim([0 length(dp)+1])
%         if (s==length(dp))
%             set(gca,'XTick',1:length(dp),'XTickLabel',sims);
%             stamp(cfile)
%         end
%         title([loc ' S'])
%         
%         subplot(2,1,2)
%         plot(s-0.25,(Fpred(2,i)),'sk',...
%             'MarkerFaceColor',cmap_ppt(3,:),...
%             'MarkerSize',15); hold on;
%         plot(s,(Ppred(2,i)),'sk',...
%             'MarkerFaceColor',cmap_ppt(1,:),...
%             'MarkerSize',15); hold on;
%         plot(s+0.25,(Dpred(2,i)),'sk',...
%             'MarkerFaceColor',cmap_ppt(2,:),...
%             'MarkerSize',15); hold on;
%         xlim([0 length(dp)+1])
%         if (s==length(dp))
%             set(gca,'XTick',1:length(dp),'XTickLabel',sims);
%         end
%         ylabel('Mean predation rate (g g^-^1 d^-^1) in final year',...
%             'HorizontalAlignment','left')
        title('M')
        
        %% Sum mean biom over stages
        sumspec = squeeze(nansum(nansum(all_mean)));
        
        f15=figure(15);
        plot(s,log10(sumspec(i)),'k.','MarkerSize',25); hold on;
        xlim([0 length(dp)+1])
        if (s==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title([loc ' All fishes and stages'])
        
    end
    
    print(f1,'-dpng',[fpath sname sname2 cfile '_' lname 'Logmean_biomass.png'])
    print(f2,'-dpng',[fpath sname sname2 cfile '_' lname 'con_level.png'])
    print(f3,'-dpng',[fpath sname sname2 cfile '_' lname 'nu.png'])
    print(f4,'-dpng',[fpath sname sname2 cfile '_' lname 'consump.png'])
    print(f5,'-dpng',[fpath sname sname2 cfile '_' lname 'frac_zoop_loss.png'])
    print(f7,'-dpng',[fpath sname sname2 cfile '_' lname 'size_spec.png'])
    print(f8,'-dpng',[fpath sname sname2 cfile '_' lname 'prod.png'])
    print(f9,'-dpng',[fpath sname sname2 cfile '_' lname 'rep.png'])
    %print(f10,'-dpng',[fpath sname sname2 cfile '_' lname 'met.png'])
    %print(f11,'-dpng',[fpath sname sname2 cfile '_' lname 'pred.png'])
    print(f12,'-dpng',[fpath sname sname2 cfile '_' lname 'zoo_con.png'])
    print(f15,'-dpng',[fpath sname sname2 cfile '_' lname 'tot_mean_biomass_spec.png'])
        
end




