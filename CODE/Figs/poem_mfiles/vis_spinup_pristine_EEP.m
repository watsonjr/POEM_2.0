% Visualize output of POEM
% Pristine historical at one location
% 145 years

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort/';
%
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
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_sameA/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_sameA/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_MFdiffA1/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_MFdiffA1/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_MFdiffA2/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_MFdiffA2/';

cfile = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP'};

cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S'};
cols=cols';

%s=7;
for s=7%1:length(spots)
    %%
    close all
    loc = spots{s};
    sname = 'Spinup_';
    lname = [loc '_'];
    %lname = ['phen_' loc '_'];
    SP = csvread([dpath sname lname 'Sml_p.csv']);
    SF = csvread([dpath sname lname 'Sml_f.csv']);
    SD = csvread([dpath sname lname 'Sml_d.csv']);
    MP = csvread([dpath sname lname 'Med_p.csv']);
    MF = csvread([dpath sname lname 'Med_f.csv']);
    MD = csvread([dpath sname lname 'Med_d.csv']);
    LP = csvread([dpath sname lname 'Lrg_p.csv']);
    LD = csvread([dpath sname lname 'Lrg_d.csv']);
    C = csvread([dpath sname lname 'Cobalt.csv']);
    z(:,1) = C(:,2);
    z(:,2) = C(:,3);
    z(:,3) = C(:,4);
    z=floor(z);
    
    %% Plots over time
    x=1:length(SP);
    y=x(1:10950)/365;
    lstd=length(SP);
    
    
    %% Mean consumption level
    c=[SF(:,21) SP(:,21) SD(:,21) MF(:,21) MP(:,21) MD(:,21) LP(:,21) LD(:,21)];
    mclev(s,:) = nanmean(c);
    
    % Zoop overconsumption
    
    Zcon(s,:) = nansum(z)/lstd;
    
    %% PLOTS
    
    %% Piscivore
    figure(1)
    stamp(cfile)
    subplot(4,1,1)
    plot(y,log10(SP(1:6570,1)),'b','Linewidth',1); hold on;
    plot(y,log10(MP(1:6570,1)),'r','Linewidth',1); hold on;
    plot(y,log10(LP(1:6570,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Pelagic Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    legend('Larvae','Juveniles','Adults')
    
    subplot(4,1,2)
    plot(y,log10(SP(1:6570,1)),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y,log10(MP(1:6570,1)),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y,log10(LP(1:6570,1)),'k','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_pisc_time_30yr.png'])
    
    %% Planktivore
    figure(2)
    subplot(3,1,1)
    plot(y,log10(SF(1:6570,1)),'b','Linewidth',1); hold on;
    plot(y,log10(MF(1:6570,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Forage Fishes ' loc])
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    legend('Immature','Adults')
    stamp(cfile)
    
    subplot(3,1,2)
    plot(y,log10(SF(1:6570,1)),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Immature')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(3,1,3)
    plot(y,log10(MF(1:6570,1)),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    print('-dpng',[fpath sname lname 'oneloc_plan_time_30yr.png'])
    
    %% Detritivore
    figure(3)
    subplot(4,1,1)
    plot(y,log10(SD(1:6570,1)),'b','Linewidth',1); hold on;
    plot(y,log10(MD(1:6570,1)),'r','Linewidth',1); hold on;
    plot(y,log10(LD(1:6570,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Demersal Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    legend('Larvae','Juveniles','Adults')
    stamp(cfile)
    
    subplot(4,1,2)
    plot(y,log10(SD(1:6570,1)),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y,log10(MD(1:6570,1)),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y,log10(LD(1:6570,1)),'k','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_detr_time_30yr.png'])
    
    %% All biomass in subplots
    %SP
    figure(4)
    subplot(3,3,2)
    plot(y,log10(SP(1:6570,1)),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title({loc; 'SP'})
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    stamp(cfile)
    
    subplot(3,3,5)
    plot(y,log10(MP(1:6570,1)),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('MP')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(3,3,8)
    plot(y,log10(LP(1:6570,1)),'k','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('LP')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    %FF
    subplot(3,3,1)
    plot(y,log10(SF(1:6570,1)),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('SF')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(3,3,4)
    plot(y,log10(MF(1:6570,1)),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('MF')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    %Detritivore
    subplot(3,3,3)
    plot(y,log10(SD(1:6570,1)),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('SD')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    stamp(cfile)
    
    subplot(3,3,6)
    plot(y,log10(MD(1:6570,1)),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('MD')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(3,3,9)
    plot(y,log10(LD(1:6570,1)),'k','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('LD')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_sub_30yr.png'])
    
    %% All size classes of all
    
    figure(5)
    plot(y,log10(SP(1:6570,1)),'Linewidth',2); hold on;
    plot(y,log10(MP(1:6570,1)),'Linewidth',2); hold on;
    plot(y,log10(LP(1:6570,1)),'Linewidth',2); hold on;
    plot(y,log10(SF(1:6570,1)),'Linewidth',2); hold on;
    plot(y,log10(MF(1:6570,1)),'Linewidth',2); hold on;
    plot(y,log10(SD(1:6570,1)),'Linewidth',2); hold on;
    plot(y,log10(MD(1:6570,1)),'Linewidth',2); hold on;
    plot(y,log10(LD(1:6570,1)),'Linewidth',2); hold on;
    legend('SP','MP','LP','SF','MF','SD','MD','LD')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    title(['Spinup ' loc])
    stamp(cfile)
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_30yr.png'])
    
    %F starts to decrease at y=0.1699, d=62
    
    %% Reproduction
    rep(1:6570,1)=MF(1:6570,1).*MF(1:6570,30);
    rep(1:6570,2)=LD(1:6570,1).*LD(1:6570,30);
    rep(1:6570,3)=LP(1:6570,1).*LP(1:6570,30);
    
    figure(7)
    subplot(4,1,1)
    plot(y,log10(rep),'Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Reproduction ' loc])
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    legend('F','D','P')
    stamp(cfile)
    
    subplot(4,1,2)
    plot(y,log10(rep(1:6570,1)),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y,log10(rep(1:6570,2)),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Demersal Piscivores')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y,log10(rep(1:6570,3)),'k','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Pelagic Piscivores')
    xlabel('Time (y)')
    ylabel('log10 Biomass (g m^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_rep_time_30yr.png'])
    
    %% Maturation
    m(:,1)=MF(:,19);
    m(:,2)=MD(:,19);
    m(:,3)=MP(:,19);
    m(:,4)=LD(:,19);
    m(:,5)=LP(:,19);
    
    figure(8)
    subplot(3,2,1)
    plot(y,log10(m(1:6570,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title([loc ' log10 Maturation Biomass (g m^-^2)'],'HorizontalAlignment','left')
    ylabel('Forage Fishes')
    stamp(cfile)
    
    subplot(3,2,3)
    plot(y,log10(m(1:6570,2)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('M')
    ylabel('Demersal Piscivores')
    
    subplot(3,2,5)
    plot(y,log10(m(1:6570,3)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('M')
    ylabel('Pelagic Piscivores')
    xlabel('Time (y)')
    
    subplot(3,2,4)
    plot(y,log10(m(1:6570,4)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('L')
    ylabel('Demersal Piscivores')
    
    subplot(3,2,6)
    plot(y,log10(m(1:6570,5)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('L')
    xlabel('Time (y)')
    ylabel('Pelagic Piscivores')
    print('-dpng',[fpath sname lname 'oneloc_matur_time_30yr.png'])
    
    %% Predation mortality
    %SP
    figure(9)
    subplot(2,3,2)
    plot(y,log10(SP(1:6570,17)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title({loc; 'SP'})
    stamp(cfile)
    
    subplot(2,3,5)
    plot(y,log10(MP(1:6570,17)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MP')
    xlabel('Time (y)')
    
    %FF
    subplot(2,3,1)
    plot(y,log10(SF(1:6570,17)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SF')
    ylabel('log10 Biomass eaten by predators (g m^-^2)','HorizontalAlignment','right')
    
    subplot(2,3,4)
    plot(y,log10(MF(1:6570,17)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MF')
    xlabel('Time (y)')
    %ylabel('Biomass eaten by predators (g m^-^2)')
    
    %Detritivore
    subplot(2,3,3)
    plot(y,log10(SD(1:6570,17)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SD')
    
    subplot(2,3,6)
    plot(y,log10(MD(1:6570,17)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MD')
    xlabel('Time (y)')
    
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_pred_sub_30yr.png'])
    
    %% All consumption in subplots
    %SP
    figure(10)
    subplot(3,3,2)
    plot(y,(SP(1:6570,14)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([0.2 0.6])
    title({loc; 'SP'})
    stamp(cfile)
    
    subplot(3,3,5)
    plot(y,(MP(1:6570,14)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([0.015 0.05])
    title('MP')
    
    subplot(3,3,8)
    plot(y,(LP(1:6570,14)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LP')
    xlabel('Time (y)')
    
    %FF
    subplot(3,3,1)
    plot(y,(SF(1:6570,14)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([0.2 0.6])
    title('SF')
    
    subplot(3,3,4)
    plot(y,(MF(1:6570,14)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([0.015 0.05])
    title('MF')
    xlabel('Time (y)')
    ylabel('Biomass Consumed (g g^-^1 m^-^2)')
    
    %Detritivore
    subplot(3,3,3)
    plot(y,(SD(1:6570,14)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([0.2 0.6])
    title('SD')
    
    subplot(3,3,6)
    plot(y,(MD(1:6570,14)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    ylim([0.015 0.05])
    title('MD')
    
    subplot(3,3,9)
    plot(y,(LD(1:6570,14)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LD')
    xlabel('Time (y)')
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_consump_sub_30yr.png'])
    
    %% Consumption per biomass
    %SP
    figure(11)
    subplot(3,3,2)
    plot(y,(SP(1:6570,8:10)),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    legend('F','P','D')
    title({loc; 'SP'})
    stamp(cfile)
    
    subplot(3,3,5)
    plot(y,(MP(1:6570,8:10)),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MP')
    
    subplot(3,3,8)
    plot(y,(LP(1:6570,8:10)),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LP')
    xlabel('Time (y)')
    
    %FF
    subplot(3,3,1)
    plot(y,(SF(1:6570,8:10)),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SF')
    
    subplot(3,3,4)
    plot(y,(MF(1:6570,8:10)),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MF')
    xlabel('Time (y)')
    ylabel('Biomass Consumed (g g^-^1 m^-^2)')
    
    %Detritivore
    subplot(3,3,3)
    plot(y,(SD(1:6570,8:10)),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SD')
    
    subplot(3,3,6)
    plot(y,(MD(1:6570,8:10)),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MD')
    
    subplot(3,3,9)
    plot(y,(LD(1:6570,8:10)),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LD')
    xlabel('Time (y)')
    print('-dpng',[fpath sname lname 'oneloc_all_sizes_conFPD_sub_30yr.png'])
end

