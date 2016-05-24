% Compare continuous vs spawning season
% Spinup pristine historical at one location
% 100 years

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ_TrefOrig/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ_TrefOrig/';

load([dpath 'Spinup_oneloc_phenol_all.mat'])
load([dpath 'Spinup_oneloc_all.mat'])
load('cmap_ppt_angles.mat')

%% Plots over time
x=1:length(psp);
lyr=x((end-365+1):end);

%% Degree days
for s=1:length(spots)
    loc = spots{s};
    figure(1)
    clf
    plot(1:365,pDDlp(lyr,s),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    plot(1:365,pDDld(lyr,s),'Linewidth',2,'color',cmap_ppt(2,:)); hold on;
    xlim([0 365])
    title(['Spinup 1980 ' loc])
    xlabel('Time (y)')
    ylabel('Cumulative degree days')
    legend('Pelagic','Demersal')
    print('-dpng',[fpath loc '_oneloc_DD_2000.png'])
end

%% Spawn flag
for s=1:length(spots)
    loc = spots{s};
    figure(2)
    clf
    plot(1:365,pKlp(lyr,s),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    plot(1:365,pKld(lyr,s),'Linewidth',2,'color',cmap_ppt(2,:)); hold on;
    xlim([0 365])
    %ylim([-0.1 1.1])
    title(['Spinup 1980 ' loc])
    xlabel('Time (y)')
    ylabel('Spawning flag')
    legend('Pelagic','Demersal')
    print('-dpng',[fpath loc '_oneloc_K_2000.png'])
end

%% Reprod
for s=1:length(spots)
    loc = spots{s};
    figure(3)
    clf
    subplot(3,1,1)
    plot(1:365,1000*pKmf(lyr,s),'--k'); hold on;
    plot(1:365,pRmf(lyr,s),'Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    xlim([0 365])
    ylim([0 max(pRmf(lyr,s))])
    title({[loc ' Spinup 1980']; ['Forage Fishes']})
    
    subplot(3,1,3)
    plot(1:365,1000*pKld(lyr,s),'--k'); hold on;
    plot(1:365,pRld(lyr,s),'Linewidth',2,'color',cmap_ppt(2,:)); hold on;
    xlim([0 365])
    if (max(pRld(lyr,s)) > 0)
        ylim([0 max(pRld(lyr,s))])
    else
        ylim([-0.1 0.1])
    end
    title('Demersal Fishes')
    
    subplot(3,1,2)
    plot(1:365,1000*pKlp(lyr,s),'--k'); hold on;
    plot(1:365,pRlp(lyr,s),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    xlim([0 365])
    if (max(pRlp(lyr,s)) > 0)
        ylim([0 max(pRlp(lyr,s))])
    else
        ylim([-0.1 0.1])
    end
    title('Pelagic Piscivores')
    ylabel('Reproductive output (g km^-^2)')
    xlabel('Time (y)')
    print('-dpng',[fpath loc '_oneloc_Rep_2000.png'])
end

%% Reprod compare
for s=1:length(spots)
    loc = spots{s};
    figure(4)
    clf
    subplot(3,1,2)
    plot(1:365,pRlp(lyr,s),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    plot(1:365,Rlp(lyr,s),'--','Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    ylabel('Reproductive output (g km^-^2)')
    xlim([0 365])
    title('Pelagic Piscivores')
    legend('phenology','constant')
    legend('location','eastoutside')
    
    subplot(3,1,3)
    plot(1:365,pRld(lyr,s),'Linewidth',2,'color',cmap_ppt(2,:)); hold on;
    plot(1:365,Rld(lyr,s),'--','Linewidth',2,'color',cmap_ppt(2,:)); hold on;
    xlim([0 365])
    legend('phenology','constant')
    legend('location','eastoutside')
    title('Demersal Fishes')
    
    subplot(3,1,1)
    plot(1:365,pRmf(lyr,s),'Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    plot(1:365,Rmf(lyr,s),'--','Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    xlim([0 365])
    title({[loc ' Spinup 1980']; ['Forage Fishes']})
    xlabel('Time (y)')
    legend('phenology','constant')
    legend('location','eastoutside')
    print('-dpng',[fpath loc '_oneloc_Rep_2000_compare.png'])
end

%% Compare locations
cmap_ppt(6,:) = [72/255 61/255 139/255]; %dark slate blue
%cmap_ppt(6,:) = [0/255 0/255 128/255]; %navy
% Degree days
for s=1:length(spots)
    loc = spots{s};
    figure(5)
    plot(1:365,pDDlp(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    xlim([0 365])
    title('Spinup 1980')
    xlabel('Time (y)')
    ylabel('Cumulative degree days')
    legend(spots)
    print('-dpng',[fpath 'All_oneloc_DD_2000.png'])
end

% Spawn flag
for s=1:length(spots)
    loc = spots{s};
    figure(6)
    plot(1:365,pKlp(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    xlim([0 365])
    ylim([0 0.05])
    title('Spinup 1980')
    xlabel('Time (y)')
    ylabel('Spawning flag')
    legend(spots)
    print('-dpng',[fpath 'All_oneloc_K_2000.png'])
end

%% Reprod
for s=1:length(spots)
    loc = spots{s};
    figure(6)
    subplot(3,1,2)
    plot(1:365,pRlp(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    ylabel('Reproductive output (g km^-^2)')
    xlim([0 365])
    legend(spots)
    legend('location','eastoutside')
    title('Pelagic Piscivores')
    
    subplot(3,1,1)
    plot(1:365,pRmf(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    xlim([0 365])
    legend(spots)
    legend('location','eastoutside')
    title({['Spinup 1980']; ['Forage Fishes']})
    xlabel('Time (y)')
    
    subplot(3,1,3)
    plot(1:365,pRld(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    xlim([0 365])
    legend(spots)
    legend('location','eastoutside')
    title('Demersal Piscivores')
    xlabel('Time (y)')
    print('-dpng',[fpath 'All_oneloc_Rep_2000.png'])
end

