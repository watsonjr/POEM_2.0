% Compare continuous vs spawning season
% Pristine historical at one location
% 145 years

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

load('Oneloc_pris_phenol_all.mat')
%load('Oneloc_pris_all.mat')
load('cmap_ppt_angles.mat')

%% Plots over time
x=1:length(psp);
yfrac=x/365;
y=1861:(1/365):(1861+yfrac(end));
%y=1861:(1/365):2005;
y=y(1:(end-1));
lstd=length(psp);
id1 = 0:365:(lstd-1);
id2 = 365:365:(lstd);
ID  = [id1 id2];

%year = 2000
lyr=50736:(50736+364);

%% Degree days
for s=1:length(spots)
    close all
    loc = spots{s};
    figure(s)
    plot(1:365,pDDlp(lyr,s),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    plot(1:365,pDDmd(lyr,s),'Linewidth',2,'color',cmap_ppt(5,:)); hold on;
    plot(1:365,pDDmf(lyr,s),'--','Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    %colormap(cmap_ppt(1:3,:))
    xlim([0 365])
    title(['Historic ' loc])
    xlabel('Time (y)')
    ylabel('Cumulative degree days')
    legend('Pisc','Demersal','Forage')
    print('-dpng',[fpath loc '_oneloc_DD_2000.png'])
end

%% Spawn flag
for s=1:length(spots)
    close all
    loc = spots{s};
    figure(2)
    plot(1:365,1-pKlp(lyr,s),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    plot(1:365,1-pKmd(lyr,s),'Linewidth',2,'color',cmap_ppt(5,:)); hold on;
    plot(1:365,1-pKmf(lyr,s),'--','Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    xlim([0 365])
    ylim([-0.1 1.1])
    title(['Historic ' loc])
    xlabel('Time (y)')
    ylabel('Spawning flag')
    legend('Pisc','Demersal','Forage')
    print('-dpng',[fpath loc '_oneloc_K_2000.png'])
end

%% Reprod
for s=1:length(spots)
    close all
    loc = spots{s};
    figure(3)
    subplot(2,1,1)
    plot(1:365,pRlp(lyr,s),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    xlim([0 365])
    title({[loc ' Historic']; ['Pelagic Piscivores']})
    ylabel('Reproductive output (g km^-^2)')
    
%     subplot(3,1,2)
%     plot(1:365,pRmd(lyr,s),'Linewidth',2,'color',cmap_ppt(5,:)); hold on;
%     xlim([0 365])
%     title('Demersal Fishes')
%     ylabel('Reproductive output (g km^-^2)')
%     
    subplot(2,1,2)
    plot(1:365,pRmf(lyr,s),'Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    xlim([0 365])
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('Reproductive output (g km^-^2)')
    print('-dpng',[fpath loc '_oneloc_Rep_2000.png'])
end