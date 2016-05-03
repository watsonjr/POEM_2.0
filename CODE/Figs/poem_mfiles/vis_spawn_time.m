% Compare continuous vs spawning season
% Pristine historical at one location
% 145 years

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

load('Oneloc_pris_phenol_all.mat')
load('Oneloc_pris_all.mat')
load('cmap_ppt_angles.mat')

%% Plots over time
x=1:length(sp);
yfrac=x/365;
y=1861:(1/365):(1861+yfrac(end));
%y=1861:(1/365):2005;
y=y(1:(end-1));
lstd=length(sp);
id1 = 0:365:(lstd-1);
id2 = 365:365:(lstd);
ID  = [id1 id2];

%last year = 2005
lyr=x((end-364):end);

%% Degree days
for s=1:length(spots)
    %close all
    loc = spots{s};
    figure(s)
    %plot(1:365,pDDlp(lyr,s),'Linewidth',2); hold on;
    plot(x,pDDlp(:,s),'Linewidth',2); hold on;
    %xlim([0 365])
    title(['Historic Pelagic Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('Cumulative degree days')
    %print('-dpng',[fpath sname 'oneloc_pisc_time.png'])
    
    figure
    plot(x,pDDmd(:,s),'Linewidth',2); hold on;
    %xlim([0 365])
    title(['Historic Demersal ' loc])
    xlabel('Time (y)')
    ylabel('Cumulative degree days')
    
    figure
    plot(x,pDDmf(:,s),'Linewidth',2); hold on;
    %xlim([0 365])
    title(['Historic Forage Fish ' loc])
    xlabel('Time (y)')
    ylabel('Cumulative degree days')
end

%% Spawn flag
for s=1%:length(spots)
    close all
    loc = spots{s};
    figure(1)
    plot(1:365,1-pKlp(lyr,s),'Linewidth',2); hold on;
    xlim([0 365])
    title(['Historic Pelagic Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('Spawning flag')
    %print('-dpng',[fpath sname 'oneloc_pisc_time.png'])
    
    figure
    plot(1:365,1-pKmd(lyr,s),'Linewidth',2); hold on;
    xlim([0 365])
    title(['Historic Demersal ' loc])
    xlabel('Time (y)')
    ylabel('Spawning flag')
    
    figure
    plot(1:365,1-pKmf(lyr,s),'Linewidth',2); hold on;
    xlim([0 365])
    title(['Historic Forage Fish ' loc])
    xlabel('Time (y)')
    ylabel('Spawning flag')
end

%% Reprod
for s=1%:length(spots)
    close all
    loc = spots{s};
    figure(1)
    plot(1:365,pRlp(lyr,s),'Linewidth',2); hold on;
    xlim([0 365])
    title(['Historic Pelagic Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('Reproductive output (g km^-^2)')
    %print('-dpng',[fpath sname 'oneloc_pisc_time.png'])
    
    figure
    plot(1:365,pRmd(lyr,s),'Linewidth',2); hold on;
    xlim([0 365])
    title(['Historic Demersal ' loc])
    xlabel('Time (y)')
    ylabel('Reproductive output (g km^-^2)')
    
    figure
    plot(1:365,pRmf(lyr,s),'Linewidth',2); hold on;
    xlim([0 365])
    title(['Historic Forage Fish ' loc])
    xlabel('Time (y)')
    ylabel('Reproductive output (g km^-^2)')
end