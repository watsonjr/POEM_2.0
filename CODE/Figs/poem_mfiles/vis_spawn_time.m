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
    subplot(3,1,3)
    plot(1:365,pRlp(lyr,s),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    xlim([0 365])
    title({[loc ' Historic']; ['Pelagic Piscivores']})
    ylabel('Reproductive output (g km^-^2)')
    
    subplot(3,1,2)
    plot(1:365,pRmd(lyr,s),'Linewidth',2,'color',cmap_ppt(5,:)); hold on;
    xlim([0 365])
    title('Demersal Fishes')
    ylabel('Reproductive output (g km^-^2)')
    
    subplot(3,1,1)
    plot(1:365,pRmf(lyr,s),'Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    xlim([0 365])
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('Reproductive output (g km^-^2)')
    print('-dpng',[fpath loc '_oneloc_Rep_2000.png'])
end

%% Reprod compare
for s=1:length(spots)
    close all
    loc = spots{s};
    figure(3)
    subplot(3,1,3)
    plot(1:365,log(pRlp(lyr,s)),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    plot(1:365,log(Rlp(lyr,s)),'Linewidth',2,'color','k'); hold on;
    xlim([0 365])
    title({[loc ' Historic']; ['Pelagic Piscivores']})
    ylabel('log Reproductive output (g km^-^2)')
    legend('phenology','constant')
    
    subplot(3,1,2)
    plot(1:365,pRmd(lyr,s),'Linewidth',2,'color',cmap_ppt(5,:)); hold on;
    xlim([0 365])
    title('Demersal Fishes')
    ylabel('Reproductive output (g km^-^2)')
    
    subplot(3,1,1)
    plot(1:365,log(pRmf(lyr,s)),'Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    plot(1:365,log(Rmf(lyr,s)),'Linewidth',2,'color',[0.5 0.5 0.5]); hold on;
    xlim([0 365])
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('log Reproductive output (g km^-^2)')
    legend('phenology','constant')
    print('-dpng',[fpath loc '_oneloc_Rep_2000_compare.png'])
end

%% Compare locations
cmap_ppt(6,:) = [72/255 61/255 139/255]; %dark slate blue
%cmap_ppt(6,:) = [0/255 0/255 128/255]; %navy
% Degree days
for s=1:length(spots)
    loc = spots{s};
    figure(4)
    plot(1:365,pDDlp(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    xlim([0 365])
    title(['Historic ' loc])
    xlabel('Time (y)')
    ylabel('Cumulative degree days')
    legend(spots)
    print('-dpng',[fpath 'All_oneloc_DD_2000.png'])
end

% Spawn flag
for s=1:length(spots)
    loc = spots{s};
    figure(5)
    plot(1:365,1-pKlp(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    xlim([0 365])
    ylim([-0.1 1.1])
    title(['Historic ' loc])
    xlabel('Time (y)')
    ylabel('Spawning flag')
    legend(spots)
    print('-dpng',[fpath 'All_oneloc_K_2000.png'])
end

%% Reprod
for s=1:length(spots)
    loc = spots{s};
    figure(6)
    subplot(3,1,1)
    plot(1:365,pRlp(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    xlim([0 365])
    %legend(spots)
    title({[loc ' Historic']; ['Pelagic Piscivores']})
    ylabel('Reproductive output (g km^-^2)')
    
    subplot(3,1,2)
    plot(1:365,pRmf(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    xlim([0 365])
    legend(spots)
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('Reproductive output (g km^-^2)')
    
    subplot(3,1,3)
    plot(1:365,pRmd(lyr,s),'Linewidth',2,'color',cmap_ppt(s,:)); hold on;
    xlim([0 365])
    legend(spots)
    title('Demersal Piscivores')
    xlabel('Time (y)')
    ylabel('Reproductive output (g km^-^2)')
    print('-dpng',[fpath 'All_oneloc_Rep_2000.png'])
end

%% Scaled Reprod compare
spRlp = pRlp(lyr,:) ./ repmat(sum(pRlp(lyr,:)),365,1);
sRlp = Rlp(lyr,:) ./ repmat(sum(Rlp(lyr,:)),365,1);
spRmf = pRmf(lyr,:) ./ repmat(sum(pRmf(lyr,:)),365,1);
sRmf = Rmf(lyr,:) ./ repmat(sum(Rmf(lyr,:)),365,1);
spRmd = pRmd(lyr,:) ./ repmat(sum(pRmd(lyr,:)),365,1);
sRmd = Rmd(lyr,:) ./ repmat(sum(Rmd(lyr,:)),365,1);

for s=1:length(spots)
    loc = spots{s};
    figure(3)
    clf
    subplot(3,1,1)
    plot(1:365,spRlp(:,s),'Linewidth',2,'color',cmap_ppt(1,:)); hold on;
    plot(1:365,sRlp(:,s),'Linewidth',2,'color','k'); hold on;
    xlim([0 365])
    title({[loc ' Historic']; ['Pelagic Piscivores']})
    ylabel('Fraction Annual Reproductive Output')
    legend('phenology','constant')
         
    subplot(3,1,2)
    plot(1:365,spRmf(:,s),'Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    plot(1:365,sRmf(:,s),'Linewidth',2,'color',[0.5 0.5 0.5]); hold on;
    xlim([0 365])
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('Fraction Annual Reproductive Output')
    legend('phenology','constant')
    
    subplot(3,1,3)
    plot(1:365,spRmd(:,s),'Linewidth',2,'color',cmap_ppt(3,:)); hold on;
    plot(1:365,sRmd(:,s),'Linewidth',2,'color',[0.5 0.5 0.5]); hold on;
    xlim([0 365])
    title('Demersal Piscivores')
    xlabel('Time (y)')
    ylabel('Fraction Annual Reproductive Output')
    legend('phenology','constant')
    print('-dpng',[fpath loc '_oneloc_Rep_2000_compare_scaled.png'])

    f20=figure(20);
    subplot(2,3,s)
    plot(1:365,sRlp(:,s),'Linewidth',2); hold on;
    plot(1:365,spRlp(:,s),'Linewidth',2); hold on;
    xlim([0 365])
    title(loc)
    if (s==1)
    ylabel('Fraction Annual Reproductive Output')
    end
    if (s==4)
    ylabel('Forage Fish')
    legend('C','P')
    legend('location','northwest')
    end
    
    f21=figure(21);
    subplot(2,3,s)
    plot(1:365,sRmf(:,s),'Linewidth',2); hold on;
    plot(1:365,spRmf(:,s),'Linewidth',2); hold on;
    xlim([0 365])
    title(loc)
    if (s==1)
    ylabel('Fraction Annual Reproductive Output')
    end
    if (s==4)
    ylabel('Pelagic Piscivore')
    legend('C','P')
    legend('location','northwest')
    end
    
    f22=figure(22);
    subplot(2,3,s)
    plot(1:365,sRmd(:,s),'Linewidth',2); hold on;
    plot(1:365,spRmd(:,s),'Linewidth',2); hold on;
    xlim([0 365])
    title(loc)
    if (s==1)
    ylabel('Fraction Annual Reproductive Output')
    end
    if (s==4)
    ylabel('Demersal Piscivore')
    legend('C','P')
    legend('location','northwest')
    end

end
print(f20,'-dpng',[fpath 'Forage_oneloc_Rep_2000_compare_scaled.png'])
print(f21,'-dpng',[fpath 'Pisc_oneloc_Rep_2000_compare_scaled.png'])
print(f22,'-dpng',[fpath 'Dem_oneloc_Rep_2000_compare_scaled.png'])
