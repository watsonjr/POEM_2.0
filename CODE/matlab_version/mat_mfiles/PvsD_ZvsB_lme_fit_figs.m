% Plot mean biomasses of zoop to compare to POEM results
% COMPARE BY LME, COLOR-CODE BY TEMP, Plot GAM fit

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

gpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
cpath = [Pdrpbx 'Princeton/POEM_other/cobalt_data/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([gpath 'esm26_area_1deg.mat']);
load([gpath 'LME_clim_temp_zoop_det.mat']);
load([gpath 'LME_depth_area.mat'],'lme_depth','lme_shal_frac');

%% POEM results
frate = 0.3;
tfish = num2str(100+int64(10*frate));
cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = ['All_fish',tfish(2:end)];
tharv = ['Harvest all fish ' num2str(frate)];
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end
load([fpath 'LME_ZBratios_clim_fished_',harv,'_' cfile '.mat']);
load([fpath 'ZlB_gam_fits.mat']);

%% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

%% Scatter plot
% Plot GAM fit
figure(1)
subplot(3,2,1)
for i=1:66
    plot(log10(RatZlB(i)),FracPD(i),'.','MarkerSize',15,'color',tmap(tid(i,2),:)); hold on;
end
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+D)')
axis([-6.6 -5 0 1])

subplot(3,2,3)
for i=1:66
    plot(log10(RatZlB(i)),FracPF(i),'.','MarkerSize',15,'color',tmap(tid(i,2),:)); hold on;
end
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+F)')
axis([-6.6 -5 0 1])

subplot(3,2,5)
for i=1:66
    plot(log10(RatZlB(i)),FracLM(i),'.','MarkerSize',15,'color',tmap(tid(i,2),:)); hold on;
end
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Bent')
ylabel('L / (L+M)')
axis([-6.6 -5 0 1])
print('-dpng',[ppath 'lme_scatter_ratio_bent_GAMfit_colorT.png'])


figure(2)
subplot(3,2,1)
plot(log10(RatZlB),FracPD,'.k','MarkerSize',15); hold on;
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+D)')
axis([-6.6 -5 0 1])

subplot(3,2,3)
plot(log10(RatZlB),FracPF,'.k','MarkerSize',15); hold on;
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+F)')
axis([-6.6 -5 0 1])

subplot(3,2,5)
plot(log10(RatZlB),FracLM,'.k','MarkerSize',15); hold on;
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Bent')
ylabel('L / (L+M)')
axis([-6.6 -5 0 1])
print('-dpng',[ppath 'lme_scatter_ratio_bent_GAMfit_BW.png'])

%% Horiztonal
figure(1)
subplot(1,3,1)
for i=1:66
    plot(log10(RatZlB(i)),FracPD(i),'.','MarkerSize',15,'color',tmap(tid(i,2),:)); hold on;
end
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+D)')
axis([-6.6 -5 0 1])

subplot(1,3,2)
for i=1:66
    plot(log10(RatZlB(i)),FracPF(i),'.','MarkerSize',15,'color',tmap(tid(i,2),:)); hold on;
end
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+F)')
axis([-6.6 -5 0 1])

subplot(1,3,3)
for i=1:66
    plot(log10(RatZlB(i)),FracLM(i),'.','MarkerSize',15,'color',tmap(tid(i,2),:)); hold on;
end
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Bent')
ylabel('L / (L+M)')
axis([-6.6 -5 0 1])
print('-dpng',[ppath 'lme_scatter_ratio_bent_GAMfit_colorTh.png'])


figure(2)
subplot(1,3,1)
plot(log10(RatZlB),FracPD,'.k','MarkerSize',15); hold on;
plot(PDgam(:,3),PDgam(:,1),'k'); hold on;
plot(PDgam(:,3),PDgam(:,1)+2*PDgam(:,2),'--k'); hold on;
plot(PDgam(:,3),PDgam(:,1)-2*PDgam(:,2),'--k'); hold on;
%title('Color = LME mean temp')
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+D)')
axis([-6.6 -5 0 1])

subplot(1,3,2)
plot(log10(RatZlB),FracPF,'.k','MarkerSize',15); hold on;
plot(PFgam(:,3),PFgam(:,1),'k'); hold on;
plot(PFgam(:,3),PFgam(:,1)+2*PFgam(:,2),'--k'); hold on;
plot(PFgam(:,3),PFgam(:,1)-2*PFgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Bent')
ylabel('P / (P+F)')
axis([-6.6 -5 0 1])

subplot(1,3,3)
plot(log10(RatZlB),FracLM,'.k','MarkerSize',15); hold on;
plot(LMgam(:,3),LMgam(:,1),'k'); hold on;
plot(LMgam(:,3),LMgam(:,1)+2*LMgam(:,2),'--k'); hold on;
plot(LMgam(:,3),LMgam(:,1)-2*LMgam(:,2),'--k'); hold on;
xlabel('log10 ZoopLoss:Bent')
ylabel('L / (L+M)')
axis([-6.6 -5 0 1])
print('-dpng',[ppath 'lme_scatter_ratio_bent_GAMfit_BWh.png'])




