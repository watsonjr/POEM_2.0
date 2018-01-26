% Plot effective TEs at LME scale
% Climatology
% 150 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));


%%
AREA_OCN = max(area,1);

cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([dpath 'TEeff_Climatol_All_fish03_' cfile '.mat']);

%% Calc LMEs
tlme = lme_mask_onedeg;

lme_te = NaN*ones(66,2);
for L=1:66
    lid = find(tlme==L);
    %TEeff
    lme_te(L,1) = nanmean(TEeffM(lid));
    lme_te(L,2) = nanmean(TEeffL(lid));
    
end

lme_m = NaN*ones(ni,nj);
lme_l = lme_m;
for L=1:66
    lid = find(tlme==L);

    lme_m(lid) = lme_te(L,1);
    lme_l(lid) = lme_te(L,2);
end

%%
save([dpath 'TEeff_Climatol_All_fish03_' cfile '.mat'],'lme_te',...
    'lme_m','lme_l','-append');

%% Figures
% M
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_m)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.005 0.03]);
hcb = colorbar('h');
ylim(hcb,[0.005 0.03])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEeff M')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffM.png'])

% L
figure(2)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lme_l)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0.0005 0.003]);
hcb = colorbar('h');
ylim(hcb,[0.0005 0.003])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('Climatology TEeff L')
stamp(cfile)
print('-dpng',[ppath 'Clim_fished_',harv,'_LME_TEeffL.png'])

