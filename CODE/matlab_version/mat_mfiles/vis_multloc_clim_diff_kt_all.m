% Visualize output of POEM
% ESM2.6 Climatology of 5 yrs
% 150 years
% Saved as nc files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%
harv = 'All_fish03';
tharv = 'Harvest all fish 0.3 yr^-^1';

close all

% plot info
[ni,nj]=size(lon);
geolat_t=lat;
geolon_t=lon;
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));


% colors
load('MyColormaps.mat')
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
    0 0 0];...      %black
    
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm9);

%% Plots in space

cfileA = 'Dc_enc70-b200_cm20_m-b175-k05_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathA=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileA '/'];
load([fpathA 'Clim_means_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Afish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Aall=NaN*ones(ni,nj);
Aall(ID)=Afish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileB = 'Dc_enc70-b200_cm20_m-b175-k07_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathB=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileB '/'];
load([fpathB 'Clim_means_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Bfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Ball=NaN*ones(ni,nj);
Ball(ID)=Bfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileC = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathC=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileC '/'];
load([fpathC 'Clim_means_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Cfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Call=NaN*ones(ni,nj);
Call(ID)=Cfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

cfileD = 'Dc_enc70-b200_cm20_m-b175-k11_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
fpathD=['/Volumes/GFDL/CSV/Matlab_new_size/' cfileD '/'];
load([fpathD 'Clim_means_All_fish03.mat'],'sf_mean','mf_mean','sp_mean',...
    'mp_mean','lp_mean','sd_mean','md_mean','ld_mean');
Dfish=sf_mean+sp_mean+sd_mean+mf_mean+mp_mean+md_mean+lp_mean+ld_mean;
Dall=NaN*ones(ni,nj);
Dall(ID)=Dfish;
clear sf_mean mf_mean sp_mean mp_mean lp_mean sd_mean md_mean ld_mean

%% 
figure(1)
%A
subplot('Position',[0 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Aall))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'A')

%B
subplot('Position',[0.5 0.51 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Ball))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'B')

%C
subplot('Position',[0 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Call))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
text(-2.75,1.75,'C')

%D
subplot('Position',[0.5 0.1 0.5 0.4])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Dall))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
colorbar('Position',[0.25 0.075 0.5 0.03],'orientation','horizontal')
set(gcf,'renderer','painters')
text(-2.75,1.75,'D')
%stamp([harv '_' cfile])
print('-dpng',[pp 'Climatol_' harv '_All_comp_kt.png'])

