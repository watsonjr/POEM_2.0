% Changes to global biomass 
% Changes to biomass in Pacific biomes

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

cfile = 'Dc_enc70_cmax-metab20_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00400';
harv = '03';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

load([gpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([gpath 'lme_mask_esm2m.mat']);
load([gpath 'NPid_esm2m.mat']);
grid = csvread([gpath 'grid_csv.csv']);

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')
cmap1(1,:)=[1 1 1];
cmap1(2,:)=cmap_ppt(1,:);
cmap1(3,:)=cmap_ppt(3,:);
cmap1(4,:)=cmap_ppt(5,:);

cmap2(1,:)=cmap_ppt(1,:);
cmap2(2,:)=cmap_ppt(3,:);
cmap2(3,:)=cmap_ppt(5,:);

cmap3(1,:)=[0 0 0];
cmap3(2,:)=cmap_ppt(1,:);
cmap3(3,:)=cmap_ppt(3,:);
cmap3(4,:)=cmap_ppt(5,:);

cmap4(1,:)=cmap_ppt(3,:);
cmap4(2,:)=cmap_ppt(1,:);
cmap4(3,:)=cmap_ppt(5,:);

% figure
% y=ones(10,1);
% plot(1:10,y,'color',cmap_ppt(1,:)); hold on;    %blue
% plot(1:10,2*y,'color',cmap_ppt(2,:)); hold on;  %grey
% plot(1:10,3*y,'color',cmap_ppt(3,:)); hold on;  %orange
% plot(1:10,4*y,'color',cmap_ppt(4,:)); hold on;  %black
% plot(1:10,5*y,'color',cmap_ppt(5,:)); hold on;  %green
% plot(1:10,6*y,'color',cmap_ppt(6,:)); hold on;  %dark grey
% plot(1:10,7*y,'color',cmap_ppt(7,:)); hold on;  %light grey

% Biomes 1 = LC; 2 = ECCS; 3 = ECSS;

%% Hist fished
load([dpath 'LME_hist_fished',harv,'_' cfile '.mat'],'lme_bio50','lme_mbio50');
load([dpath 'Biomes_hist_fished',harv,'_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_hist_fished',harv,'_' cfile '.mat'],'sf_mean5000','sp_mean5000',...
    'sd_mean5000','mf_mean5000','mp_mean5000','md_mean5000','b_mean5000',...
    'lp_mean5000','ld_mean5000');

histf_g = lme_bio00;
histf_gm2 = lme_mbio00;

histf_biome = biome_mbio50;

%Mean global biomass in g/m2
all_histf(:,1) = sf_mean5000+sp_mean5000+sd_mean5000+mf_mean5000+...
    mp_mean5000+md_mean5000+lp_mean5000+ld_mean5000;
all_histf(:,2) = sf_mean5000+mf_mean5000; %all f
all_histf(:,3) = sp_mean5000+mp_mean5000+lp_mean5000; %all p
all_histf(:,4) = sd_mean5000+md_mean5000+ld_mean5000; %all d
all_histf(:,5) = sf_mean5000+sp_mean5000+sd_mean5000; %all s
all_histf(:,6) = mf_mean5000+mp_mean5000+md_mean5000; %all m
all_histf(:,7) = lp_mean5000+ld_mean5000; %all l

%Mean total global biomass in g
tot_all_histf = all_histf .* repmat(grid(:,5),1,7);

clear lme_bio00 lme_mbio00 biome_mbio50 sf_mean5000 sp_mean5000 sd_mean5000
clear mf_mean5000 mp_mean5000 md_mean5000 b_mean5000 lp_mean5000 ld_mean5000

%% Fore fished
load([dpath 'LME_fore_fished_' cfile '.mat'],'lme_bio00','lme_mbio00');
load([dpath 'Biomes_fore_fished_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_fore_fished_' cfile '.mat'],'sf_mean5000','sp_mean5000',...
    'sd_mean5000','mf_mean5000','mp_mean5000','md_mean5000','b_mean5000',...
    'lp_mean5000','ld_mean5000');

foref_g = lme_bio00;
foref_gm2 = lme_mbio00;

foref_biome = biome_mbio50;

%Mean global biomass in g/m2
all_foref(:,1) = sf_mean5000+sp_mean5000+sd_mean5000+mf_mean5000+...
    mp_mean5000+md_mean5000+lp_mean5000+ld_mean5000;
all_foref(:,2) = sf_mean5000+mf_mean5000; %all f
all_foref(:,3) = sp_mean5000+mp_mean5000+lp_mean5000; %all p
all_foref(:,4) = sd_mean5000+md_mean5000+ld_mean5000; %all d
all_foref(:,5) = sf_mean5000+sp_mean5000+sd_mean5000; %all s
all_foref(:,6) = mf_mean5000+mp_mean5000+md_mean5000; %all m
all_foref(:,7) = lp_mean5000+ld_mean5000; %all l

%Mean total global biomass in g
tot_all_foref = all_foref .* repmat(grid(:,5),1,7);

clear lme_bio00 lme_mbio00 biome_mbio50 sf_mean5000 sp_mean5000 sd_mean5000
clear mf_mean5000 mp_mean5000 md_mean5000 b_mean5000 lp_mean5000 ld_mean5000

%% Global

Msum_histf = nansum(all_histf);
Msum_foref = nansum(all_foref);

Tsum_histf = nansum(tot_all_histf);
Tsum_foref = nansum(tot_all_foref);

%Future fishing vs. historic fishing
Mdiff_foref_histf = (Msum_foref - Msum_histf) ./ Msum_histf;
Tdiff_foref_histf = (Tsum_foref - Tsum_histf) ./ Tsum_histf;

%% Time series
load([dpath 'All_biom_mcatch_hist_pristine.mat']);
all_bio_ts = nansum(all_bio);
ts_bio_histp = all_bio_ts(2:16);
F_bio_ts = nansum(F_bio);
Fts_bio_histp = F_bio_ts(2:16);
P_bio_ts = nansum(P_bio);
Pts_bio_histp = P_bio_ts(2:16);
D_bio_ts = nansum(D_bio);
Dts_bio_histp = D_bio_ts(2:16);
clear all_bio all_bio_ts F_bio F_bio_ts P_bio P_bio_ts D_bio D_bio_ts

load([dpath 'All_biom_mcatch_hist_fished.mat']);
all_bio_ts = nansum(all_bio);
ts_bio_histf = all_bio_ts(2:16);
all_mcatch_ts = nansum(all_mcatch);
ts_mcatch_histf = all_mcatch_ts(2:16);
F_bio_ts = nansum(F_bio);
Fts_bio_histf = F_bio_ts(2:16);
P_bio_ts = nansum(P_bio);
Pts_bio_histf = P_bio_ts(2:16);
D_bio_ts = nansum(D_bio);
Dts_bio_histf = D_bio_ts(2:16);
F_mcatch_ts = nansum(F_mcatch);
Fts_mcatch_histf = F_mcatch_ts(2:16);
P_mcatch_ts = nansum(P_mcatch);
Pts_mcatch_histf = P_mcatch_ts(2:16);
D_mcatch_ts = nansum(D_mcatch);
Dts_mcatch_histf = D_mcatch_ts(2:16);
clear all_bio all_bio_ts F_bio F_bio_ts P_bio P_bio_ts D_bio D_bio_ts
clear all_mcatch all_mcatch_ts F_mcatch F_mcatch_ts P_mcatch P_mcatch_ts D_mcatch D_mcatch_ts

%Only fcrit30 can't use
% load([dpath 'All_biom_fore_pristine.mat']);
% all_bio_ts = nansum(all_bio);
% ts_bio_forep = all_bio_ts(2:11);
% F_bio_ts = nansum(F_bio);
% Fts_bio_forep = F_bio_ts(2:11);
% P_bio_ts = nansum(P_bio);
% Pts_bio_forep = P_bio_ts(2:11);
% D_bio_ts = nansum(D_bio);
% Dts_bio_forep = D_bio_ts(2:11);
% clear all_bio all_bio_ts F_bio F_bio_ts P_bio P_bio_ts D_bio D_bio_ts

load([dpath 'All_biom_mcatch_fore_fished.mat']);
all_bio_ts = nansum(all_bio);
ts_bio_foref = all_bio_ts(2:11);
all_mcatch_ts = nansum(all_mcatch);
ts_mcatch_foref = all_mcatch_ts(2:11);
F_bio_ts = nansum(F_bio);
Fts_bio_foref = F_bio_ts(2:11);
P_bio_ts = nansum(P_bio);
Pts_bio_foref = P_bio_ts(2:11);
D_bio_ts = nansum(D_bio);
Dts_bio_foref = D_bio_ts(2:11);
F_mcatch_ts = nansum(F_mcatch);
Fts_mcatch_foref = F_mcatch_ts(2:11);
P_mcatch_ts = nansum(P_mcatch);
Pts_mcatch_foref = P_mcatch_ts(2:11);
D_mcatch_ts = nansum(D_mcatch);
Dts_mcatch_foref = D_mcatch_ts(2:11);
clear all_bio all_bio_ts F_bio F_bio_ts P_bio P_bio_ts D_bio D_bio_ts
clear all_mcatch all_mcatch_ts F_mcatch F_mcatch_ts P_mcatch P_mcatch_ts D_mcatch D_mcatch_ts

%%
yrh=1865:10:2005;
yrf=2005:10:2100;

load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

figure(1)
plot(yrh,Fts_bio_histf,'color',cmap_ppt(3,:),'LineWidth',2); hold on;
plot(yrh,Pts_bio_histf,'color',cmap_ppt(1,:),'LineWidth',2); hold on;
plot(yrh,Dts_bio_histf,'color',cmap_ppt(5,:),'LineWidth',2); hold on;
plot(yrf,Fts_bio_foref,'color',cmap_ppt(3,:),'LineWidth',2); hold on;
plot(yrf,Pts_bio_foref,'color',cmap_ppt(1,:),'LineWidth',2); hold on;
plot(yrf,Dts_bio_foref,'color',cmap_ppt(5,:),'LineWidth',2); hold on;
xlim([1860 2100])
xlabel('Year')
ylabel('Mean biomass (g/m^2)')
legend('F','P','D')
%legend('location','southwest')
title('Mean biomass')
print('-dpng',[ppath 'ts_mbio_type.png'])
%%
figure(2)
plot(yrh,Fts_mcatch_histf*365,'color',cmap_ppt(3,:),'LineWidth',2); hold on;
plot(yrh,Pts_mcatch_histf*365,'color',cmap_ppt(1,:),'LineWidth',2); hold on;
plot(yrh,Dts_mcatch_histf*365,'color',cmap_ppt(5,:),'LineWidth',2); hold on;
plot(yrf,Fts_mcatch_foref*365,'color',cmap_ppt(3,:),'LineWidth',2); hold on;
plot(yrf,Pts_mcatch_foref*365,'color',cmap_ppt(1,:),'LineWidth',2); hold on;
plot(yrf,Dts_mcatch_foref*365,'color',cmap_ppt(5,:),'LineWidth',2); hold on;
xlim([1860 2100])
ylim([1500 6500])
xlabel('Year')
ylabel('Mean annual catch (g/m^2)')
legend('F','P','D')
%legend('location','southwest')
title('Mean catch')
print('-dpng',[ppath 'ts_mcatch_type.png'])



%% N Pac
[npid,ix,grid_np] = intersect(NPid,grid(:,1));

Msum_histf_np = nansum(all_histf(grid_np,:));
Msum_foref_np = nansum(all_foref(grid_np,:));

Tsum_histf_np = nansum(tot_all_histf(grid_np,:));
Tsum_foref_np = nansum(tot_all_foref(grid_np,:));

%Future fishing vs. historic fishing
Mdiff_foref_histf_np = (Msum_foref_np - Msum_histf_np) ./ Msum_histf_np;
Tdiff_foref_histf_np = (Tsum_foref_np - Tsum_histf_np) ./ Tsum_histf_np;

%% Maps
%plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

%Global grid
All_histf_gg = griddata(grid(:,2),grid(:,3),all_histf(:,1),X,Y);
All_foref_gg = griddata(grid(:,2),grid(:,3),all_foref(:,1),X,Y);
F_histf_gg = griddata(grid(:,2),grid(:,3),all_histf(:,2),X,Y);
F_foref_gg = griddata(grid(:,2),grid(:,3),all_foref(:,2),X,Y);
P_histf_gg = griddata(grid(:,2),grid(:,3),all_histf(:,3),X,Y);
P_foref_gg = griddata(grid(:,2),grid(:,3),all_foref(:,3),X,Y);
D_histf_gg = griddata(grid(:,2),grid(:,3),all_histf(:,4),X,Y);
D_foref_gg = griddata(grid(:,2),grid(:,3),all_foref(:,4),X,Y);

%MOM grid
All_histf_mg = NaN*ones(size(geolon_t));
All_foref_mg = NaN*ones(size(geolon_t));
All_histf_mg(grid(:,1)) = all_histf(:,1);
All_foref_mg(grid(:,1)) = all_foref(:,1);

F_histf_mg = NaN*ones(size(geolon_t));
F_foref_mg = NaN*ones(size(geolon_t));
F_histf_mg(grid(:,1)) = all_histf(:,2);
F_foref_mg(grid(:,1)) = all_foref(:,2);

P_histf_mg = NaN*ones(size(geolon_t));
P_foref_mg = NaN*ones(size(geolon_t));
P_histf_mg(grid(:,1)) = all_histf(:,3);
P_foref_mg(grid(:,1)) = all_foref(:,3);

D_histf_mg = NaN*ones(size(geolon_t));
D_foref_mg = NaN*ones(size(geolon_t));
D_histf_mg(grid(:,1)) = all_histf(:,4);
D_foref_mg(grid(:,1)) = all_foref(:,4);

Diff_foref_histf_gg = (All_foref_gg - All_histf_gg) ./All_histf_gg;
Diff_foref_histf_mg = (All_foref_mg - All_histf_mg) ./All_histf_mg;
FDiff_foref_histf_gg = (F_foref_gg - F_histf_gg) ./F_histf_gg;
FDiff_foref_histf_mg = (F_foref_mg - F_histf_mg) ./F_histf_mg;
PDiff_foref_histf_gg = (P_foref_gg - P_histf_gg) ./P_histf_gg;
PDiff_foref_histf_mg = (P_foref_mg - P_histf_mg) ./P_histf_mg;
DDiff_foref_histf_gg = (D_foref_gg - D_histf_gg) ./D_histf_gg;
DDiff_foref_histf_mg = (D_foref_mg - D_histf_mg) ./D_histf_mg;

%% Global
%All
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,Diff_foref_histf_gg*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference from 1951-2000 mean biomass of all fishes')
colormap(cmap_color_rb)
colorbar('h')
caxis([-100 100])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_all_foref_histf.png'])

%F
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,FDiff_foref_histf_gg*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference from 1951-2000 mean biomass of forage fishes')
colormap(cmap_color_rb)
colorbar('h')
caxis([-100 100])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_F_foref_histf.png'])

%P
figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,PDiff_foref_histf_gg*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference from 1951-2000 mean biomass of large pelagic fishes')
colormap(cmap_color_rb)
colorbar('h')
caxis([-100 100])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_P_foref_histf.png'])

%D
figure(6)
m_proj('miller','lat',82);
m_pcolor(X,Y,DDiff_foref_histf_gg*100); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2051-2100 % difference from 1951-2000 mean biomass of demersal fishes')
colormap(cmap_color_rb)
colorbar('h')
caxis([-100 100])
stamp(cfile)
print('-dpng',[ppath 'Diff_global_D_foref_histf.png'])

%% N Pac
%All
figure(7)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,Diff_foref_histf_mg*100)
colormap(cmap_color_rb)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
hcb = colorbar('h');
ylim(hcb,[-100 100])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1951-2000 mean biomass of all fishes')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_all_foref_histf.png'])

%F
figure(8)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,FDiff_foref_histf_mg*100)
colormap(cmap_color_rb)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
hcb = colorbar('h');
ylim(hcb,[-100 100])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1951-2000 mean biomass of forage fishes')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_F_foref_histf.png'])

%P
figure(9)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,PDiff_foref_histf_mg*100)
colormap(cmap_color_rb)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
hcb = colorbar('h');
ylim(hcb,[-100 100])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1951-2000 mean biomass of large pelagic fishes')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_P_foref_histf.png'])

%D
figure(10)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,DDiff_foref_histf_mg*100)
colormap(cmap_color_rb)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-100 100]);
hcb = colorbar('h');
ylim(hcb,[-100 100])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2051-2100 % difference from 1951-2000 mean biomass of demersal fishes')
stamp(cfile)
print('-dpng',[ppath 'Diff_NPac_D_foref_histf.png'])




