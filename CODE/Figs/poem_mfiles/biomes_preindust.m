% Calc biome biomass of POEM
% MOM grid cells in 3 biomes
% LC: surf chl < 0.15 mg/m3
% ECCS: surf chl > 0.15 mg/m3; max MLD < 75 m
% ECSS: surf chl > 0.15 mg/m3; max MLD > 75 m

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Means_preindust_' cfile '.mat']);
%load([dpath 'Preindust_1800-1850_means.mat']); %fcrit30
load([cpath 'COBALT_biomes.mat'],'biome_pre');
load([gpath 'hindcast_gridspec.mat'],'geolat_t','geolon_t');
grid = csvread([gpath 'grid_csv.csv']);

%% Calc biome biomass avgs
gid = grid(:,1);

sf_biome_bio = NaN*ones(3,1);
sf_biome_mbio = sf_biome_bio;
sp_biome_mbio = sf_biome_bio;
sd_biome_mbio = sf_biome_bio;
mf_biome_mbio = sf_biome_bio;
mp_biome_mbio = sf_biome_bio;
md_biome_mbio = sf_biome_bio;
lp_biome_mbio = sf_biome_bio;
ld_biome_mbio = sf_biome_bio;
b_biome_mbio = sf_biome_bio;


sf = NaN*ones(size(geolon_t));
sp = sf;
sd = sf;
mf = sf;
mp = sf;
md = sf;
lp = sf;
ld = sf;
b = sf;

sf(gid) = sf_smean(:,1);
sp(gid) = sp_smean(:,1);
sd(gid) = sd_smean(:,1);
mf(gid) = mf_smean(:,1);
mp(gid) = mp_smean(:,1);
md(gid) = md_smean(:,1);
lp(gid) = lp_smean(:,1);
ld(gid) = ld_smean(:,1);
b(gid) = b_smean(:,1);

for L=1:3
    lid = find(biome_pre==L);
    sf_biome_mbio(L,1) = nanmean(sf(lid));
    sp_biome_mbio(L,1) = nanmean(sp(lid));
    sd_biome_mbio(L,1) = nanmean(sd(lid));
    mf_biome_mbio(L,1) = nanmean(mf(lid));
    mp_biome_mbio(L,1) = nanmean(mp(lid));
    md_biome_mbio(L,1) = nanmean(md(lid));
    lp_biome_mbio(L,1) = nanmean(lp(lid));
    ld_biome_mbio(L,1) = nanmean(ld(lid));
    b_biome_mbio(L,1) = nanmean(b(lid));
end

biome_mbio(:,1) = sf_biome_mbio;
biome_mbio(:,2) = sp_biome_mbio;
biome_mbio(:,3) = sd_biome_mbio;
biome_mbio(:,4) = mf_biome_mbio;
biome_mbio(:,5) = mp_biome_mbio;
biome_mbio(:,6) = md_biome_mbio;
biome_mbio(:,7) = lp_biome_mbio;
biome_mbio(:,8) = ld_biome_mbio;
biome_mbio(:,9) = b_biome_mbio;

save([dpath 'Biomes_preindust_' cfile '.mat'],'biome_mbio',...
    'sf_biome_mbio','sp_biome_mbio','sd_biome_mbio','mf_biome_mbio','mp_biome_mbio',...
    'md_biome_mbio','b_biome_mbio','lp_biome_mbio','ld_biome_mbio');


%% Figures
biome_sf = NaN*ones(360,200);
biome_sp = biome_sf;
biome_sd = biome_sf;
biome_mf = biome_sf;
biome_mp = biome_sf;
biome_md = biome_sf;
biome_lp = biome_sf;
biome_ld = biome_sf;
biome_b = biome_sf;

for L=1:3
    lid = find(biome_pre==L);
    
    biome_sf(lid) = sf_biome_mbio(L,1);
    biome_sp(lid) = sp_biome_mbio(L,1);
    biome_sd(lid) = sd_biome_mbio(L,1);
    biome_mf(lid) = mf_biome_mbio(L,1);
    biome_mp(lid) = mp_biome_mbio(L,1);
    biome_md(lid) = md_biome_mbio(L,1);
    biome_lp(lid) = lp_biome_mbio(L,1);
    biome_ld(lid) = ld_biome_mbio(L,1);
    biome_b(lid) = b_biome_mbio(L,1);
end

biome_All = biome_sf+biome_sp+biome_sd+biome_mf+biome_mp+biome_md+biome_lp+biome_ld;
biome_AllF = biome_sf+biome_mf;
biome_AllP = biome_sp+biome_mp+biome_lp;
biome_AllD = biome_sd+biome_md+biome_ld;

%% Check map
% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];
% ENTER -100 TO MAP ORIGIN LONG

%% ALL
figure(1)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_All)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
hcb = colorbar('h');
ylim(hcb,[-0.5 0.5])                    %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
axesmui
stamp(cfile)
print('-dpng',[ppath 'Preindust_biomes_All.png'])

%% all F
figure(2)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_AllF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Forage Fishes (g m^-^2)')
axesmui
stamp(cfile)
print('-dpng',[ppath 'Preindust_biomes_AllF.png'])

% all P
figure(3)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_AllP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Pelagic Fishes (g m^-^2)')
axesmui
stamp(cfile)
print('-dpng',[ppath 'Preindust_biomes_AllP.png'])

% All D
figure(4)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_AllD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Demersal Fishes (g m^-^2)')
axesmui
stamp(cfile)
print('-dpng',[ppath 'Preindust_biomes_AllD.png'])

%% Pac only
lonlim=[-255 -60];
latlim=[0 80];

% ALL
figure(41)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_All)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.5 0.5]);
hcb = colorbar('h');
ylim(hcb,[-0.5 0.5])                    %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Preindust_Pac_biomes_All.png'])

%% all F
figure(42)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_AllF)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Forage Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Preindust_Pac_biomes_AllF.png'])

% all P
figure(43)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_AllP)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Pelagic Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Preindust_Pac_biomes_AllP.png'])

% All D
figure(44)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(biome_AllD)))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Preindust_Pac_biomes_AllD.png'])
