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

load([dpath 'Means_hist_pristine_' cfile '.mat']);
load([cpath 'COBALT_biomes_last50yr.mat'],'biome_hist');
load([gpath 'hindcast_gridspec.mat'],'geolat_t','geolon_t');
grid = csvread([gpath 'grid_csv.csv']);

%% Add 50 yr time periods 
% 1950-2000
sp_mean(:,17)=sp_mean5000;
sf_mean(:,17)=sf_mean5000;
sd_mean(:,17)=sd_mean5000;
mp_mean(:,17)=mp_mean5000;
mf_mean(:,17)=mf_mean5000;
md_mean(:,17)=md_mean5000;
lp_mean(:,17)=lp_mean5000;
ld_mean(:,17)=ld_mean5000;
b_mean(:,17)=b_mean5000;
% 1955-2005
% sp_mean(:,18)=sp_mean5505;
% sf_mean(:,18)=sf_mean5505;
% sd_mean(:,18)=sd_mean5505;
% mp_mean(:,18)=mp_mean5505;
% mf_mean(:,18)=mf_mean5505;
% md_mean(:,18)=md_mean5505;
% lp_mean(:,18)=lp_mean5505;
% ld_mean(:,18)=ld_mean5505;
% b_mean(:,18)=b_mean5505;

%% Calc biome biomass avgs
gid = grid(:,1);

sf_biome_bio = NaN*ones(3,17);
sf_biome_mbio = sf_biome_bio;
sp_biome_mbio = sf_biome_bio;
sd_biome_mbio = sf_biome_bio;
mf_biome_mbio = sf_biome_bio;
mp_biome_mbio = sf_biome_bio;
md_biome_mbio = sf_biome_bio;
lp_biome_mbio = sf_biome_bio;
ld_biome_mbio = sf_biome_bio;
b_biome_mbio = sf_biome_bio;

for n = 1:17
    sf = NaN*ones(size(geolon_t));
    sp = sf;
    sd = sf;
    mf = sf;
    mp = sf;
    md = sf;
    lp = sf;
    ld = sf;
    b = sf;
    
    sf(gid) = sf_mean(:,n);
    sp(gid) = sp_mean(:,n);
    sd(gid) = sd_mean(:,n);
    mf(gid) = mf_mean(:,n);
    mp(gid) = mp_mean(:,n);
    md(gid) = md_mean(:,n);
    lp(gid) = lp_mean(:,n);
    ld(gid) = ld_mean(:,n);
    b(gid) = b_mean(:,n);
    
    for L=1:3
        lid = find(biome_hist==L);
        sf_biome_mbio(L,n) = nanmean(sf(lid));
        sp_biome_mbio(L,n) = nanmean(sp(lid));
        sd_biome_mbio(L,n) = nanmean(sd(lid));
        mf_biome_mbio(L,n) = nanmean(mf(lid));
        mp_biome_mbio(L,n) = nanmean(mp(lid));
        md_biome_mbio(L,n) = nanmean(md(lid));
        lp_biome_mbio(L,n) = nanmean(lp(lid));
        ld_biome_mbio(L,n) = nanmean(ld(lid));
        b_biome_mbio(L,n) = nanmean(b(lid));
    end
    
end

%%
biome_mbio50(:,1) = sf_biome_mbio(:,17);
biome_mbio50(:,2) = sp_biome_mbio(:,17);
biome_mbio50(:,3) = sd_biome_mbio(:,17);
biome_mbio50(:,4) = mf_biome_mbio(:,17);
biome_mbio50(:,5) = mp_biome_mbio(:,17);
biome_mbio50(:,6) = md_biome_mbio(:,17);
biome_mbio50(:,7) = lp_biome_mbio(:,17);
biome_mbio50(:,8) = ld_biome_mbio(:,17);
biome_mbio50(:,9) = b_biome_mbio(:,17);

% biome_mbio55(:,1) = sf_biome_mbio(:,18);
% biome_mbio55(:,2) = sp_biome_mbio(:,18);
% biome_mbio55(:,3) = sd_biome_mbio(:,18);
% biome_mbio55(:,4) = mf_biome_mbio(:,18);
% biome_mbio55(:,5) = mp_biome_mbio(:,18);
% biome_mbio55(:,6) = md_biome_mbio(:,18);
% biome_mbio55(:,7) = lp_biome_mbio(:,18);
% biome_mbio55(:,8) = ld_biome_mbio(:,18);
% biome_mbio55(:,9) = b_biome_mbio(:,18);

save([dpath 'Biomes_hist_pristine_' cfile '.mat'],'biome_mbio50',...
    'sf_biome_mbio','sp_biome_mbio','sd_biome_mbio','mf_biome_mbio','mp_biome_mbio',...
    'md_biome_mbio','b_biome_mbio','lp_biome_mbio','ld_biome_mbio'); %,'biome_mbio55');


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

%Pick which time period
yr=17; %1950-2000

for L=1:3
    lid = find(biome_hist==L);
    
    biome_sf(lid) = sf_biome_mbio(L,yr);
    biome_sp(lid) = sp_biome_mbio(L,yr);
    biome_sd(lid) = sd_biome_mbio(L,yr);
    biome_mf(lid) = mf_biome_mbio(L,yr);
    biome_mp(lid) = mp_biome_mbio(L,yr);
    biome_md(lid) = md_biome_mbio(L,yr);
    biome_lp(lid) = lp_biome_mbio(L,yr);
    biome_ld(lid) = ld_biome_mbio(L,yr);
    biome_b(lid) = b_biome_mbio(L,yr);
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
print('-dpng',[ppath 'Hist_pristine_biomes_All.png'])

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
print('-dpng',[ppath 'Hist_pristine_biomes_AllF.png'])

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
print('-dpng',[ppath 'Hist_pristine_biomes_AllP.png'])

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
print('-dpng',[ppath 'Hist_pristine_biomes_AllD.png'])

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
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                    %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_Pac_biomes_All.png'])

% all F
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
print('-dpng',[ppath 'Hist_pristine_Pac_biomes_AllF.png'])

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
print('-dpng',[ppath 'Hist_pristine_Pac_biomes_AllP.png'])

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
print('-dpng',[ppath 'Hist_pristine_Pac_biomes_AllD.png'])
