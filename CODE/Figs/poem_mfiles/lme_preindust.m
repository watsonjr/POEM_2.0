% Calc LME biomass of POEM
% Pre-industrial spinup at all locations
% 100 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Preindust_1800-1850_means.mat']);
load([cpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([cpath 'lme_mask_esm2m.mat']);
grid = csvread([cpath 'grid_csv.csv']);

%%
[loc,days]=size(SP.bio);
x=1:days;

%lyr=x((end-365+1):end);  %daily
lyr=x((end-12+1):end);   %monthly means

% sp_mean=nanmean(SP.bio(:,lyr),2);
% sf_mean=nanmean(SF.bio(:,lyr),2);
% sd_mean=nanmean(SD.bio(:,lyr),2);
% mp_mean=nanmean(MP.bio(:,lyr),2);
% mf_mean=nanmean(MF.bio(:,lyr),2);
% md_mean=nanmean(MD.bio(:,lyr),2);
% lp_mean=nanmean(LP.bio(:,lyr),2);
% ld_mean=nanmean(LD.bio(:,lyr),2);
% b_mean=nanmean(BENT.bio(:,lyr),2);

%% g/m2 --> total g
area = grid(:,5);
Asf_mean = sf_smean .* area;
Asp_mean = sp_smean .* area;
Asd_mean = sd_smean .* area;
Amf_mean = mf_smean .* area;
Amp_mean = mp_smean .* area;
Amd_mean = md_smean .* area;
Alp_mean = lp_smean .* area;
Ald_mean = ld_smean .* area;
Ab_mean = b_smean .* area;

%% Plots in space
lat = geolat_t;
lon = geolon_t;
%fix lon shift
id=find(lon(:)<-180);
lon(id)=lon(id)+360;

gid = grid(:,1);

sf = NaN*ones(size(lon));
sp = sf;
sd = sf;
mf = sf;
mp = sf;
md = sf;
lp = sf;
ld = sf;
b = sf;
Asf = sf;
Asp = sf;
Asd = sf;
Amf = sf;
Amp = sf;
Amd = sf;
Alp = sf;
Ald = sf;
Ab = sf;

sf(gid) = sf_smean;
sp(gid) = sp_smean;
sd(gid) = sd_smean;
mf(gid) = mf_smean;
mp(gid) = mp_smean;
md(gid) = md_smean;
lp(gid) = lp_smean;
ld(gid) = ld_smean;
b(gid) = b_smean;

Asf(gid) = Asf_mean;
Asp(gid) = Asp_mean;
Asd(gid) = Asd_mean;
Amf(gid) = Amf_mean;
Amp(gid) = Amp_mean;
Amd(gid) = Amd_mean;
Alp(gid) = Alp_mean;
Ald(gid) = Ald_mean;
Ab(gid) = Ab_mean;

%% Calc LMEs
tlme = lme_mask_esm2m';

figure(1)
m_proj('miller','lat',82);
m_pcolor(mask_lon_esm2m,mask_lat_esm2m,lme_mask_esm2m); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('LMEs')
colormap('prism')
colorbar('h')

figure(2)
m_proj('miller','lat',82);
m_pcolor(lon,lat,tlme); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('LMEs')
colormap('prism')
colorbar('h')

figure(3)
m_proj('miller','lat',82);
m_pcolor(lon',lat',lme_mask_esm2m); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('LMEs')
colormap('prism')
colorbar('h')

%%
lme_bio = NaN*ones(66,9);
lme_mbio = lme_bio;
for L=1:66
    lid = find(tlme==L);
    lme_bio(L,1) = nansum(Asf(lid));
    lme_bio(L,2) = nansum(Asp(lid));
    lme_bio(L,3) = nansum(Asd(lid));
    lme_bio(L,4) = nansum(Amf(lid));
    lme_bio(L,5) = nansum(Amp(lid));
    lme_bio(L,6) = nansum(Amd(lid));
    lme_bio(L,7) = nansum(Alp(lid));
    lme_bio(L,8) = nansum(Ald(lid));
    lme_bio(L,9) = nansum(Ab(lid));
    
    lme_mbio(L,1) = nanmean(sf(lid));
    lme_mbio(L,2) = nanmean(sp(lid));
    lme_mbio(L,3) = nanmean(sd(lid));
    lme_mbio(L,4) = nanmean(mf(lid));
    lme_mbio(L,5) = nanmean(mp(lid));
    lme_mbio(L,6) = nanmean(md(lid));
    lme_mbio(L,7) = nanmean(lp(lid));
    lme_mbio(L,8) = nanmean(ld(lid));
    lme_mbio(L,9) = nanmean(b(lid));
end
%%
save([dpath 'LME_preindust_Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05.mat'],...
    'lme_bio','lme_mbio');


%% Diff maps of all fish
All = nansum(lme_bio(:,1:8),2);
AllF = lme_bio(:,1) + lme_bio(:,4);
AllP = lme_bio(:,2) + lme_bio(:,5) + lme_bio(:,7);
AllD = lme_bio(:,3) + lme_bio(:,6) + lme_bio(:,8);
AllS = nansum(lme_bio(:,1:3));
AllM = nansum(lme_bio(:,4:6));
AllL = nansum(lme_bio(:,7:8));
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);

%grams to tonnes/metric tons 1g = 1e-6MT
nansum(All) % = 1.7206e+14
nansum(All) * 1e-6 % = 1.7206e+08 = 0.17e9 MT
% My max fish size if 2.5kg technically rep fish 80g to 80kg
% Jennings 10-10^6 g = 0.34-26.12 10^9 MT
% Jennings 10^2-10^4 g = 0.09-8.89 10^9 MT

%% ALL
lname = 1:66;
[sort_bio,ix] = sort(All);
slname = lname(ix);

%%
figure(4)
plot(log10(sort_bio),1:66,'.k','MarkerSize',15); hold on;
set(gca,'YTickLabel',[])
for i=1:66
    text(11,i,num2str(ix(i)))
end
title('log10 mean biomass All Fishes (g)')
ylabel('LME')
stamp(cfile)
print('-dpng',[ppath 'Preindust_LME_dot_areal_sum.png'])

%%
mAll = nansum(lme_mbio(:,1:8),2);

[sort_mbio,ix] = sort(mAll,'descend');
mslname = lname(ix);

figure(5)
plot(log10(sort_mbio),1:66,'.k','MarkerSize',15); hold on;
set(gca,'YTickLabel',[])
for i=1:66
    text(-1.4,i,num2str(ix(i)))
end
title('log10 mean biomass All Fishes (g m^-^2)')
ylabel('LME')
stamp(cfile)
print('-dpng',[ppath 'Preindust_LME_dot_mean.png'])


%% Figures
lme_sf = NaN*ones(360,200);
lme_sp = lme_sf;
lme_sd = lme_sf;
lme_mf = lme_sf;
lme_mp = lme_sf;
lme_md = lme_sf;
lme_lp = lme_sf;
lme_ld = lme_sf;
lme_b = lme_sf;

for L=1:66
    lid = find(tlme==L);
    
    lme_sf(lid) = lme_mbio(L,1);
    lme_sp(lid) = lme_mbio(L,2);
    lme_sd(lid) = lme_mbio(L,3);
    lme_mf(lid) = lme_mbio(L,4);
    lme_mp(lid) = lme_mbio(L,5);
    lme_md(lid) = lme_mbio(L,6);
    lme_lp(lid) = lme_mbio(L,7);
    lme_ld(lid) = lme_mbio(L,8);
    lme_b(lid) = lme_mbio(L,9);
end

lme_All = lme_sf+lme_sp+lme_sd+lme_mf+lme_mp+lme_md+lme_lp+lme_ld;
lme_AllF = lme_sf+lme_mf;
lme_AllP = lme_sp+lme_mp+lme_lp;
lme_AllD = lme_sd+lme_md+lme_ld;

% plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

%% ALL
figure(31)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_All)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
axesmui
stamp(cfile)
print('-dpng',[ppath 'Preindust_LME_All.png'])

% all F
figure(32)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllF)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Forage Fishes (g m^-^2)')
axesmui
stamp(cfile)
print('-dpng',[ppath 'Preindust_LME_AllF.png'])

% all P
figure(33)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllP)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Pelagic Fishes (g m^-^2)')
axesmui
stamp(cfile)
print('-dpng',[ppath 'Preindust_LME_AllP.png'])

% All D
figure(34)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllD)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Demersal Fishes (g m^-^2)')
axesmui
stamp(cfile)
print('-dpng',[ppath 'Preindust_LME_AllD.png'])

%% NPac only
% ALL
figure(10)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_All)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Preindust_Pac_LME_All.png'])

% all F
figure(11)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllF)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Forage Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Preindust_Pac_LME_AllF.png'])

% all P
figure(12)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllP)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Pelagic Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Preindust_Pac_LME_AllP.png'])

% All D
figure(13)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllD)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Preindust_Pac_LME_AllD.png'])

