% Calc LME biomass of POEM
% Historic time period (1861-2005) at all locations
% 145 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Means_hist_pristine_' cfile '.mat']);
load([cpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([cpath 'lme_mask_esm2m.mat']);
grid = csvread([cpath 'grid_csv.csv']);

%% g/m2 --> total g
area = grid(:,5);
Asf_mean = sf_mean .* repmat(area,1,16);
Asp_mean = sp_mean .* repmat(area,1,16);
Asd_mean = sd_mean .* repmat(area,1,16);
Amf_mean = mf_mean .* repmat(area,1,16);
Amp_mean = mp_mean .* repmat(area,1,16);
Amd_mean = md_mean .* repmat(area,1,16);
Alp_mean = lp_mean .* repmat(area,1,16);
Ald_mean = ld_mean .* repmat(area,1,16);
Ab_mean = b_mean .* repmat(area,1,16);

%% Calc LMEs
lat = geolat_t;
lon = geolon_t;
%fix lon shift
id=find(lon(:)<-180);
lon(id)=lon(id)+360;

gid = grid(:,1);

tlme = lme_mask_esm2m';

sf_lme_bio = NaN*ones(66,16);
sp_lme_bio = sf_lme_bio;
sd_lme_bio = sf_lme_bio;
mf_lme_bio = sf_lme_bio;
mp_lme_bio = sf_lme_bio;
md_lme_bio = sf_lme_bio;
lp_lme_bio = sf_lme_bio;
ld_lme_bio = sf_lme_bio;
b_lme_bio = sf_lme_bio;
sf_lme_mbio = sf_lme_bio;
sp_lme_mbio = sf_lme_bio;
sd_lme_mbio = sf_lme_bio;
mf_lme_mbio = sf_lme_bio;
mp_lme_mbio = sf_lme_bio;
md_lme_mbio = sf_lme_bio;
lp_lme_mbio = sf_lme_bio;
ld_lme_mbio = sf_lme_bio;
b_lme_mbio = sf_lme_bio;

for n = 1:16
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
    
    sf(gid) = sf_mean(:,n);
    sp(gid) = sp_mean(:,n);
    sd(gid) = sd_mean(:,n);
    mf(gid) = mf_mean(:,n);
    mp(gid) = mp_mean(:,n);
    md(gid) = md_mean(:,n);
    lp(gid) = lp_mean(:,n);
    ld(gid) = ld_mean(:,n);
    b(gid) = b_mean(:,n);
    
    Asf(gid) = Asf_mean(:,n);
    Asp(gid) = Asp_mean(:,n);
    Asd(gid) = Asd_mean(:,n);
    Amf(gid) = Amf_mean(:,n);
    Amp(gid) = Amp_mean(:,n);
    Amd(gid) = Amd_mean(:,n);
    Alp(gid) = Alp_mean(:,n);
    Ald(gid) = Ald_mean(:,n);
    Ab(gid) = Ab_mean(:,n);
    
    
    
    for L=1:66
        lid = find(tlme==L);
        sf_lme_bio(L,n) = nansum(Asf(lid));
        sp_lme_bio(L,n) = nansum(Asp(lid));
        sd_lme_bio(L,n) = nansum(Asd(lid));
        mf_lme_bio(L,n) = nansum(Amf(lid));
        mp_lme_bio(L,n) = nansum(Amp(lid));
        md_lme_bio(L,n) = nansum(Amd(lid));
        lp_lme_bio(L,n) = nansum(Alp(lid));
        ld_lme_bio(L,n) = nansum(Ald(lid));
        b_lme_bio(L,n) = nansum(Ab(lid));
        
        sf_lme_mbio(L,n) = nanmean(sf(lid));
        sp_lme_mbio(L,n) = nanmean(sp(lid));
        sd_lme_mbio(L,n) = nanmean(sd(lid));
        mf_lme_mbio(L,n) = nanmean(mf(lid));
        mp_lme_mbio(L,n) = nanmean(mp(lid));
        md_lme_mbio(L,n) = nanmean(md(lid));
        lp_lme_mbio(L,n) = nanmean(lp(lid));
        ld_lme_mbio(L,n) = nanmean(ld(lid));
        b_lme_mbio(L,n) = nanmean(b(lid));
    end
    
end
%% Comps by LME to J&C15
All = sf_lme_bio(:,1)+sp_lme_bio(:,1)+sd_lme_bio(:,1)+mf_lme_bio(:,1)+...
    mp_lme_bio(:,1)+md_lme_bio(:,1)+lp_lme_bio(:,1)+ld_lme_bio(:,1);
ML = mf_lme_bio(:,1)+mp_lme_bio(:,1)+md_lme_bio(:,1)+lp_lme_bio(:,1)+ld_lme_bio(:,1);

%grams to tonnes/metric tons 1g = 1e-6MT
nansum(All) % = 1.7118e+14
nansum(All) * 1e-6 % = 1.7118e+08 = 0.17e9 MT
% My max fish size if 2.5kg technically rep fish 80g to 80kg
% Jennings 10-10^6 g = 0.34-26.12 10^9 MT
% Jennings 10^2-10^4 g = 0.09-8.89 10^9 MT
% My med & lrg 10^2-10^4 = 0.14299 10^9 MT

lname = 1:66;
[sort_bio,ix] = sort(All);
slname = lname(ix);

figure(4)
plot(log10(sort_bio),1:66,'.k','MarkerSize',15); hold on;
set(gca,'YTickLabel',[])
for i=1:66
    text(11,i,num2str(ix(i)))
end
title('log10 mean biomass All Fishes (g)')
ylabel('LME')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_LME_dot_areal_sum.png'])

%%
mAll = nansum(lme_mbio(:,1:8),2);

[sort_mbio,ix] = sort(mAll,'descend');
mslname = lname(ix);

figure(5)
plot(log10(sort_mbio),1:66,'.k','MarkerSize',15); hold on;
set(gca,'YTickLabel',[])
for i=1:66
    text(-1.15,i,num2str(ix(i)))
end
title('log10 mean biomass All Fishes (g m^-^2)')
ylabel('LME')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_LME_dot_mean.png'])

%%
save([fpath 'LME_hist_pristine_Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05.mat'],...
    'sf_lme_bio','sp_lme_bio','sd_lme_bio','mf_lme_bio','mp_lme_bio','md_lme_bio',...
    'b_lme_bio','lp_lme_bio','ld_lme_bio','sf_lme_mbio','sp_lme_mbio','sd_lme_mbio',...
    'mf_lme_mbio','mp_lme_mbio','md_lme_mbio','b_lme_mbio','lp_lme_mbio','ld_lme_mbio');


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
    
    lme_sf(lid) = sf_lme_mbio(L,1);
    lme_sp(lid) = sp_lme_mbio(L,1);
    lme_sd(lid) = sd_lme_mbio(L,1);
    lme_mf(lid) = mf_lme_mbio(L,1);
    lme_mp(lid) = mp_lme_mbio(L,1);
    lme_md(lid) = md_lme_mbio(L,1);
    lme_lp(lid) = lp_lme_mbio(L,1);
    lme_ld(lid) = ld_lme_mbio(L,1);
    lme_b(lid) = b_lme_mbio(L,1);
end

lme_All = lme_sf+lme_sp+lme_sd+lme_mf+lme_mp+lme_md+lme_lp+lme_ld;
lme_AllF = lme_sf+lme_mf;
lme_AllP = lme_sp+lme_mp+lme_lp;
lme_AllD = lme_sd+lme_md+lme_ld;

% plot info
GEOLON_T = double(GEOLON_T);
GEOLAT_T = double(GEOLAT_T);
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
surfm(GEOLAT_T,GEOLON_T,real(log10(lme_All)))
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
print('-dpng',[ppath 'Hist_pristine_LME_All.png'])

% all F
figure(32)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(GEOLAT_T,GEOLON_T,real(log10(lme_AllF)))
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
print('-dpng',[ppath 'Hist_pristine_LME_AllF.png'])

% all P
figure(33)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(GEOLAT_T,GEOLON_T,real(log10(lme_AllP)))
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
print('-dpng',[ppath 'Hist_pristine_LME_AllP.png'])

% All D
figure(34)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(GEOLAT_T,GEOLON_T,real(log10(lme_AllD)))
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
print('-dpng',[ppath 'Hist_pristine_LME_AllD.png'])

%% NPac only
% ALL
figure(10)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(GEOLAT_T,GEOLON_T,real(log10(lme_All)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_Pac_LME_All.png'])

% all F
figure(11)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(GEOLAT_T,GEOLON_T,real(log10(lme_AllF)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Forage Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_Pac_LME_AllF.png'])

% all P
figure(12)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(GEOLAT_T,GEOLON_T,real(log10(lme_AllP)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Pelagic Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_Pac_LME_AllP.png'])

% All D
figure(13)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(GEOLAT_T,GEOLON_T,real(log10(lme_AllD)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1.5 1.5]);
hcb = colorbar('h');
ylim(hcb,[-1.5 1.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_pristine_Pac_LME_AllD.png'])

