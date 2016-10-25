% Calc LME biomass of POEM
% Pre-industrial spinup at all locations
% 100 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Data_spinup_pristine_' cfile '.mat']);
load([cpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([cpath 'lme_mask_esm2m.mat']);
grid = csvread([cpath 'grid_csv.csv']);

%%
[loc,days]=size(SP.bio);
x=1:days;

%lyr=x((end-365+1):end);  %daily
lyr=x((end-12+1):end);   %monthly means

sp_mean=nanmean(SP.bio(:,lyr),2);
sf_mean=nanmean(SF.bio(:,lyr),2);
sd_mean=nanmean(SD.bio(:,lyr),2);
mp_mean=nanmean(MP.bio(:,lyr),2);
mf_mean=nanmean(MF.bio(:,lyr),2);
md_mean=nanmean(MD.bio(:,lyr),2);
lp_mean=nanmean(LP.bio(:,lyr),2);
ld_mean=nanmean(LD.bio(:,lyr),2);
b_mean=nanmean(BENT.bio(:,lyr),2);

%% g/m2 --> total g
area = grid(:,5);
Asf_mean = sf_mean .* area;
Asp_mean = sp_mean .* area;
Asd_mean = sd_mean .* area;
Amf_mean = mf_mean .* area;
Amp_mean = mp_mean .* area;
Amd_mean = md_mean .* area;
Alp_mean = lp_mean .* area;
Ald_mean = ld_mean .* area;
Ab_mean = b_mean .* area;

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

sf(gid) = sf_mean;
sp(gid) = sp_mean;
sd(gid) = sd_mean;
mf(gid) = mf_mean;
mp(gid) = mp_mean;
md(gid) = md_mean;
lp(gid) = lp_mean;
ld(gid) = ld_mean;
b(gid) = b_mean;

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
nansum(All) % = 1.7102e+14
nansum(All) * 1e-6 % = 1.7102e+08 = 0.17e9 MT
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

