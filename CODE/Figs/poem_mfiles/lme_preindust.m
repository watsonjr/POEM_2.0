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

lyr=x((end-365+1):end); %daily
%lyr=x((end-12+1):end)   %monthly means

sp_mean=mean(SP.bio(:,lyr),2);
sf_mean=mean(SF.bio(:,lyr),2);
sd_mean=mean(SD.bio(:,lyr),2);
mp_mean=mean(MP.bio(:,lyr),2);
mf_mean=mean(MF.bio(:,lyr),2);
md_mean=mean(MD.bio(:,lyr),2);
lp_mean=mean(LP.bio(:,lyr),2);
ld_mean=mean(LD.bio(:,lyr),2);
b_mean=mean(BENT.bio(:,lyr),2);

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
for L=1:66
    lid = find(tlme==L);
    lme_bio(L,1) = sum(Asf(lid));
    lme_bio(L,2) = sum(Asp(lid));
    lme_bio(L,3) = sum(Asd(lid));
    lme_bio(L,4) = sum(Amf(lid));
    lme_bio(L,5) = sum(Amp(lid));
    lme_bio(L,6) = sum(Amd(lid));
    lme_bio(L,7) = sum(Alp(lid));
    lme_bio(L,8) = sum(Ald(lid));
    lme_bio(L,9) = sum(Ab(lid));
end


%% Diff maps of all fish
All = sum(lme_bio(:,1:8),2);
AllF = lme_bio(:,1) + lme_bio(:,4);
AllP = lme_bio(:,2) + lme_bio(:,5) + lme_bio(:,7);
AllD = lme_bio(:,3) + lme_bio(:,6) + lme_bio(:,8);
AllS = sum(lme_bio(:,1:3));
AllM = sum(lme_bio(:,4:6));
AllL = sum(lme_bio(:,7:8));
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);

%grams to tonnes/metric tons 1g = 1e-6MT
sum(All)
sum(All) * 1e-6
% My max fish size if 2.5kg technically rep fish 80g to 80kg
% Jennings 10-10^6 g = 0.34-26.12 MT
% Jennings 10^2-10^4 g = 0.09-8.89 MT

%% ALL
lname = 1:66;
[sort_bio,ix] = sort(All);
slname = lname(ix);

%%
figure(4)
plot(sort_bio,1:66,'.k','MarkerSize',15); hold on;
set(gca,'YTickLabel',[])
for i=1:66
    text(-1e11,i,num2str(ix(i)))
end
title('mean biomass All Fishes (g)')
%stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_All.png'])

