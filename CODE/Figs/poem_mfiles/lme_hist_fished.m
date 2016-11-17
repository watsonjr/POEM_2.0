% Calc LME biomass of POEM
% Historic time period (1861-2005) at all locations
% 145 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Means_hist_fished_' cfile '.mat']);
load([cpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([cpath 'lme_mask_esm2m.mat']);
grid = csvread([cpath 'grid_csv.csv']);

%%
sf_mean(:,17:18) = [sf_mean5000, sf_mean5505];
sp_mean(:,17:18) = [sp_mean5000, sp_mean5505];
sd_mean(:,17:18) = [sd_mean5000, sd_mean5505];
mf_mean(:,17:18) = [mf_mean5000, mf_mean5505];
mp_mean(:,17:18) = [mp_mean5000, mp_mean5505];
md_mean(:,17:18) = [md_mean5000, md_mean5505];
lp_mean(:,17:18) = [lp_mean5000, lp_mean5505];
ld_mean(:,17:18) = [ld_mean5000, ld_mean5505];
b_mean(:,17:18) = [b_mean5000, b_mean5505];

mf_mcatch(:,17:18) = [mf_mcatch5000, mf_mcatch5505];
mp_mcatch(:,17:18) = [mp_mcatch5000, mp_mcatch5505];
md_mcatch(:,17:18) = [md_mcatch5000, md_mcatch5505];
lp_mcatch(:,17:18) = [lp_mcatch5000, lp_mcatch5505];
ld_mcatch(:,17:18) = [ld_mcatch5000, ld_mcatch5505];

%% g/m2 --> total g
area = grid(:,5);
Asf_mean = sf_mean .* repmat(area,1,18);
Asp_mean = sp_mean .* repmat(area,1,18);
Asd_mean = sd_mean .* repmat(area,1,18);
Amf_mean = mf_mean .* repmat(area,1,18);
Amp_mean = mp_mean .* repmat(area,1,18);
Amd_mean = md_mean .* repmat(area,1,18);
Alp_mean = lp_mean .* repmat(area,1,18);
Ald_mean = ld_mean .* repmat(area,1,18);
Ab_mean = b_mean .* repmat(area,1,18);

Amf_mcatch = mf_mcatch .* repmat(area,1,18) * 365; %mean fish catch per yr
Amp_mcatch = mp_mcatch .* repmat(area,1,18) * 365;
Amd_mcatch = md_mcatch .* repmat(area,1,18) * 365;
Alp_mcatch = lp_mcatch .* repmat(area,1,18) * 365;
Ald_mcatch = ld_mcatch .* repmat(area,1,18) * 365;

%% Calc LMEs
lat = geolat_t;
lon = geolon_t;
%fix lon shift
id=find(lon(:)<-180);
lon(id)=lon(id)+360;

gid = grid(:,1);

tlme = lme_mask_esm2m';

sf_lme_bio = NaN*ones(66,18);
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
mf_lme_mcatch = sf_lme_bio;
mp_lme_mcatch = sf_lme_bio;
md_lme_mcatch = sf_lme_bio;
lp_lme_mcatch = sf_lme_bio;
ld_lme_mcatch = sf_lme_bio;

for n = 1:18
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
    Cmf = sf;
    Cmp = sf;
    Cmd = sf;
    Clp = sf;
    Cld = sf;
    
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
    
    Cmf(gid) = Amf_mcatch(:,n);
    Cmp(gid) = Amp_mcatch(:,n);
    Cmd(gid) = Amd_mcatch(:,n);
    Clp(gid) = Alp_mcatch(:,n);
    Cld(gid) = Ald_mcatch(:,n);
    
    for L=1:66
        lid = find(tlme==L);
        %total biomass g
        sf_lme_bio(L,n) = nansum(Asf(lid));
        sp_lme_bio(L,n) = nansum(Asp(lid));
        sd_lme_bio(L,n) = nansum(Asd(lid));
        mf_lme_bio(L,n) = nansum(Amf(lid));
        mp_lme_bio(L,n) = nansum(Amp(lid));
        md_lme_bio(L,n) = nansum(Amd(lid));
        lp_lme_bio(L,n) = nansum(Alp(lid));
        ld_lme_bio(L,n) = nansum(Ald(lid));
        b_lme_bio(L,n) = nansum(Ab(lid));
        %mean biomass g/m2
        sf_lme_mbio(L,n) = nanmean(sf(lid));
        sp_lme_mbio(L,n) = nanmean(sp(lid));
        sd_lme_mbio(L,n) = nanmean(sd(lid));
        mf_lme_mbio(L,n) = nanmean(mf(lid));
        mp_lme_mbio(L,n) = nanmean(mp(lid));
        md_lme_mbio(L,n) = nanmean(md(lid));
        lp_lme_mbio(L,n) = nanmean(lp(lid));
        ld_lme_mbio(L,n) = nanmean(ld(lid));
        b_lme_mbio(L,n) = nanmean(b(lid));
        %total catch g
        mf_lme_mcatch(L,n) = nansum(Cmf(lid));
        mp_lme_mcatch(L,n) = nansum(Cmp(lid));
        md_lme_mcatch(L,n) = nansum(Cmd(lid));
        lp_lme_mcatch(L,n) = nansum(Clp(lid));
        ld_lme_mcatch(L,n) = nansum(Cld(lid));
    end
    
end
%% Comps by LME to J&C15

%Pick which time period
yr=17; %1950-2000
%yr=18; %1955-2005

All = sf_lme_bio(:,yr)+sp_lme_bio(:,yr)+sd_lme_bio(:,yr)+mf_lme_bio(:,yr)+...
    mp_lme_bio(:,yr)+md_lme_bio(:,yr)+lp_lme_bio(:,yr)+ld_lme_bio(:,yr);
ML = mf_lme_bio(:,yr)+mp_lme_bio(:,yr)+md_lme_bio(:,yr)+lp_lme_bio(:,yr)+ld_lme_bio(:,yr);

% JUST LMEs, NOT WHOLE GLOBE!!!
%grams to tonnes/metric tons 1g = 1e-6MT
nansum(All) % = 1.0424e+14
nansum(All) * 1e-6 % = 1.0424e+08 = 0.15e9 MT
% My max fish size if 2.5kg technically rep fish 80g to 80kg
% Jennings 10-10^6 g = 0.34-26.12 10^9 MT
% Jennings 10^2-10^4 g = 0.09-8.89 10^9 MT
% My med & lrg 10^2-10^4 = 0.12819 10^9 MT
% If I use their wet weight conversion, which is 3x mine, 0.31271e9 MT

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
print('-dpng',[ppath 'Hist_fished_LME_dot_areal_sum.png'])

%%
mAll = sf_lme_mbio(:,yr)+sp_lme_mbio(:,yr)+sd_lme_mbio(:,yr)+mf_lme_mbio(:,yr)+...
    mp_lme_mbio(:,yr)+md_lme_mbio(:,yr)+lp_lme_mbio(:,yr)+ld_lme_mbio(:,yr);

[sort_mbio,ix] = sort(mAll,'descend');
mslname = lname(ix);

figure(5)
plot(log10(sort_mbio),1:66,'.k','MarkerSize',15); hold on;
set(gca,'YTickLabel',[])
for i=1:66
    text(-1.35,i,num2str(ix(i)))
end
title('log10 mean biomass All Fishes (g m^-^2)')
ylabel('LME')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_dot_mean.png'])

T = table(ix,sort_mbio,'VariableNames',{'LME','mbiomgm2'});
writetable(T,[dpath 'LME_mbio_hist_fished_' cfile '.csv']);

%% save 1950-2000 and 1955-2005 for comps
lme_bio00(:,1) = sf_lme_bio(:,17);
lme_bio00(:,2) = sp_lme_bio(:,17);
lme_bio00(:,3) = sd_lme_bio(:,17);
lme_bio00(:,4) = mf_lme_bio(:,17);
lme_bio00(:,5) = mp_lme_bio(:,17);
lme_bio00(:,6) = md_lme_bio(:,17);
lme_bio00(:,7) = lp_lme_bio(:,17);
lme_bio00(:,8) = ld_lme_bio(:,17);
lme_bio00(:,9) = b_lme_bio(:,17);

lme_mbio00(:,1) = sf_lme_mbio(:,17);
lme_mbio00(:,2) = sp_lme_mbio(:,17);
lme_mbio00(:,3) = sd_lme_mbio(:,17);
lme_mbio00(:,4) = mf_lme_mbio(:,17);
lme_mbio00(:,5) = mp_lme_mbio(:,17);
lme_mbio00(:,6) = md_lme_mbio(:,17);
lme_mbio00(:,7) = lp_lme_mbio(:,17);
lme_mbio00(:,8) = ld_lme_mbio(:,17);
lme_mbio00(:,9) = b_lme_mbio(:,17);

lme_mcatch00(:,1) = mf_lme_mcatch(:,17);
lme_mcatch00(:,2) = mp_lme_mcatch(:,17);
lme_mcatch00(:,3) = md_lme_mcatch(:,17);
lme_mcatch00(:,4) = lp_lme_mcatch(:,17);
lme_mcatch00(:,5) = ld_lme_mcatch(:,17);

lme_bio05(:,1) = sf_lme_bio(:,18);
lme_bio05(:,2) = sp_lme_bio(:,18);
lme_bio05(:,3) = sd_lme_bio(:,18);
lme_bio05(:,4) = mf_lme_bio(:,18);
lme_bio05(:,5) = mp_lme_bio(:,18);
lme_bio05(:,6) = md_lme_bio(:,18);
lme_bio05(:,7) = lp_lme_bio(:,18);
lme_bio05(:,8) = ld_lme_bio(:,18);
lme_bio05(:,9) = b_lme_bio(:,18);

lme_mbio05(:,1) = sf_lme_mbio(:,18);
lme_mbio05(:,2) = sp_lme_mbio(:,18);
lme_mbio05(:,3) = sd_lme_mbio(:,18);
lme_mbio05(:,4) = mf_lme_mbio(:,18);
lme_mbio05(:,5) = mp_lme_mbio(:,18);
lme_mbio05(:,6) = md_lme_mbio(:,18);
lme_mbio05(:,7) = lp_lme_mbio(:,18);
lme_mbio05(:,8) = ld_lme_mbio(:,18);
lme_mbio05(:,9) = b_lme_mbio(:,18);

lme_mcatch05(:,1) = mf_lme_mcatch(:,18);
lme_mcatch05(:,2) = mp_lme_mcatch(:,18);
lme_mcatch05(:,3) = md_lme_mcatch(:,18);
lme_mcatch05(:,4) = lp_lme_mcatch(:,18);
lme_mcatch05(:,5) = ld_lme_mcatch(:,18);


%%
save([dpath 'LME_hist_fished_' cfile '.mat'],...
    'sf_lme_bio','sp_lme_bio','sd_lme_bio','mf_lme_bio','mp_lme_bio','md_lme_bio',...
    'b_lme_bio','lp_lme_bio','ld_lme_bio','sf_lme_mbio','sp_lme_mbio','sd_lme_mbio',...
    'mf_lme_mbio','mp_lme_mbio','md_lme_mbio','b_lme_mbio','lp_lme_mbio','ld_lme_mbio',...
    'lme_bio00','lme_mbio00','lme_bio05','lme_mbio05','lme_mcatch00','lme_mcatch05',...
    'mf_lme_mcatch','mp_lme_mcatch','md_lme_mcatch','lp_lme_mcatch','ld_lme_mcatch');


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
clme_mf = lme_sf;
clme_mp = lme_sf;
clme_md = lme_sf;
clme_lp = lme_sf;
clme_ld = lme_sf;

for L=1:66
    lid = find(tlme==L);
    
    lme_sf(lid) = sf_lme_mbio(L,yr);
    lme_sp(lid) = sp_lme_mbio(L,yr);
    lme_sd(lid) = sd_lme_mbio(L,yr);
    lme_mf(lid) = mf_lme_mbio(L,yr);
    lme_mp(lid) = mp_lme_mbio(L,yr);
    lme_md(lid) = md_lme_mbio(L,yr);
    lme_lp(lid) = lp_lme_mbio(L,yr);
    lme_ld(lid) = ld_lme_mbio(L,yr);
    lme_b(lid) = b_lme_mbio(L,yr);
    
    clme_mf(lid) = mf_lme_mcatch(L,yr);
    clme_mp(lid) = mp_lme_mcatch(L,yr);
    clme_md(lid) = md_lme_mcatch(L,yr);
    clme_lp(lid) = lp_lme_mcatch(L,yr);
    clme_ld(lid) = ld_lme_mcatch(L,yr);
end

lme_All = lme_sf+lme_sp+lme_sd+lme_mf+lme_mp+lme_md+lme_lp+lme_ld;
lme_AllF = lme_sf+lme_mf;
lme_AllP = lme_sp+lme_mp+lme_lp;
lme_AllD = lme_sd+lme_md+lme_ld;

clme_All = clme_mf+clme_mp+clme_md+clme_lp+clme_ld;
clme_AllF = clme_mf;
clme_AllP = clme_mp+clme_lp;
clme_AllD = clme_md+clme_ld;
clme_AllM = clme_mf+clme_mp+clme_md;
clme_AllL = clme_lp+clme_ld;

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
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_All.png'])

% all F
figure(32)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllF)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 log10 mean biomass Forage Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_AllF.png'])

% all P
figure(33)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllP)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 log10 mean biomass Pelagic Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_AllP.png'])

% All D
figure(34)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllD)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_AllD.png'])

%% Catch
% all
figure(51)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_All*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([4.5 6.5]);
hcb = colorbar('h');
ylim(hcb,[4.5 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 mean log10 total annual catch (MT) All Fishes')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_catch_All.png'])

%% all F
figure(52)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllF*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 6.5]);
hcb = colorbar('h');
ylim(hcb,[3 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 mean log10 total annual catch (MT) Forage Fishes')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_catch_AllF.png'])

% all P
figure(53)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllP*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 6.5]);
hcb = colorbar('h');
ylim(hcb,[3 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 mean log10 total annual catch (MT) Pelagic Fishes')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_catch_AllP.png'])

% All D
figure(54)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllD*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 6.5]);
hcb = colorbar('h');
ylim(hcb,[3 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 mean log10 total annual catch (MT) Demersal Fishes')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_catch_AllD.png'])

% all M
figure(55)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllM*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 6.5]);
hcb = colorbar('h');
ylim(hcb,[3 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 mean log10 total annual catch (MT) Medium Fishes')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_catch_AllM.png'])

% all L
figure(56)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(clme_AllL*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([3 6.5]);
hcb = colorbar('h');
ylim(hcb,[3 6.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 mean log10 total annual catch (MT) Large Fishes')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_LME_catch_AllL.png'])

%% NPac only
% ALL
figure(10)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_All)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_Pac_LME_All.png'])

% all F
figure(11)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllF)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 log10 mean biomass Forage Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_Pac_LME_AllF.png'])

% all P
figure(12)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllP)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 log10 mean biomass Pelagic Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_Pac_LME_AllP.png'])

% All D
figure(13)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(lme_AllD)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('1951-2000 log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Hist_fished_Pac_LME_AllD.png'])
