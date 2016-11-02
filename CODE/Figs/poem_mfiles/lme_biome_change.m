% Changes to biomass in Pacific LMEs and biomes

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([cpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([cpath 'lme_mask_esm2m.mat']);
grid = csvread([cpath 'grid_csv.csv']);

% Biomes 1 = LC; 2 = ECCS; 3 = ECSS;

%% Preindust
load([dpath 'LME_preindust_' cfile '.mat']);
load([dpath 'Biomes_preindust_' cfile '.mat'],'biome_mbio');
load([dpath 'Preindust_1800-1850_means.mat']); %fcrit30

pre_g = lme_bio;
pre_gm2 = lme_mbio;

pre_biome = biome_mbio;

clear lme_bio lme_mbio biome_mbio

%% Hist prist
load([dpath 'LME_hist_pristine_' cfile '.mat'],'lme_bio00','lme_mbio00');
load([dpath 'Biomes_hist_pristine_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_hist_pristine_' cfile '.mat'],'sf_mean5000','sp_mean5000',...
    'sd_mean5000','mf_mean5000','mp_mean5000','md_mean5000','b_mean5000',...
    'lp_mean5000','ld_mean5000');

histp_g = lme_bio00;
histp_gm2 = lme_mbio00;

histp_biome = biome_mbio50;

clear lme_bio00 lme_mbio00 biome_mbio50

%% Global

%Mean global biomass in g/m2
all_pre(:,1) = sf_smean+sp_smean+sd_smean+mf_smean+mp_smean+md_smean+...
    lp_smean+ld_smean;
all_pre(:,2) = sf_smean+mf_smean; %all f
all_pre(:,3) = sp_smean+mp_smean+lp_smean; %all p
all_pre(:,4) = sd_smean+md_smean+ld_smean; %all d
all_pre(:,5) = sf_smean+sp_smean+sd_smean; %all s
all_pre(:,6) = mf_smean+mp_smean+md_smean; %all m
all_pre(:,7) = lp_smean+ld_smean; %all l

all_histp(:,1) = sf_mean5000+sp_mean5000+sd_mean5000+mf_mean5000+...
    mp_mean5000+md_mean5000+lp_mean5000+ld_mean5000;
all_histp(:,2) = sf_mean5000+mf_mean5000; %all f
all_histp(:,3) = sp_mean5000+mp_mean5000+lp_mean5000; %all p
all_histp(:,4) = sd_mean5000+md_mean5000+ld_mean5000; %all d
all_histp(:,5) = sf_mean5000+sp_mean5000+sd_mean5000; %all s
all_histp(:,6) = mf_mean5000+mp_mean5000+md_mean5000; %all m
all_histp(:,7) = lp_mean5000+ld_mean5000; %all l

%% NPac biomes

%Global biomes mean biomass in g/m2
pre_b(:,1) = nansum(pre_biome(:,1:8),2);
pre_b(:,2) = pre_biome(:,1) + pre_biome(:,4); %all f
pre_b(:,3) = pre_biome(:,2) + pre_biome(:,5) + pre_biome(:,7); %all p
pre_b(:,4) = pre_biome(:,3) + pre_biome(:,6) + pre_biome(:,8); %all d
pre_b(:,5) = nansum(pre_biome(:,1:3),2); %all s
pre_b(:,6) = nansum(pre_biome(:,4:6),2); %all m
pre_b(:,7) = nansum(pre_biome(:,7:8),2); %all l

histp_b(:,1) = nansum(histp_biome(:,1:8),2);
histp_b(:,2) = histp_biome(:,1) + histp_biome(:,4);
histp_b(:,3) = histp_biome(:,2) + histp_biome(:,5) + histp_biome(:,7);
histp_b(:,4) = histp_biome(:,3) + histp_biome(:,6) + histp_biome(:,8);
histp_b(:,5) = nansum(histp_biome(:,1:3),2);
histp_b(:,6) = nansum(histp_biome(:,4:6),2);
histp_b(:,7) = nansum(histp_biome(:,7:8),2);

%NPac biomes
lat=grid(:,3);
lon=grid(:,2);
id=find(lon<-180);
lon(id)=lon(id)+360;

%%
npid1 = find(geolon_t(:)>=-245 & geolon_t(:)<=-75);
npid2 = find(geolat_t(:)<=70 & geolat_t(:)>0);
npid=intersect(npid1,npid2);
test=zeros(size(geolon_t));
r=rand(size(npid));
test(npid)=r;

lonlim=[-255 -60];
latlim=[0 80];
figure(1)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','on','FLineWidth',1)
surfm(geolat_t,geolon_t,test)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hcb = colorbar('h');

%% N Pac LMEs

%LME areal biomass in g
pre(:,1) = nansum(pre_g(:,1:8),2);
pre(:,2) = pre_g(:,1) + pre_g(:,4); %all f
pre(:,3) = pre_g(:,2) + pre_g(:,5) + pre_g(:,7); %all p
pre(:,4) = pre_g(:,3) + pre_g(:,6) + pre_g(:,8); %all d
pre(:,5) = nansum(pre_g(:,1:3),2); %all s
pre(:,6) = nansum(pre_g(:,4:6),2); %all m
pre(:,7) = nansum(pre_g(:,7:8),2); %all l

histp(:,1) = nansum(histp_g(:,1:8),2);
histp(:,2) = histp_g(:,1) + histp_g(:,4);
histp(:,3) = histp_g(:,2) + histp_g(:,5) + histp_g(:,7);
histp(:,4) = histp_g(:,3) + histp_g(:,6) + histp_g(:,8);
histp(:,5) = nansum(histp_g(:,1:3),2);
histp(:,6) = nansum(histp_g(:,4:6),2);
histp(:,7) = nansum(histp_g(:,7:8),2);

%NPac LMEs = 47-57; 1-4; 10-11; 65

np = [1:4 10:11 47:57 65];

NPpre = pre(np,:);
NPhistp = histp(np,:);

diff_histp_pre = (NPhistp - NPpre) ./ NPpre;
diff_histp_pre_tot = (sum(NPhistp) - sum(NPpre)) ./ sum(NPpre);

% Figures
test(:,1)=histp_g(:,1);
test(:,2)=pre_g(:,2);
figure(100)
bar(test)
legend('histp','pre')

% By fish type
figure(1)
bar(diff_histp_pre_tot(2:4)*100)
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('Change in total fish abundance (g) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'NPac_Fish_biom_preindust_histprist.png'])

% By LME
figure(2)
bar(diff_histp_pre(:,1)*100)
set(gca,'XTick',1:length(np),'XTickLabel',np)
ylabel('Percent change')
xlabel('LME')
title('Change in total fish abundance (g) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'NPac_LME_biom_preindust_histprist.png'])






