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


%% Preindust
load([dpath 'LME_preindust_' cfile '.mat']);

pre_g = lme_bio;
pre_gm2 = lme_mbio;

clear lme_bio lme_mbio

pre(:,1) = nansum(pre_g(:,1:8),2);
pre(:,2) = pre_g(:,1) + pre_g(:,4); %all f
pre(:,3) = pre_g(:,2) + pre_g(:,5) + pre_g(:,7); %all p
pre(:,4) = pre_g(:,3) + pre_g(:,6) + pre_g(:,8); %all d
pre(:,5) = nansum(pre_g(:,1:3),2); %all s
pre(:,6) = nansum(pre_g(:,4:6),2); %all m
pre(:,7) = nansum(pre_g(:,7:8),2); %all l

%% Hist prist
load([dpath 'LME_hist_pristine_' cfile '.mat'],'lme_bio05','lme_mbio05');

histp_g = lme_bio05;
histp_gm2 = lme_mbio05;

clear lme_bio05 lme_mbio05

histp(:,1) = nansum(histp_g(:,1:8),2);
histp(:,2) = histp_g(:,1) + histp_g(:,4);
histp(:,3) = histp_g(:,2) + histp_g(:,5) + histp_g(:,7);
histp(:,4) = histp_g(:,3) + histp_g(:,6) + histp_g(:,8);
histp(:,5) = nansum(histp_g(:,1:3),2);
histp(:,6) = nansum(histp_g(:,4:6),2);
histp(:,7) = nansum(histp_g(:,7:8),2);

%% N Pac
%NPac LMEs = 47-57; 1-4; 10-11; 65

np = [1:4 10:11 47:57 65];

NPpre = pre(np,:);
NPhistp = histp(np,:);

diff_histp_pre = (NPhistp - NPpre) ./ NPpre;
diff_histp_pre_tot = (sum(NPhistp) - sum(NPpre)) ./ sum(NPpre);

%% Figures

% By fish type
figure(1)
bar(diff_histp_pre_tot(2:4)*100)
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('Change in total fish abundance (g) from 1800-1850 to 1955-2005')

% By LME
figure(2)
bar(diff_histp_pre(:,1)*100)
set(gca,'XTick',1:length(np),'XTickLabel',np)
ylabel('Percent change')
xlabel('LME')
title('Change in total fish abundance (g) from 1800-1850 to 1955-2005')







