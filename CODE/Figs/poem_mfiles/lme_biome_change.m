% Changes to biomass in Pacific LMEs and biomes

clear all
close all

gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([gpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([gpath 'lme_mask_esm2m.mat']);
load([gpath 'NPid_esm2m.mat']);
load([cpath 'COBALT_biomes.mat']);
grid = csvread([gpath 'grid_csv.csv']);

%Colormap
load('cmap_ppt_angles.mat')
cmap1(1,:)=[1 1 1];
cmap1(2,:)=cmap_ppt(1,:);
cmap1(3,:)=cmap_ppt(3,:);
cmap1(4,:)=cmap_ppt(5,:);

cmap3(1,:)=[0 0 0];
cmap3(2,:)=cmap_ppt(1,:);
cmap3(3,:)=cmap_ppt(3,:);
cmap3(4,:)=cmap_ppt(5,:);

cmap2(1,:)=cmap_ppt(1,:);
cmap2(2,:)=cmap_ppt(3,:);
cmap2(3,:)=cmap_ppt(5,:);
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

%% Preindust
load([dpath 'LME_preindust_' cfile '.mat']);
load([dpath 'Biomes_preindust_' cfile '.mat'],'biome_mbio');
load([dpath 'Preindust_1800-1850_means.mat']); %fcrit30

pre_g = lme_bio;
pre_gm2 = lme_mbio;

pre_biome = biome_mbio;

%Mean global biomass in g/m2
all_pre(:,1) = sf_smean+sp_smean+sd_smean+mf_smean+mp_smean+md_smean+...
    lp_smean+ld_smean;
all_pre(:,2) = sf_smean+mf_smean; %all f
all_pre(:,3) = sp_smean+mp_smean+lp_smean; %all p
all_pre(:,4) = sd_smean+md_smean+ld_smean; %all d
all_pre(:,5) = sf_smean+sp_smean+sd_smean; %all s
all_pre(:,6) = mf_smean+mp_smean+md_smean; %all m
all_pre(:,7) = lp_smean+ld_smean; %all l

%Mean total global biomass in g
tot_all_pre = all_pre .* repmat(grid(:,5),1,7);

clear lme_bio lme_mbio biome_mbio sf_smean sp_smean sd_smean
clear mf_smean mp_smean md_smean b_smean lp_smean ld_smean

%% Hist prist
load([dpath 'LME_hist_pristine_' cfile '.mat'],'lme_bio00','lme_mbio00');
load([dpath 'Biomes_hist_pristine_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_hist_pristine_' cfile '.mat'],'sf_mean5000','sp_mean5000',...
    'sd_mean5000','mf_mean5000','mp_mean5000','md_mean5000','b_mean5000',...
    'lp_mean5000','ld_mean5000');

histp_g = lme_bio00;
histp_gm2 = lme_mbio00;

histp_biome = biome_mbio50;

%Mean global biomass in g/m2
all_histp(:,1) = sf_mean5000+sp_mean5000+sd_mean5000+mf_mean5000+...
    mp_mean5000+md_mean5000+lp_mean5000+ld_mean5000;
all_histp(:,2) = sf_mean5000+mf_mean5000; %all f
all_histp(:,3) = sp_mean5000+mp_mean5000+lp_mean5000; %all p
all_histp(:,4) = sd_mean5000+md_mean5000+ld_mean5000; %all d
all_histp(:,5) = sf_mean5000+sp_mean5000+sd_mean5000; %all s
all_histp(:,6) = mf_mean5000+mp_mean5000+md_mean5000; %all m
all_histp(:,7) = lp_mean5000+ld_mean5000; %all l

%Mean total global biomass in g
tot_all_histp = all_histp .* repmat(grid(:,5),1,7);

clear lme_bio00 lme_mbio00 biome_mbio50 sf_mean5000 sp_mean5000 sd_mean5000
clear mf_mean5000 mp_mean5000 md_mean5000 b_mean5000 lp_mean5000 ld_mean5000

%% Hist fished
load([dpath 'LME_hist_fished_' cfile '.mat'],'lme_bio00','lme_mbio00');
load([dpath 'Biomes_hist_fished_' cfile '.mat'],'biome_mbio50');
load([dpath 'Means_hist_fished_' cfile '.mat'],'sf_mean5000','sp_mean5000',...
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

%% Global

Msum_pre = nansum(all_pre);
Msum_histp = nansum(all_histp);
Msum_histf = nansum(all_histf);

Tsum_pre = nansum(tot_all_pre);
Tsum_histp = nansum(tot_all_histp);
Tsum_histf = nansum(tot_all_histf);

%Indust vs. Preindust
Mdiff_histp_pre = (Msum_histp-Msum_pre) ./ Msum_pre;
Tdiff_histp_pre = (Tsum_histp-Tsum_pre) ./ Tsum_pre;
%Fishing vs. no fishing
Mdiff_histf_histp = (Msum_histf-Msum_histp) ./ Msum_histp;
Tdiff_histf_histp = (Tsum_histf-Tsum_histp) ./ Tsum_histp;

% By fish type
figure(1)
bar(Mdiff_histp_pre(2:4)*100)
colormap(cmap2)
ylim([-40 40])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
title('Global change in mean fish abundance (g/m^2) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'Global_fish_mean_biom_preindust_histp_type.png'])

figure(2)
bar(Mdiff_histf_histp(2:4)*100)
colormap(cmap2)
ylim([-40 40])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
title('Global change in mean fish abundance (g/m^2) in 1950-2000 with fishing')
print('-dpng',[ppath 'Global_fish_mean_biom_histp_histf_type.png'])

figure(3)
bar(Mdiff_histp_pre(5:7)*100)
colormap(cmap2)
ylim([-40 40])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('Global change in mean fish abundance (g/m^2) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'Global_fish_mean_biom_preindust_histp_size.png'])

figure(4)
bar(Mdiff_histf_histp(5:7)*100)
colormap(cmap2)
ylim([-40 40])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('Global change in mean fish abundance (g/m^2) in 1950-2000 with fishing')
print('-dpng',[ppath 'Global_fish_mean_biom_histp_histf_size.png'])

figure(5)
bar(Tdiff_histp_pre(2:4)*100)
colormap(cmap2)
ylim([-40 40])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
title('Global change in mean fish total abundance (g) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'Global_fish_tot_biom_preindust_histp_type.png'])

figure(6)
bar(Tdiff_histf_histp(2:4)*100)
colormap(cmap2)
ylim([-40 40])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
title('Global change in mean fish total abundance (g) in 1950-2000 with fishing')
print('-dpng',[ppath 'Global_fish_tot_biom_histp_histf_type.png'])

figure(7)
bar(Tdiff_histp_pre(5:7)*100)
colormap(cmap2)
ylim([-40 40])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('Global change in mean fish total abundance (g) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'Global_fish_tot_biom_preindust_histp_size.png'])

figure(8)
bar(Tdiff_histf_histp(5:7)*100)
colormap(cmap2)
ylim([-40 40])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('Global change in mean fish total abundance (g) in 1950-2000 with fishing')
print('-dpng',[ppath 'Global_fish_tot_biom_histp_histf_size.png'])

%% N Pac
[npid,ix,grid_np] = intersect(NPid,grid(:,1));

Msum_pre_np = nansum(all_pre(grid_np,:));
Msum_histp_np = nansum(all_histp(grid_np,:));
Msum_histf_np = nansum(all_histf(grid_np,:));

Tsum_pre_np = nansum(tot_all_pre(grid_np,:));
Tsum_histp_np = nansum(tot_all_histp(grid_np,:));
Tsum_histf_np = nansum(tot_all_histf(grid_np,:));

%Indust vs. Preindust
Mdiff_histp_pre_np = (Msum_histp_np - Msum_pre_np) ./ Msum_pre_np;
Tdiff_histp_pre_np = (Tsum_histp_np - Tsum_pre_np) ./ Tsum_pre_np;
%Fishing vs. no fishing
Mdiff_histf_histp_np = (Msum_histf_np - Msum_histp_np) ./ Msum_histp_np;
Tdiff_histf_histp_np = (Tsum_histf_np - Tsum_histp_np) ./ Tsum_histp_np;

% By fish type
figure(9)
bar(Mdiff_histp_pre_np(2:4)*100)
colormap(cmap2)
ylim([-50 50])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
title('N Pac change in mean fish abundance (g/m^2) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'NPac_fish_mean_biom_preindust_histp_type.png'])

figure(10)
bar(Mdiff_histf_histp_np(2:4)*100)
colormap(cmap2)
ylim([-50 50])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
title('N Pac change in mean fish abundance (g/m^2) in 1950-2000 with fishing')
print('-dpng',[ppath 'NPac_fish_mean_biom_histp_histf_type.png'])

figure(11)
bar(Mdiff_histp_pre_np(5:7)*100)
colormap(cmap2)
ylim([-50 50])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('N Pac change in mean fish abundance (g/m^2) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'NPac_fish_mean_biom_preindust_histp_size.png'])

figure(12)
bar(Mdiff_histf_histp_np(5:7)*100)
colormap(cmap2)
ylim([-50 50])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('N Pac change in mean fish abundance (g/m^2) in 1950-2000 with fishing')
print('-dpng',[ppath 'NPac_fish_mean_biom_histp_histf_size.png'])

figure(13)
bar(Tdiff_histp_pre_np(2:4)*100)
colormap(cmap2)
ylim([-50 50])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
title('N Pac change in mean fish total abundance (g) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'NPac_fish_tot_biom_preindust_histp_type.png'])

figure(14)
bar(Tdiff_histf_histp_np(2:4)*100)
colormap(cmap2)
ylim([-50 50])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
title('N Pac change in mean fish total abundance (g) in 1950-2000 with fishing')
print('-dpng',[ppath 'NPac_fish_tot_biom_histp_histf_type.png'])

figure(15)
bar(Tdiff_histp_pre_np(5:7)*100)
colormap(cmap2)
ylim([-50 50])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('N Pac change in mean fish total abundance (g) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'NPac_fish_tot_biom_preindust_histp_size.png'])

figure(16)
bar(Tdiff_histf_histp_np(5:7)*100)
colormap(cmap2)
ylim([-50 50])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
title('N Pac change in mean fish total abundance (g) in 1950-2000 with fishing')
print('-dpng',[ppath 'NPac_fish_tot_biom_histp_histf_size.png'])

%% Global biomes

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

histf_b(:,1) = nansum(histf_biome(:,1:8),2);
histf_b(:,2) = histf_biome(:,1) + histf_biome(:,4);
histf_b(:,3) = histf_biome(:,2) + histf_biome(:,5) + histf_biome(:,7);
histf_b(:,4) = histf_biome(:,3) + histf_biome(:,6) + histf_biome(:,8);
histf_b(:,5) = nansum(histf_biome(:,1:3),2);
histf_b(:,6) = nansum(histf_biome(:,4:6),2);
histf_b(:,7) = nansum(histf_biome(:,7:8),2);

%% NPac biomes grid cells
%Preindust
pb1 = find(biome_pre==1);
NPpb1 = intersect(NPid,pb1);
[nppb1,ia,grid_nppb1] = intersect(NPpb1,grid(:,1));
pb2 = find(biome_pre==2);
NPpb2 = intersect(NPid,pb2);
[nppb2,ia2,grid_nppb2] = intersect(NPpb2,grid(:,1));
pb3 = find(biome_pre==3);
NPpb3 = intersect(NPid,pb3);
[nppb3,ia3,grid_nppb3] = intersect(NPpb3,grid(:,1));
%Hist
hb1 = find(biome_hist==1);
NPhb1 = intersect(NPid,hb1);
[nphb1,ia4,grid_nphb1] = intersect(NPhb1,grid(:,1));
hb2 = find(biome_hist==2);
NPhb2 = intersect(NPid,hb2);
[nphb2,ia5,grid_nphb2] = intersect(NPhb2,grid(:,1));
hb3 = find(biome_hist==3);
NPhb3 = intersect(NPid,hb3);
[nphb3,ia6,grid_nphb3] = intersect(NPhb3,grid(:,1));

NPbiome_pre = zeros(360,200);
NPbiome_hist = NPbiome_pre;
%NPbiome_fore = NPbiome_pre;
NPbiome_pre(NPpb1) = ones(size(NPpb1));
NPbiome_hist(NPhb1) = ones(size(NPhb1));
%NPbiome_fore(ilc_fore) = ones(size(ilc_fore));
NPbiome_pre(NPpb2) = 2*ones(size(NPpb2));
NPbiome_hist(NPhb2) = 2*ones(size(NPhb2));
%NPbiome_fore(ieccs_fore) = 2*ones(size(ieccs_fore));
NPbiome_pre(NPpb3) = 3*ones(size(NPpb3));
NPbiome_hist(NPhb3) = 3*ones(size(NPhb3));
%NPbiome_fore(iecss_fore) = 3*ones(size(iecss_fore));

%% test correct on map
%plot info
geolon_t = double(geolon_t);
geolat_t = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

figure(50)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,NPbiome_pre)
colormap(cmap1)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gcf,'renderer','painters')
title('N Pac Biomes Pre-Industrial')
print('-dpng',[cpath 'NPac_biomes_preindust.png'])

figure(51)
axesm ('Robinson','MapLatLimit',[0 80],'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,NPbiome_hist)
colormap(cmap1)          
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gcf,'renderer','painters')
title('N Pac Biomes Historic')
print('-dpng',[cpath 'NPac_biomes_historic.png'])

figure(52)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,biome_pre)
colormap(cmap1)              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gcf,'renderer','painters')
title('Ocean Biomes Pre-Industrial')
axesmui
print('-dpng',[cpath 'biomes_preindust.png'])

figure(53)
axesm ('mollweid','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,biome_hist)
colormap(cmap1)
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
set(gcf,'renderer','painters')
title('Ocean Biomes Historic')
axesmui
print('-dpng',[cpath 'biomes_historic.png'])

%% calcs
%Only total g to see effect of change in biome size
sum_pre_np1 = nansum(tot_all_pre(grid_nppb1,:));
sum_pre_np2 = nansum(tot_all_pre(grid_nppb2,:));
sum_pre_np3 = nansum(tot_all_pre(grid_nppb3,:));

sum_histp_np1 = nansum(tot_all_histp(grid_nphb1,:));
sum_histp_np2 = nansum(tot_all_histp(grid_nphb2,:));
sum_histp_np3 = nansum(tot_all_histp(grid_nphb3,:));

sum_histf_np1 = nansum(tot_all_histf(grid_nphb1,:));
sum_histf_np2 = nansum(tot_all_histf(grid_nphb2,:));
sum_histf_np3 = nansum(tot_all_histf(grid_nphb3,:));

%Indust vs. Preindust
diff_histp_pre_np(1,:) = (sum_histp_np1 - sum_pre_np1) ./ sum_pre_np1;
diff_histp_pre_np(2,:) = (sum_histp_np2 - sum_pre_np2) ./ sum_pre_np2;
diff_histp_pre_np(3,:) = (sum_histp_np3 - sum_pre_np3) ./ sum_pre_np3;
%Fishing vs. no fishing
diff_histf_histp_np(1,:) = (sum_histf_np1 - sum_histp_np1) ./ sum_histp_np1;
diff_histf_histp_np(2,:) = (sum_histf_np2 - sum_histp_np2) ./ sum_histp_np2;
diff_histf_histp_np(3,:) = (sum_histf_np3 - sum_histp_np3) ./ sum_histp_np3;

%% By fish type
figure(17)
bar(diff_histp_pre_np(:,2:4)*100)
colormap(cmap2)
ylim([-60 60])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
legend('LC','ECCS','ECSS')
legend('location','northwest')
title('N Pac change in mean fish total abundance (g) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'NPacBiome_fish_tot_biom_preindust_histp_type.png'])

figure(18)
bar(diff_histf_histp_np(:,2:4)*100)
colormap(cmap2)
ylim([-60 60])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
xlabel('Fish type')
legend('LC','ECCS','ECSS')
legend('location','southeast')
title('N Pac change in mean fish total abundance (g) in 1950-2000 with fishing')
print('-dpng',[ppath 'NPacBiome_fish_tot_biom_histp_histf_type.png'])

figure(19)
bar(diff_histp_pre_np(:,5:7)*100)
colormap(cmap2)
ylim([-100 100])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
legend('LC','ECCS','ECSS')
legend('location','southeast')
title('N Pac change in mean fish total abundance (g) from 1800-1850 to 1950-2000')
print('-dpng',[ppath 'NPacBiome_fish_tot_biom_preindust_histp_size.png'])

figure(20)
bar(diff_histf_histp_np(:,5:7)*100)
colormap(cmap2)
ylim([-100 100])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
xlabel('Fish size')
legend('LC','ECCS','ECSS')
legend('location','northeast')
title('N Pac change in mean fish total abundance (g) in 1950-2000 with fishing')
print('-dpng',[ppath 'NPacBiome_fish_tot_biom_histp_histf_size.png'])

%% Put things together

Mdiff(1,:) = Mdiff_histp_pre;
Mdiff(2,:) = Mdiff_histf_histp;
Mdiff(3,:) = Mdiff_histp_pre_np;
Mdiff(4,:) = Mdiff_histf_histp_np;

Tdiff(1,:) = Tdiff_histp_pre;
Tdiff(2,:) = Tdiff_histf_histp;
Tdiff(3,:) = Tdiff_histp_pre_np;
Tdiff(4:6,:) = diff_histp_pre_np;
Tdiff(7,:) = Tdiff_histf_histp_np;
Tdiff(8:10,:) = diff_histf_histp_np;

%% Global pre-indust & indust
figure(21)
bar(Mdiff(1:2,1)*100)
colormap(cmap2)
xlim([0.5 2.5])
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
title('Global change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'Global_fish_mean_biom.png'])

figure(22)
bar(Mdiff(1:2,2:4)*100)
colormap(cmap2)
%ylim([-40 40])
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
legend('Forage','Pelagic','Demersal')
legend('location','southwest')
title('Global change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'Global_fish_mean_biom_type.png'])

figure(23)
bar(Mdiff(1:2,5:7)*100)
colormap(cmap2)
%ylim([-40 40])
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
legend('S','M','L')
legend('location','southwest')
title('Global change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'Global_fish_mean_biom_size.png'])

figure(24)
bar(Mdiff(1:2,:)*100)
colormap(cmap_ppt)
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
legend('All','Forage','Pelagic','Demersal','S','M','L')
legend('location','southwest')
title('Global change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'Global_fish_mean_biom_all.png'])
%%
figure(25)
bar([Mdiff(1:2,2:4)]'*100)
colormap(cmap2)
%ylim([-40 40])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
%xlabel('Driver')
legend('Industrial CO_2','Fishing')
legend('location','southwest')
title('Global change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'Global_fish_mean_biom_type_flip.png'])

figure(26)
bar([Mdiff(1:2,5:7)]'*100)
colormap(cmap2)
%ylim([-40 40])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
%xlabel('Driver')
legend('Industrial CO_2','Fishing')
legend('location','southwest')
title('Global change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'Global_fish_mean_biom_size_flip.png'])

figure(27)
bar([Mdiff(1:2,:)]'*100)
colormap(cmap2)
set(gca,'XTickLabel',{'All','Forage','Pelagic','Demersal','S','M','L'})
ylabel('Percent change')
%xlabel('Driver')
legend('Industrial CO_2','Fishing')
legend('location','southwest')
title('Global change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'Global_fish_mean_biom_all_flip.png'])

%% NPac pre-indust & indust
MdiffA(1,1) = Mdiff(1,1);
MdiffA(1,2) = Mdiff(3,1);
MdiffA(2,1) = Mdiff(2,1);
MdiffA(2,2) = Mdiff(4,1);

figure(28)
bar(MdiffA*100)
colormap(cmap2)
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
legend('Global','N Pacific')
legend('location','southwest')
title('Change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'GlobalvsNPac_fish_mean_biom.png'])

figure(29)
bar(Mdiff(3:4,1)*100)
colormap(cmap2)
xlim([0.5 2.5])
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
title('N Pac change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'NPac_fish_mean_biom.png'])

figure(30)
bar(Mdiff(3:4,2:4)*100)
colormap(cmap2)
%ylim([-40 40])
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
legend('Forage','Pelagic','Demersal')
legend('location','southwest')
title('N Pac change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'NPac_fish_mean_biom_type.png'])

figure(31)
bar(Mdiff(3:4,5:7)*100)
colormap(cmap2)
%ylim([-40 40])
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
legend('S','M','L')
legend('location','southwest')
title('N Pac change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'NPac_fish_mean_biom_size.png'])

figure(32)
bar(Mdiff(3:4,:)*100)
colormap(cmap_ppt)
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
legend('All','Forage','Pelagic','Demersal','S','M','L')
legend('location','southwest')
title('N Pac change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'NPac_fish_mean_biom_all.png'])
%%
figure(33)
bar([Mdiff(3:4,2:4)]'*100)
colormap(cmap2)
%ylim([-40 40])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
%xlabel('Driver')
legend('Industrial CO_2','Fishing')
legend('location','southwest')
title('N Pac change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'NPac_fish_mean_biom_type_flip.png'])

figure(34)
bar([Mdiff(3:4,5:7)]'*100)
colormap(cmap2)
%ylim([-40 40])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
%xlabel('Driver')
legend('Industrial CO_2','Fishing')
legend('location','southwest')
title('N Pac change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'NPac_fish_mean_biom_size_flip.png'])

figure(35)
bar([Mdiff(3:4,:)]'*100)
colormap(cmap2)
set(gca,'XTickLabel',{'All','Forage','Pelagic','Demersal','S','M','L'})
ylabel('Percent change')
%xlabel('Driver')
legend('Industrial CO_2','Fishing')
legend('location','southwest')
title('N Pac change in mean fish abundance (g/m^2)')
print('-dpng',[ppath 'NPac_fish_mean_biom_all_flip.png'])

%% N Pac biomes
TdiffA(1,1) = Tdiff(3,1);
TdiffA(1,2) = Tdiff(4,1);
TdiffA(1,3) = Tdiff(5,1);
TdiffA(1,4) = Tdiff(6,1);
TdiffA(2,1) = Tdiff(7,1);
TdiffA(2,2) = Tdiff(8,1);
TdiffA(2,3) = Tdiff(9,1);
TdiffA(2,4) = Tdiff(10,1);

figure(36)
bar(TdiffA*100)
colormap(cmap3)
set(gca,'XTickLabel',{'Industrial CO_2','Fishing'})
ylabel('Percent change')
%xlabel('Driver')
legend('All','LC','ECCS','ECSS')
legend('location','southwest')
title('N Pac change in mean fish total abundance (g)')
print('-dpng',[ppath 'NPacBiome_tot_mean_biom.png'])

figure(37)
subplot(1,2,1)
bar([Tdiff(3:6,2:4)]'*100)
colormap(cmap3)
ylim([-60 30])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
legend('All','LC','ECCS','ECSS')
legend('location','southwest')
text(4.5,37,'N Pac change in mean fish total abundance (g)','HorizontalAlignment','center')
title('Industrial CO_2')
subplot(1,2,2)
bar([Tdiff(7:10,2:4)]'*100)
colormap(cmap3)
ylim([-60 30])
set(gca,'XTickLabel',{'Forage','Pelagic','Demersal'})
ylabel('Percent change')
title('Fishing')
print('-dpng',[ppath 'NPacBiome_fish_tot_biom_type.png'])

figure(38)
subplot(1,2,1)
bar([Tdiff(3:6,5:7)]'*100)
colormap(cmap3)
ylim([-80 100])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
legend('All','LC','ECCS','ECSS')
legend('location','southwest')
text(4.5,112,'N Pac change in mean fish total abundance (g)','HorizontalAlignment','center')
title('Industrial CO_2')
subplot(1,2,2)
bar([Tdiff(7:10,5:7)]'*100)
colormap(cmap3)
ylim([-80 100])
set(gca,'XTickLabel',{'S','M','L'})
ylabel('Percent change')
title('Fishing')
print('-dpng',[ppath 'NPacBiome_fish_tot_biom_size.png'])



%% N Pac LMEs

%LME areal biomass in g
% pre(:,1) = nansum(pre_g(:,1:8),2);
% pre(:,2) = pre_g(:,1) + pre_g(:,4); %all f
% pre(:,3) = pre_g(:,2) + pre_g(:,5) + pre_g(:,7); %all p
% pre(:,4) = pre_g(:,3) + pre_g(:,6) + pre_g(:,8); %all d
% pre(:,5) = nansum(pre_g(:,1:3),2); %all s
% pre(:,6) = nansum(pre_g(:,4:6),2); %all m
% pre(:,7) = nansum(pre_g(:,7:8),2); %all l
% 
% histp(:,1) = nansum(histp_g(:,1:8),2);
% histp(:,2) = histp_g(:,1) + histp_g(:,4);
% histp(:,3) = histp_g(:,2) + histp_g(:,5) + histp_g(:,7);
% histp(:,4) = histp_g(:,3) + histp_g(:,6) + histp_g(:,8);
% histp(:,5) = nansum(histp_g(:,1:3),2);
% histp(:,6) = nansum(histp_g(:,4:6),2);
% histp(:,7) = nansum(histp_g(:,7:8),2);
% 
% %NPac LMEs = 47-57; 1-4; 10-11; 65
% 
% np = [1:4 10:11 47:57 65];
% 
% NPpre = pre(np,:);
% NPhistp = histp(np,:);
% 
% NPdiff_histp_pre = (NPhistp - NPpre) ./ NPpre;
% NPdiff_histp_pre_tot = (sum(NPhistp) - sum(NPpre)) ./ sum(NPpre);
% 
% % Figures
% test(:,1)=histp_g(:,1);
% test(:,2)=pre_g(:,2);
% figure(100)
% bar(test)
% legend('histp','pre')
% 
% % By fish type
% figure(1)
% bar(NPdiff_histp_pre_tot(2:4)*100)
% set(gca,'XTickLabel',{'S','M','L'})
% ylabel('Percent change')
% xlabel('Fish size')
% title('Change in total fish abundance (g) from 1800-1850 to 1950-2000')
% %print('-dpng',[ppath 'NPac_Fish_biom_preindust_histprist.png'])
% 
% % By LME
% figure(2)
% bar(diff_histp_pre(:,1)*100)
% set(gca,'XTick',1:length(np),'XTickLabel',np)
% ylabel('Percent change')
% xlabel('LME')
% title('Change in total fish abundance (g) from 1800-1850 to 1950-2000')
% %print('-dpng',[ppath 'NPac_LME_biom_preindust_histprist.png'])

