% Visualize output of POEM
% Spinup at 100 locations
% 50 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
% dpath = '/Volumes/GFDL/NC/MFbetterMP4_fcrit05/';
% ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit05/';
% dpath = '/Volumes/GFDL/NC/MFbetterMP4_fcrit10/';
% ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/';
% dpath = '/Volumes/GFDL/NC/MFeqMP4_fcrit10_Tmort/';
% ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/';
% dpath = '/Volumes/GFDL/NC/fcrit10_FdiffA2_Tmort/';
% ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2_Tmort/';
% dpath = '/Volumes/GFDL/NC/Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% dpath = '/Volumes/GFDL/NC/Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50/';
% ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

%load([dpath 'Data_spinup_pristine.mat']);
%load([dpath 'Data_spinup_pristine_fcrit10.mat']);
%load([dpath 'Data_spinup_pristine_fcrit10_Tmort_30yr.mat']);
%load([dpath 'Data_spinup_pristine_fcrit10_FdiffA2_Tmort.mat']);
% load([dpath 'Data_spinup_pristine_Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50.mat']);
load([dpath 'Data_spinup_pristine_' cfile '.mat']);


%%
[loc,days]=size(SP.bio);
x=1:days;

lyr=x((end-365+1):end);

sp_mean=mean(SP.bio(:,lyr),2);
sf_mean=mean(SF.bio(:,lyr),2);
sd_mean=mean(SD.bio(:,lyr),2);
mp_mean=mean(MP.bio(:,lyr),2);
mf_mean=mean(MF.bio(:,lyr),2);
md_mean=mean(MD.bio(:,lyr),2);
lp_mean=mean(LP.bio(:,lyr),2);
ld_mean=mean(LD.bio(:,lyr),2);
b_mean=mean(BENT.bio(:,lyr),2);

% Plots in space
grid = csvread([cpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

Zsp=griddata(grid(:,2),grid(:,3),sp_mean,X,Y);
Zsf=griddata(grid(:,2),grid(:,3),sf_mean,X,Y);
Zsd=griddata(grid(:,2),grid(:,3),sd_mean,X,Y);
Zmp=griddata(grid(:,2),grid(:,3),mp_mean,X,Y);
Zmf=griddata(grid(:,2),grid(:,3),mf_mean,X,Y);
Zmd=griddata(grid(:,2),grid(:,3),md_mean,X,Y);
Zlp=griddata(grid(:,2),grid(:,3),lp_mean,X,Y);
Zld=griddata(grid(:,2),grid(:,3),ld_mean,X,Y);
Zb=griddata(grid(:,2),grid(:,3),b_mean,X,Y);

%% bent
figure(50)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zb))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass of Benthic Invertebrates (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_BENT_NF.png'])

%% bent mg C
mgZb = (Zb/9)*1e3;
figure(51)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(mgZb))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean benthic biomass (mg C m^-^2)')
colormap('jet')
colorbar('h')
caxis([-1 3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_BENT_mgC_NF.png'])


%% Diff maps of all fish
All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
AllF = Zsf+Zmf;
AllP = Zsp+Zmp+Zlp;
AllD = Zsd+Zmd+Zld;
AllS = Zsp+Zsf+Zsd;
AllM = Zmp+Zmf+Zmd;
AllL = Zlp+Zld;
FracPD = AllP ./ (AllP+AllD);
FracPF = AllP ./ (AllP+AllF);
FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);
FracPDs = Zsp ./ (Zsp+Zsd);
FracPDm = Zmp ./ (Zmp+Zmd);
FracPDl = Zlp ./ (Zlp+Zld);
FracPFs = Zsp ./ (Zsp+Zsf);
FracPFm = Zmp ./ (Zmp+Zmf);
FracPFvDs = (Zsp+Zsf) ./ (Zsp+Zsf+Zsd);
FracPFvDm = (Zmp+Zmf) ./ (Zmp+Zmf+Zmd);

%% ALL
figure(21)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(All))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass of All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-1 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_All_NF1.png'])
%%
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(All))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass of All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_All_NF2.png'])

%% all F
figure(22)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(AllF))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass of Forage Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllF_NF.png'])

% all D
figure(23)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(AllD))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass of Demersal Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllD_NF.png'])

% All P
figure(24)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(AllP))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass of Large Pelagics (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllP_NF.png'])

