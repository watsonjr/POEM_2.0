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
%load([dpath 'Means_spinup_' cfile '.mat']);


%% Feeding level
[loc,days]=size(SP.clev);
x=1:days;

lyr=x((end-365+1):end);

sp_mean=mean(SP.clev(:,lyr),2);
sf_mean=mean(SF.clev(:,lyr),2);
sd_mean=mean(SD.clev(:,lyr),2);
mp_mean=mean(MP.clev(:,lyr),2);
mf_mean=mean(MF.clev(:,lyr),2);
md_mean=mean(MD.clev(:,lyr),2);
lp_mean=mean(LP.clev(:,lyr),2);
ld_mean=mean(LD.clev(:,lyr),2);

%Plots in space
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

%% sp
figure(1)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zsp); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Larval P feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_clev_SP.png'])

%% sf
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zsf); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Larval F feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_clev_SF.png'])

% sd
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zsd); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Larval D feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_clev_SD.png'])

% mp
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zmp); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Juvenile P feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_clev_MP.png'])

% mf
figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zmf); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','red','HorizontalAlignment','center');
title('mean Adult F feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_clev_MF.png'])

% md
figure(6)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zmd); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Juvenile D feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_clev_MD.png'])

% lp
figure(7)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zlp); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('mean Adult P feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_clev_LP.png'])

% ld
figure(8)
m_proj('miller','lat',82);
m_pcolor(X,Y,Zld); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('mean Adult D feeding level')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_clev_LD.png'])

%% Consumption
L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

% convert g/g/d to g/d and sum over year to give g/yr
sp_con=sum(SP.con(:,lyr)*M_s,2);
sf_con=sum(SF.con(:,lyr)*M_s,2);
sd_con=sum(SD.con(:,lyr)*M_s,2);
mp_con=sum(MP.con(:,lyr)*M_m,2);
mf_con=sum(MF.con(:,lyr)*M_m,2);
md_con=sum(MD.con(:,lyr)*M_m,2);
lp_con=sum(LP.con(:,lyr)*M_l,2);
ld_con=sum(LD.con(:,lyr)*M_l,2);

Csp=griddata(grid(:,2),grid(:,3),sp_con,X,Y);
Csf=griddata(grid(:,2),grid(:,3),sf_con,X,Y);
Csd=griddata(grid(:,2),grid(:,3),sd_con,X,Y);
Cmp=griddata(grid(:,2),grid(:,3),mp_con,X,Y);
Cmf=griddata(grid(:,2),grid(:,3),mf_con,X,Y);
Cmd=griddata(grid(:,2),grid(:,3),md_con,X,Y);
Clp=griddata(grid(:,2),grid(:,3),lp_con,X,Y);
Cld=griddata(grid(:,2),grid(:,3),ld_con,X,Y);

%% sp
figure(1)
m_proj('miller','lat',82);
m_pcolor(X,Y,Csp); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Larval P consumption rate (g yr^-^1)')
colormap('jet')
colorbar('h')
%caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_con_SP.png'])

%% sf
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,Csf); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Larval F consumption rate (g yr^-^1)')
colormap('jet')
colorbar('h')
%caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_con_SF.png'])

% sd
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,Csd); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Larval D consumption rate (g yr^-^1)')
colormap('jet')
colorbar('h')
%caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_con_SD.png'])

% mp
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,Cmp); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Juvenile P consumption rate (g yr^-^1)')
colormap('jet')
colorbar('h')
%caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_con_MP.png'])

% mf
figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,Cmf); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','red','HorizontalAlignment','center');
title('mean Adult F consumption rate (g yr^-^1)')
colormap('jet')
colorbar('h')
%caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_con_MF.png'])

% md
figure(6)
m_proj('miller','lat',82);
m_pcolor(X,Y,Cmd); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('mean Juvenile D consumption rate (g yr^-^1)')
colormap('jet')
colorbar('h')
%caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_con_MD.png'])

% lp
figure(7)
m_proj('miller','lat',82);
m_pcolor(X,Y,Clp); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('mean Adult P consumption rate (g yr^-^1)')
colormap('jet')
colorbar('h')
%caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_con_LP.png'])

% ld
figure(8)
m_proj('miller','lat',82);
m_pcolor(X,Y,Cld); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('mean Adult D consumption rate (g yr^-^1)')
colormap('jet')
colorbar('h')
%caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_con_LD.png'])

