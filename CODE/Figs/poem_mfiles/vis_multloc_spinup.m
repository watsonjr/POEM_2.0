% Visualize output of POEM
% Spinup at 100 locations
% 30 years
% Saved as mat files

%clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
% dpath = '/Volumes/GFDL/NC/fcrit05/';
% ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4/';
% dpath = '/Volumes/GFDL/NC/fcrit10/';
% ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/';
dpath = '/Volumes/GFDL/NC/Tmort/';
ppath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/';

cfile = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort';

%load([dpath 'Data_spinup_pristine.mat']);
%load([dpath 'Data_spinup_pristine_fcrit10.mat']);
%load([dpath 'Data_spinup_pristine_fcrit10_Tmort.mat']);

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

%% Plots in space
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
m_pcolor(X,Y,real(log10(Zsp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SP.png'])

%% sf
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Zsf)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SF.png'])

%% sd
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Zsd)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SD.png'])

%% mp
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zmp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Juvenile P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MP.png'])

%% mf
figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Zmf)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MF.png'])

%% md
figure(6)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Zmd)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Juvenile D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MD.png'])

%% lp
figure(7)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Zlp)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LP.png'])

%% ld
figure(8)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Zld)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LD.png'])

%% Production
% Production = gamma + die + rep + egg
SP_prod=nanmean(SP.gamma(:,lyr)+SP.die(:,lyr)+SP.rep(:,lyr)+SP.egg(:,lyr),2);
SF_prod=nanmean(SF.gamma(:,lyr)+SF.die(:,lyr)+SF.rep(:,lyr)+SF.egg(:,lyr),2);
SD_prod=nanmean(SD.gamma(:,lyr)+SD.die(:,lyr)+SD.rep(:,lyr)+SD.egg(:,lyr),2);
MP_prod=nanmean(MP.gamma(:,lyr)+MP.die(:,lyr)+MP.rep(:,lyr)+MP.egg(:,lyr),2);
MF_prod=nanmean(MF.gamma(:,lyr)+MF.die(:,lyr)+MF.rep(:,lyr)+MF.egg(:,lyr),2);
MD_prod=nanmean(MD.gamma(:,lyr)+MD.die(:,lyr)+MD.rep(:,lyr)+MD.egg(:,lyr),2);
LP_prod=nanmean(LP.gamma(:,lyr)+LP.die(:,lyr)+LP.rep(:,lyr)+LP.egg(:,lyr),2);
LD_prod=nanmean(LD.gamma(:,lyr)+LD.die(:,lyr)+LD.rep(:,lyr)+LD.egg(:,lyr),2);

%%
Psp=griddata(grid(:,2),grid(:,3),SP_prod,X,Y);
Psf=griddata(grid(:,2),grid(:,3),SF_prod,X,Y);
Psd=griddata(grid(:,2),grid(:,3),SD_prod,X,Y);
Pmp=griddata(grid(:,2),grid(:,3),MP_prod,X,Y);
Pmf=griddata(grid(:,2),grid(:,3),MF_prod,X,Y);
Pmd=griddata(grid(:,2),grid(:,3),MD_prod,X,Y);
Plp=griddata(grid(:,2),grid(:,3),LP_prod,X,Y);
Pld=griddata(grid(:,2),grid(:,3),LD_prod,X,Y);

%% sp
figure(11)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Psp)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval P production (g g^-^1 m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 -2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_prod_SP.png'])

% sf
figure(12)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Psf)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval F production (g g^-^1 m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 -2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_prod_SF.png'])

% sd
figure(13)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Psd)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval D production (g g^-^1 m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 -2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_prod_SD.png'])

% mp
figure(14)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Pmp)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Juvenile P production (g g^-^1 m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 -2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_prod_MP.png'])

% mf
figure(15)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Pmf)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Adult F production (g g^-^1 m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 -2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_prod_MF.png'])

% md
figure(16)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Pmd)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Juvenile D production (g g^-^1 m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 -2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_prod_MD.png'])

%% lp
figure(17)
stamp(cfile)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Plp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Adult P production (g g^-^1 m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 -2])
print('-dpng',[ppath 'Spinup_global_prod_LP.png'])

% ld
figure(18)
stamp(cfile)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Pld)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Adult D production (g g^-^1 m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 -2])
print('-dpng',[ppath 'Spinup_global_prod_LD.png'])

