% Visualize output of POEM
% Spinup at 100 locations
% 50 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/';
% pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_OG_sizes/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC100_RE1000';

fpath=['/Volumes/GFDL/NC/Matlab_big_size/' cfile '/'];
% fpath=['/Volumes/GFDL/NC/Matlab_og_size/' cfile '/'];
ppath = [pp cfile '/'];

load([fpath 'Means_spinup_' cfile '.mat']);

%% colors
cm9=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...  %b
    0 0 0];...      %black
    
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm9);

%% Plots in time
y = time;

% Piscivore
figure(1)
subplot(4,1,1)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup Pelagic Piscivores')
ylabel('log10 Biomass (g m^-^2)')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Juveniles')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,4)
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Spinup_P_time.png'])

% Planktivore
sf_tmean=sf_tmean(1:length(y));
figure(2)
subplot(3,1,1)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup Forage Fishes')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
legend('Immature','Adults')
legend('location','southeast')
stamp(cfile)

subplot(3,1,2)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Immature')
ylabel('log10 Biomass (g m^-^2)')

subplot(3,1,3)
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Spinup_F_time.png'])

% Detritivore
figure(3)
subplot(4,1,1)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup Demersal Piscivores')
ylabel('log10 Biomass (g m^-^2)')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Juveniles')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,4)
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Spinup_D_time.png'])

% All size classes of all

figure(5)
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
legend('SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title('Spinup')
stamp(cfile)
print('-dpng',[ppath 'Spinup_all_sizes.png'])


%% Plots in space
grid = csvread([cpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');

[ni,nj]=size(geolon_t);

Zsf=NaN*ones(ni,nj);
Zsp=NaN*ones(ni,nj);
Zsd=NaN*ones(ni,nj);
Zmf=NaN*ones(ni,nj);
Zmp=NaN*ones(ni,nj);
Zmd=NaN*ones(ni,nj);
Zlp=NaN*ones(ni,nj);
Zld=NaN*ones(ni,nj);
Zb=NaN*ones(ni,nj);

ocean=NaN*ones(ni,nj);
ocean(grid(:,1))=ones(size(sf_mean));

Zsf(grid(:,1))=sf_mean;
Zsp(grid(:,1))=sp_mean;
Zsd(grid(:,1))=sd_mean;
Zmf(grid(:,1))=mf_mean;
Zmp(grid(:,1))=mp_mean;
Zmd(grid(:,1))=md_mean;
Zlp(grid(:,1))=lp_mean;
Zld(grid(:,1))=ld_mean;
Zb(grid(:,1))=b_mean;

%% ocean cells
figure(55)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,ocean); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('Water cells')
colormap('jet')
colorbar('h')
caxis([1 2])
stamp(cfile)
print('-dpng',[ppath 'Ocean_cells.png'])

% bent
figure(50)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zb)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean benthic biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_BENT.png'])

%
mgZb = (Zb/9)*1e3;
figure(51)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(mgZb)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean benthic biomass (mg C m^-^2)')
colormap('jet')
colorbar('h')
caxis([-0.8 2.3])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_BENT_mgC.png'])

% sp
figure(11)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zsp)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Larval P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SP.png'])

% sf
figure(12)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zsf)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Larval F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SF.png'])

% sd
figure(13)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zsd)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Larval D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_SD.png'])

% mp
figure(14)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zmp)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Juvenile P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MP.png'])

% mf
figure(15)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zmf)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MF.png'])

% md
figure(16)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zmd)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Juvenile D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MD.png'])

% lp
figure(17)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zlp)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LP.png'])

% ld
figure(18)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(Zld)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LD.png'])

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
surf(geolon_t,geolat_t,log10(All)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean biomass All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-1 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_All.png'])

% all F
figure(22)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(AllF)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean biomass All F (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllF.png'])

% all D
figure(23)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(AllD)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean biomass All D (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllD.png'])

% All P
figure(24)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,log10(AllP)); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('log10 mean biomass All P (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllP.png'])

% FracPD
figure(25)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPD); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('P:D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPD.png'])

% FracPF
figure(26)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPF); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('P:F mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPF.png'])

%% FracPFvD
figure(27)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFvD); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('(P+F):D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFvD.png'])

% FracPDs
figure(28)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPDs); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('SP:SD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPDs.png'])

% FracPFs
figure(29)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFs); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('SP:SF mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFs.png'])

% FracPFvDs
figure(30)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFvDs); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('(SP+SF):SD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFvDs.png'])

% FracPDm
figure(31)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPDm); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('MP:MD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPDm.png'])

% FracPFm
figure(32)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFm); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('MP:MF mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFm.png'])

% FracPFvDm
figure(33)
m_proj('miller','lat',82);
surf(geolon_t,geolat_t,FracPFvDm); view(2); hold on;
shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
title('(MP+MF):MD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFvDm.png'])

