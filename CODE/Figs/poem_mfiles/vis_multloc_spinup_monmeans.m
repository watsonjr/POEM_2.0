% Visualize output of POEM
% Spinup at 100 locations
% 50 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Jul_og_sizes/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Julia_OG_sizes/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC050_RE0100';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

%load([dpath 'Data_spinup_pristine_' cfile '.mat']);
load([dpath 'Means_spinup_' cfile '_5yr.mat']);

%%
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
print('-dpng',[ppath 'Spinup_pisc_time.png'])

%% Planktivore
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

print('-dpng',[ppath 'Spinup_plan_time.png'])

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
print('-dpng',[ppath 'Spinup_detr_time.png'])

%% All size classes of all

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
print('-dpng',[ppath 'Spinup_all_sizes_nanmean.png'])

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
% figure(55)
% m_proj('miller','lat',82);
% surf(geolon_t,geolat_t,ocean); view(2); hold on;
% shading flat
% % m_coast('patch',[.5 .5 .5],'edgecolor','none');
% % m_grid;
% title('Water cells')
% colormap('jet')
% colorbar('h')
% caxis([1 2])
% stamp(cfile)
% print('-dpng',[ppath 'Ocean_cells.png'])

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



%% Plots in space
% grid = csvread([cpath 'grid_csv.csv']);
% %fix lon shift
% id=find(grid(:,2)<-180);
% grid(id,2)=grid(id,2)+360;
% 
% x=-180:180;
% y=-90:90;
% [X,Y]=meshgrid(x,y);
% 
% Zsp=griddata(grid(:,2),grid(:,3),sp_mean,X,Y);
% Zsf=griddata(grid(:,2),grid(:,3),sf_mean,X,Y);
% Zsd=griddata(grid(:,2),grid(:,3),sd_mean,X,Y);
% Zmp=griddata(grid(:,2),grid(:,3),mp_mean,X,Y);
% Zmf=griddata(grid(:,2),grid(:,3),mf_mean,X,Y);
% Zmd=griddata(grid(:,2),grid(:,3),md_mean,X,Y);
% Zlp=griddata(grid(:,2),grid(:,3),lp_mean,X,Y);
% Zld=griddata(grid(:,2),grid(:,3),ld_mean,X,Y);
% Zb=griddata(grid(:,2),grid(:,3),b_mean,X,Y);
% 
% %% bent
% figure(50)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(Zb))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean benthic biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2.5 0.5])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_BENT.png'])
% 
% %
% mgZb = (Zb/9)*1e3;
% figure(51)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(mgZb))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean benthic biomass (mg C m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-0.8 2.3])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_BENT_mgC.png'])
% 
% % sp
% figure(6)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(Zsp))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean Larval P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_SP.png'])
% 
% % sf
% figure(7)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(Zsf))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean Larval F biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_SF.png'])
% 
% % sd
% figure(8)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(Zsd))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean Larval D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_SD.png'])
% 
% % mp
% figure(9)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(Zmp))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean Juvenile P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_MP.png'])
% 
% % mf
% figure(10)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(Zmf))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% %m_plot(-111,5,'o','w','MarkerSize',10);
% m_text(-111,5,'EEP','Color','red','HorizontalAlignment','center');
% title('log10 mean Adult F biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_MF.png'])
% 
% % md
% figure(11)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(Zmd))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean Juvenile D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_MD.png'])
% 
% % lp
% figure(12)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(Zlp))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
% title('log10 mean Adult P biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_LP.png'])
% 
% % ld
% figure(13)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(Zld))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
% title('log10 mean Adult D biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_LD.png'])
% 
% %% Diff maps of all fish
% All = Zsp+Zsf+Zsd+Zmp+Zmf+Zmd+Zlp+Zld;
% AllF = Zsf+Zmf;
% AllP = Zsp+Zmp+Zlp;
% AllD = Zsd+Zmd+Zld;
% AllS = Zsp+Zsf+Zsd;
% AllM = Zmp+Zmf+Zmd;
% AllL = Zlp+Zld;
% FracPD = AllP ./ (AllP+AllD);
% FracPF = AllP ./ (AllP+AllF);
% FracPFvD = (AllP+AllF) ./ (AllP+AllF+AllD);
% FracPDs = Zsp ./ (Zsp+Zsd);
% FracPDm = Zmp ./ (Zmp+Zmd);
% FracPDl = Zlp ./ (Zlp+Zld);
% FracPFs = Zsp ./ (Zsp+Zsf);
% FracPFm = Zmp ./ (Zmp+Zmf);
% FracPFvDs = (Zsp+Zsf) ./ (Zsp+Zsf+Zsd);
% FracPFvDm = (Zmp+Zmf) ./ (Zmp+Zmf+Zmd);
% 
% %% ALL
% figure(21)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(All))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean biomass All Fishes (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-1 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_All.png'])
% 
% % all F
% figure(22)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(AllF))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean biomass All F (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_AllF.png'])
% 
% % all D
% figure(23)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(AllD))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean biomass All D (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_AllD.png'])
% 
% % All P
% figure(24)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(log10(AllP))); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('log10 mean biomass All P (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([-2 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_AllP.png'])
% 
% % FracPD
% figure(25)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(FracPD)); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('P:D mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_FracPD.png'])
% 
% % FracPF
% figure(26)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(FracPF)); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('P:F mean biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_FracPF.png'])
% 
% %% FracPFvD
% figure(27)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(FracPFvD)); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('(P+F):D mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_FracPFvD.png'])
% 
% % FracPDs
% figure(28)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(FracPDs)); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('SP:SD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_FracPDs.png'])
% 
% % FracPFs
% figure(29)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(FracPFs)); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('SP:SF mean biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_FracPFs.png'])
% 
% % FracPFvDs
% figure(30)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(FracPFvDs)); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('(SP+SF):SD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_FracPFvDs.png'])
% 
% % FracPDm
% figure(31)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(FracPDm)); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('MP:MD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_FracPDm.png'])
% 
% % FracPFm
% figure(32)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(FracPFm)); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('MP:MF mean biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_FracPFm.png'])
% 
% % FracPFvDm
% figure(33)
% m_proj('miller','lat',82);
% m_pcolor(X,Y,real(FracPFvDm)); hold on;
% shading flat
% m_coast('patch',[.5 .5 .5],'edgecolor','none');
% m_grid;
% title('(MP+MF):MD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_global_FracPFvDm.png'])
% 
% %% DEBUG ---------------------------------------------------
% % Find NaNs
% mos=291:294;
% tend=length(mos);
% bio=NaN*ones(9,tend);
% % con=NaN*ones(8,tend);
% % die=NaN*ones(8,tend);
% % gamma=NaN*ones(8,tend);
% % rec=NaN*ones(8,tend);
% % nu=NaN*ones(8,tend);
% 
% for k=1:tend
%     n=mos(k);
%     %bio
%     id = isnan(SP.bio(:,n));
%     bio(1,k) = sum(id);
%     id = isnan(SF.bio(:,n));
%     bio(2,k) = sum(id);
%     id = isnan(SD.bio(:,n));
%     bio(3,k) = sum(id);
%     id = isnan(MP.bio(:,n));
%     bio(4,k) = sum(id);
%     id = isnan(MF.bio(:,n));
%     bio(5,k) = sum(id);
%     id = isnan(MD.bio(:,n));
%     bio(6,k) = sum(id);
%     id = isnan(LP.bio(:,n));
%     bio(7,k) = sum(id);
%     id = isnan(LD.bio(:,n));
%     bio(8,k) = sum(id);
%     id = isnan(BENT.bio(:,n));
%     bio(9,k) = sum(id);
% %     %con
% %     id = isnan(SP.con(:,n));
% %     con(1,k) = sum(id);
% %     id = isnan(SF.con(:,n));
% %     con(2,k) = sum(id);
% %     id = isnan(SD.con(:,n));
% %     con(3,k) = sum(id);
% %     id = isnan(MP.con(:,n));
% %     con(4,k) = sum(id);
% %     id = isnan(MF.con(:,n));
% %     con(5,k) = sum(id);
% %     id = isnan(MD.con(:,n));
% %     con(6,k) = sum(id);
% %     id = isnan(LP.con(:,n));
% %     con(7,k) = sum(id);
% %     id = isnan(LD.con(:,n));
% %     con(8,k) = sum(id);
% %     %die
% %     id = isnan(SP.die(:,n));
% %     die(1,k) = sum(id);
% %     id = isnan(SF.die(:,n));
% %     die(2,k) = sum(id);
% %     id = isnan(SD.die(:,n));
% %     die(3,k) = sum(id);
% %     id = isnan(MP.die(:,n));
% %     die(4,k) = sum(id);
% %     id = isnan(MF.die(:,n));
% %     die(5,k) = sum(id);
% %     id = isnan(MD.die(:,n));
% %     die(6,k) = sum(id);
% %     id = isnan(LP.die(:,n));
% %     die(7,k) = sum(id);
% %     id = isnan(LD.die(:,n));
% %     die(8,k) = sum(id);
% %     %gamma/rep
% %     id = isnan(SP.gamma(:,n));
% %     gamma(1,k) = sum(id);
% %     id = isnan(SF.gamma(:,n));
% %     gamma(2,k) = sum(id);
% %     id = isnan(SD.gamma(:,n));
% %     gamma(3,k) = sum(id);
% %     id = isnan(MP.gamma(:,n));
% %     gamma(4,k) = sum(id);
% %     id = isnan(MF.rep(:,n));
% %     gamma(5,k) = sum(id);
% %     id = isnan(MD.gamma(:,n));
% %     gamma(6,k) = sum(id);
% %     id = isnan(LP.rep(:,n));
% %     gamma(7,k) = sum(id);
% %     id = isnan(LD.rep(:,n));
% %     gamma(8,k) = sum(id);
% %     %rec
% %     id = isnan(SP.rec(:,n));
% %     rec(1,k) = sum(id);
% %     id = isnan(SF.rec(:,n));
% %     rec(2,k) = sum(id);
% %     id = isnan(SD.rec(:,n));
% %     rec(3,k) = sum(id);
% %     id = isnan(MP.rec(:,n));
% %     rec(4,k) = sum(id);
% %     id = isnan(MF.rec(:,n));
% %     rec(5,k) = sum(id);
% %     id = isnan(MD.rec(:,n));
% %     rec(6,k) = sum(id);
% %     id = isnan(LP.rec(:,n));
% %     rec(7,k) = sum(id);
% %     id = isnan(LD.rec(:,n));
% %     rec(8,k) = sum(id);
% %     %nu
% %     id = isnan(SP.nu(:,n));
% %     nu(1,k) = sum(id);
% %     id = isnan(SF.nu(:,n));
% %     nu(2,k) = sum(id);
% %     id = isnan(SD.nu(:,n));
% %     nu(3,k) = sum(id);
% %     id = isnan(MP.nu(:,n));
% %     nu(4,k) = sum(id);
% %     id = isnan(MF.nu(:,n));
% %     nu(5,k) = sum(id);
% %     id = isnan(MD.nu(:,n));
% %     nu(6,k) = sum(id);
% %     id = isnan(LP.nu(:,n));
% %     nu(7,k) = sum(id);
% %     id = isnan(LD.nu(:,n));
% %     nu(8,k) = sum(id);
%     
% end
% %%
% figure
% bar(bio')
% set(gca,'XTickLabel',mos)
% legend('SP','SF','SD','MP','MF','MD','LP','LD','B')
% legend('location','northwest')
% print('-dpng',[ppath 'nan_biotest.png'])
% 
% % figure
% % bar(con')
% % set(gca,'XTickLabel',mos)
% % legend('SP','SF','SD','MP','MF','MD','LP','LD')
% % legend('location','northwest')
% % print('-dpng',[ppath 'nan_contest.png'])
% % 
% % figure
% % bar(die')
% % set(gca,'XTickLabel',mos)
% % legend('SP','SF','SD','MP','MF','MD','LP','LD')
% % legend('location','northwest')
% % print('-dpng',[ppath 'nan_dietest.png'])
% % 
% % figure
% % bar(gamma')
% % set(gca,'XTickLabel',mos)
% % legend('SP','SF','SD','MP','MF','MD','LP','LD')
% % legend('location','northwest')
% % print('-dpng',[ppath 'nan_gammatest.png'])
% % 
% % figure
% % bar(rec')
% % legend('SP','SF','SD','MP','MF','MD','LP','LD')
% % legend('location','northwest')
% % print('-dpng',[ppath 'nan_rectest.png'])
% % 
% % figure
% % bar(nu')
% % set(gca,'XTickLabel',mos)
% % legend('SP','SF','SD','MP','MF','MD','LP','LD')
% % legend('location','northwest')
% % print('-dpng',[ppath 'nan_nutest.png'])
% 
% %%
% n=293;
% id=find(isnan(SP.bio(:,n)));    
% grid(id,2:3);                   
% % id2=find(isnan(SP.die(:,n)));
% % grid(id2,2:3);
% % id3=find(isnan(SP.con(:,n)));   
% % grid(id3,2:3);
% % id4=find(isnan(MP.con(:,n)));
% % grid(id4,2:3);
% % id5=find(isnan(SP.rec(:,n)))
% % grid(id5,2:3)
% 
