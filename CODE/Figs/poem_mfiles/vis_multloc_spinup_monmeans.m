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

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC275_RE0500';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

%load([dpath 'Data_spinup_pristine_' cfile '.mat']);
load([dpath 'Means_spinup_' cfile '.mat']);

%% Plots in time
y = SP.time;

%% Piscivore
figure(1)
subplot(4,1,1)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([1 12])
title('Spinup Pelagic Piscivores')
ylabel('log10 Biomass (g m^-^2)')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
xlim([1 12])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
xlim([1 12])
title('Juveniles')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,4)
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([1 12])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
%print('-dpng',[ppath 'Spinup_pisc_time.png'])

% Planktivore
figure(2)
subplot(3,1,1)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([1 12])
title('Spinup Forage Fishes')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
legend('Immature','Adults')
legend('location','southeast')
stamp(cfile)

subplot(3,1,2)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
xlim([1 12])
title('Immature')
ylabel('log10 Biomass (g m^-^2)')

subplot(3,1,3)
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([1 12])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')

%print('-dpng',[ppath 'Spinup_plan_time.png'])

% Detritivore
figure(3)
subplot(4,1,1)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([1 12])
title('Spinup Demersal Piscivores')
ylabel('log10 Biomass (g m^-^2)')
legend('Larvae','Juveniles','Adults')
legend('location','southeast')
stamp(cfile)

subplot(4,1,2)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
xlim([1 12])
title('Larvae')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,3)
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
xlim([1 12])
title('Juveniles')
ylabel('log10 Biomass (g m^-^2)')

subplot(4,1,4)
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([1 12])
title('Adults')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
%print('-dpng',[ppath 'Spinup_detr_time.png'])

% All size classes of all

figure(5)
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','eastoutside')
xlim([1 12])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title('Spinup')
stamp(cfile)
%print('-dpng',[ppath 'Spinup_all_sizes.png'])


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
Zb=griddata(grid(:,2),grid(:,3),b_mean,X,Y);

%% bent
figure(50)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zb))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean benthic biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_BENT.png'])

%
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
caxis([-0.8 2.3])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_BENT_mgC.png'])

% sp
figure(1)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zsp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_SP.png'])

% sf
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zsf))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_SF.png'])

% sd
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zsd))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_SD.png'])

% mp
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zmp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Juvenile P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_MP.png'])

% mf
figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zmf))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_plot(-111,5,'o','w','MarkerSize',10);
m_text(-111,5,'EEP','Color','red','HorizontalAlignment','center');
title('log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_MF.png'])

% md
figure(6)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zmd))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Juvenile D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_MD.png'])

% lp
figure(7)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zlp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_LP.png'])

% ld
figure(8)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zld))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
%print('-dpng',[ppath 'Spinup_global_LD.png'])

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
title('log10 mean biomass All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-1 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_All.png'])

% all F
figure(22)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(AllF))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass All F (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllF.png'])

% all D
figure(23)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(AllD))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass All D (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllD.png'])

% All P
figure(24)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(AllP))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass All P (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_AllP.png'])

% FracPD
figure(25)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPD)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('P:D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPD.png'])

% FracPF
figure(26)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPF)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('P:F mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPF.png'])

%% FracPFvD
figure(27)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPFvD)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('(P+F):D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFvD.png'])

% FracPDs
figure(28)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPDs)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('SP:SD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPDs.png'])

% FracPFs
figure(29)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPFs)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('SP:SF mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFs.png'])

% FracPFvDs
figure(30)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPFvDs)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('(SP+SF):SD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFvDs.png'])

% FracPDm
figure(31)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPDm)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('MP:MD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPDm.png'])

% FracPFm
figure(32)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPFm)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('MP:MF mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFm.png'])

% FracPFvDm
figure(33)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPFvDm)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('(MP+MF):MD mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_FracPFvDm.png'])

