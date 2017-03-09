% Visualize output of POEM
% Pre-industrial spinup at all locations
% 100 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/';

cfile = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0100';

fpath=['/Volumes/GFDL/NC/Matlab_big_size/' cfile '/'];
ppath = [pp cfile '/'];

load([fpath 'Means_preindust_' cfile '.mat']);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

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
y = 1760+(1/12):(1/12):1860;

% Piscivore
figure(1)
subplot(4,1,1)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Preindust Pelagic Piscivores')
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
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Preindust_P_time.png'])

% Planktivore
sf_tmean=sf_tmean(1:length(y));
figure(2)
subplot(3,1,1)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Preindust Forage Fishes')
xlabel('Year')
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
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Preindust_F_time.png'])

% Detritivore
figure(3)
subplot(4,1,1)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Preindust Demersal Piscivores')
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
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Preindust_D_time.png'])

%% All size classes of all

figure(4)
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
ylim([-5 2])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
title('Preindust')
stamp(cfile)
print('-dpng',[ppath 'Preindust_all_sizes.png'])

figure(5)
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;

plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Year')
ylabel('log10 Biomass (g m^-^2)')
title(['Preindust'])
print('-dpng',[ppath 'Preindust_all_types.png'])

%% Recruitment
st=1:12:length(y);
en=12:12:length(y);
MPy = NaN*ones(100,1);
SFy = MPy;
MDy = MPy;

for n=1:100
    MPy(n) = nansum(mp_trec(st(n):en(n)));
    SFy(n) = nansum(sf_trec(st(n):en(n)));
    MDy(n) = nansum(md_trec(st(n):en(n)));
end

%%
figure(6)
subplot(3,1,2)
plot(1761:1860,log10(MPy),'b','Linewidth',2); hold on;
xlim([1761 1860])
ylabel('log10 annual Recruitment (g m^-^2)')
title('Pelagic piscivores')
stamp(cfile)

subplot(3,1,1)
plot(1761:1860,log10(SFy),'r','Linewidth',2); hold on;
xlim([1761 1860])
title('Forage fishes')

subplot(3,1,3)
plot(1761:1860,log10(MDy),'k','Linewidth',2); hold on;
xlim([1761 1860])
title('Demersal piscivores')
print('-dpng',[ppath 'Preindust_recruitment.png'])

%% Plots in space
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

Zsf(grid(:,1))=sf_mean;
Zsp(grid(:,1))=sp_mean;
Zsd(grid(:,1))=sd_mean;
Zmf(grid(:,1))=mf_mean;
Zmp(grid(:,1))=mp_mean;
Zmd(grid(:,1))=md_mean;
Zlp(grid(:,1))=lp_mean;
Zld(grid(:,1))=ld_mean;
Zb(grid(:,1))=b_mean;

%% bent
figure(50)
surf(geolon_t,geolat_t,log10(Zb)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean benthic biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2.5 0.5])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_BENT.png'])

%
mgZb = (Zb/9)*1e3;
figure(51)
surf(geolon_t,geolat_t,log10(mgZb)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean benthic biomass (mg C m^-^2)')
colormap('jet')
colorbar('h')
caxis([-0.8 2.3])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_BENT_mgC.png'])

% sp
figure(11)
surf(geolon_t,geolat_t,log10(Zsp)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Larval P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_SP.png'])

% sf
figure(12)
surf(geolon_t,geolat_t,log10(Zsf)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Larval F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_SF.png'])

% sd
figure(13)
surf(geolon_t,geolat_t,log10(Zsd)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Larval D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_SD.png'])

% mp
figure(14)
surf(geolon_t,geolat_t,log10(Zmp)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Juvenile P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_MP.png'])

% mf
figure(15)
surf(geolon_t,geolat_t,log10(Zmf)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_MF.png'])

% md
figure(16)
surf(geolon_t,geolat_t,log10(Zmd)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Juvenile D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_MD.png'])

% lp
figure(17)
surf(geolon_t,geolat_t,log10(Zlp)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_LP.png'])

% ld
figure(18)
surf(geolon_t,geolat_t,log10(Zld)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_LD.png'])

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

%% ALL
figure(21)
surf(geolon_t,geolat_t,log10(All)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean biomass All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-1 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_All.png'])

% all F
figure(22)
surf(geolon_t,geolat_t,log10(AllF)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean biomass All F (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_AllF.png'])

% all D
figure(23)
surf(geolon_t,geolat_t,log10(AllD)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean biomass All D (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_AllD.png'])

% All P
figure(24)
surf(geolon_t,geolat_t,log10(AllP)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean biomass All P (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_AllP.png'])

% FracPD
figure(25)
surf(geolon_t,geolat_t,FracPD); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 P:D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_FracPD.png'])

% FracPF
figure(26)
surf(geolon_t,geolat_t,FracPF); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 P:F mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_FracPF.png'])

%% Production
Psf=NaN*ones(ni,nj);
Psp=NaN*ones(ni,nj);
Psd=NaN*ones(ni,nj);
Pmf=NaN*ones(ni,nj);
Pmp=NaN*ones(ni,nj);
Pmd=NaN*ones(ni,nj);
Plp=NaN*ones(ni,nj);
Pld=NaN*ones(ni,nj);

Psf(grid(:,1))=sf_prod;
Psp(grid(:,1))=sp_prod;
Psd(grid(:,1))=sd_prod;
Pmf(grid(:,1))=mf_prod;
Pmp(grid(:,1))=mp_prod;
Pmd(grid(:,1))=md_prod;
Plp(grid(:,1))=lp_prod;
Pld(grid(:,1))=ld_prod;

%
Psp(Psp<=0)=NaN;
Psf(Psf<=0)=NaN;
Psd(Psd<=0)=NaN;
Pmp(Pmp<=0)=NaN;
Pmf(Pmf<=0)=NaN;
Pmd(Pmd<=0)=NaN;
Plp(Plp<=0)=NaN;
Pld(Pld<=0)=NaN;

%% sp
figure(31)
surf(geolon_t,geolat_t,log10(Psp)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Larval P production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_SP.png'])

% sf
figure(32)
surf(geolon_t,geolat_t,log10(Psf)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Larval F production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_SF.png'])

% sd
figure(33)
surf(geolon_t,geolat_t,log10(Psd)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Larval D production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_SD.png'])

% mp
figure(34)
surf(geolon_t,geolat_t,log10(Pmp)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Juvenile P production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_MP.png'])

% mf
figure(35)
surf(geolon_t,geolat_t,log10(Pmf)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Adult F production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_MF.png'])

% md
figure(36)
surf(geolon_t,geolat_t,log10(Pmd)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Juvenile D production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_MD.png'])

% lp
figure(37)
surf(geolon_t,geolat_t,log10(Plp)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Adult P production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_LP.png'])

% ld
figure(38)
surf(geolon_t,geolat_t,log10(Pld)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean Adult D production (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_LD.png'])

% Diff maps of all fish
PAll = Psp+Psf+Psd+Pmp+Pmf+Pmd+Plp+Pld;
PAllF = Psf+Pmf;
PAllP = Psp+Pmp+Plp;
PAllD = Psd+Pmd+Pld;
PAllS = Psp+Psf+Psd;
PAllM = Pmp+Pmf+Pmd;
PAllL = Plp+Pld;

% ALL
figure(39)
surf(geolon_t,geolat_t,log10(PAll)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean production All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
%caxis([-1 1])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_All.png'])

% all F
figure(40)
surf(geolon_t,geolat_t,log10(PAllF)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean production All F (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_AllF.png'])

% all D
figure(41)
surf(geolon_t,geolat_t,log10(PAllD)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean production All D (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_AllD.png'])

% All P
figure(42)
surf(geolon_t,geolat_t,log10(PAllP)); view(2); hold on;
shading flat
title('Preindustrial 1811-1860 log10 mean production All P (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-5 0])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_AllP.png'])
