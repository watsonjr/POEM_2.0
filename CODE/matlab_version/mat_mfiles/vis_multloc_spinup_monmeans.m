% Visualize output of POEM
% Spinup at 100 locations
% 150 years
% Saved as mat files

clear all
close all

Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';

cpath = [Pdrpbx 'Princeton/POEM_other/grid_cobalt/'];
pp = [Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];

cfile = 'Diff_Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE10_CC100_RE00100';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
ppath = [pp cfile '/'];
if (~isdir(ppath))
    mkdir(ppath)
end

load([fpath 'Means_spinup_fished_' cfile '.mat']);

load([Pdrpbx 'Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat'],...
    'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

[ni,nj]=size(geolon_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

%% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
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

set(groot,'defaultAxesColorOrder',cm10);

%% Plots in time
y = time;

% Piscivore
figure(1)
subplot(4,1,1)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup fished Pelagic Piscivores')
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
print('-dpng',[ppath 'Spinup_fished_P_time.png'])

% Planktivore
sf_tmean=sf_tmean(1:length(y));
figure(2)
subplot(3,1,1)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup fished Forage Fishes')
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
print('-dpng',[ppath 'Spinup_fished_F_time.png'])

% Detritivore
figure(3)
subplot(4,1,1)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup fished Demersal Piscivores')
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
print('-dpng',[ppath 'Spinup_fished_D_time.png'])

%% Benthic inverts
figure(92)
subplot(2,1,1)
plot(y,log10(b_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Spinup fished Benthic Inverts')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')

subplot(2,1,2)
plot(y,(b_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
ylabel('Biomass (g m^-^2)')
print('-dpng',[ppath 'Spinup_fished_B_time.png'])

%% All size classes of all

figure(4)
plot(y,log10(b_tmean),'Linewidth',1); hold on;
plot(y,log10(sf_tmean),'Linewidth',1); hold on;
plot(y,log10(mf_tmean),'Linewidth',1); hold on;
plot(y,log10(sp_tmean),'Linewidth',1); hold on;
plot(y,log10(mp_tmean),'Linewidth',1); hold on;
plot(y,log10(lp_tmean),'Linewidth',1); hold on;
plot(y,log10(sd_tmean),'Linewidth',1); hold on;
plot(y,log10(md_tmean),'Linewidth',1); hold on;
plot(y,log10(ld_tmean),'Linewidth',1); hold on;
legend('B','SF','MF','SP','MP','LP','SD','MD','LD')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title('Spinup')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_all_sizes.png'])

figure(5)
F = sf_tmean+mf_tmean;
P = sp_tmean+mp_tmean+lp_tmean;
D = sd_tmean+md_tmean+ld_tmean;
B = b_tmean;

plot(y,log10(B),'color',[0.5 0.5 0.5],'Linewidth',2); hold on;
plot(y,log10(F),'r','Linewidth',2); hold on;
plot(y,log10(P),'b','Linewidth',2); hold on;
plot(y,log10(D),'k','Linewidth',2); hold on;
legend('B','F','P','D')
legend('location','eastoutside')
xlim([y(1) y(end)])
ylim([-5 2])
xlabel('Time (y)')
ylabel('log10 Biomass (g m^-^2)')
title(['Spinup'])
print('-dpng',[ppath 'Spinup_fished_all_types.png'])


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

ocean=NaN*ones(ni,nj);
ocean(grid(:,1))=ones(size(sf_mean));

%% ocean cells
% figure(55)
% surf(geolon_t,geolat_t,ocean); view(2); hold on;
% shading flat
% title('Water cells')
% colormap('jet')
% colorbar('h')
% caxis([1 2])
% stamp(cfile)
% print('-dpng',[ppath 'Ocean_cells.png'])

% bent
figure(50)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(Zb))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2.5 0.5]);
hcb = colorbar('h');
ylim(hcb,[-2.5 0.5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean benthic biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_BENT.png'])

%
mgZb = (Zb/9)*1e3;
figure(51)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(mgZb))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-0.8 2.3]);
hcb = colorbar('h');
ylim(hcb,[-0.8 2.3])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean benthic biomass (mg C m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_BENT_mgC.png'])

% sp
figure(11)
surf(geolon_t,geolat_t,log10(Zsp)); view(2); hold on;
shading flat
title('log10 mean Larval P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_SP.png'])

% sf
figure(12)
surf(geolon_t,geolat_t,log10(Zsf)); view(2); hold on;
shading flat
title('log10 mean Larval F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_SF.png'])

% sd
figure(13)
surf(geolon_t,geolat_t,log10(Zsd)); view(2); hold on;
shading flat
title('log10 mean Larval D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_SD.png'])

% mp
figure(14)
surf(geolon_t,geolat_t,log10(Zmp)); view(2); hold on;
shading flat
title('log10 mean Juvenile P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_MP.png'])

% mf
figure(15)
surf(geolon_t,geolat_t,log10(Zmf)); view(2); hold on;
shading flat
title('log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_MF.png'])

% md
figure(16)
surf(geolon_t,geolat_t,log10(Zmd)); view(2); hold on;
shading flat
title('log10 mean Juvenile D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_MD.png'])

% lp
figure(17)
surf(geolon_t,geolat_t,log10(Zlp)); view(2); hold on;
shading flat
title('log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_LP.png'])

% ld
figure(18)
surf(geolon_t,geolat_t,log10(Zld)); view(2); hold on;
shading flat
title('log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_LD.png'])

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
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(All))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 2]);
hcb = colorbar('h');
ylim(hcb,[-1 2])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_All.png'])

% all F
figure(22)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllF))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All F (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_AllF.png'])

% all D
figure(23)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllD))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All D (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_AllD.png'])

% All P
figure(24)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,log10(AllP))
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('log10 mean biomass All P (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_AllP.png'])

% FracPD
figure(25)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPD)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
hcb = colorbar('h');
ylim(hcb,[0 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('P:D mean biomass(g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_FracPD.png'])

% FracPF
figure(26)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,FracPF)
colormap('jet')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
hcb = colorbar('h');
ylim(hcb,[0 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('P:F mean biomass (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Spinup_fished_global_FracPF.png'])

% %% FracPFvD
% figure(27)
% surf(geolon_t,geolat_t,FracPFvD); view(2); hold on;
% shading flat
% title('(P+F):D mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_fished_global_FracPFvD.png'])
% 
% % FracPDs
% figure(28)
% surf(geolon_t,geolat_t,FracPDs); view(2); hold on;
% shading flat
% title('SP:SD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_fished_global_FracPDs.png'])
% 
% % FracPFs
% figure(29)
% surf(geolon_t,geolat_t,FracPFs); view(2); hold on;
% shading flat
% title('SP:SF mean biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_fished_global_FracPFs.png'])
% 
% % FracPFvDs
% figure(30)
% surf(geolon_t,geolat_t,FracPFvDs); view(2); hold on;
% shading flat
% title('(SP+SF):SD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_fished_global_FracPFvDs.png'])
% 
% % FracPDm
% figure(31)
% surf(geolon_t,geolat_t,FracPDm); view(2); hold on;
% shading flat
% title('MP:MD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_fished_global_FracPDm.png'])
% 
% % FracPFm
% figure(32)
% surf(geolon_t,geolat_t,FracPFm); view(2); hold on;
% shading flat
% title('MP:MF mean biomass (g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_fished_global_FracPFm.png'])
% 
% % FracPFvDm
% figure(33)
% surf(geolon_t,geolat_t,FracPFvDm); view(2); hold on;
% shading flat
% title('(MP+MF):MD mean biomass(g m^-^2)')
% colormap('jet')
% colorbar('h')
% caxis([0 1])
% stamp(cfile)
% print('-dpng',[ppath 'Spinup_fished_global_FracPFvDm.png'])
% 
