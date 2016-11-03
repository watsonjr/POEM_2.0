% Visualize output of POEM
% Future time period (2006-2100) at all locations
% 95 years
% Means saved in mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Means_fore_fished_' cfile '.mat']);
load([cpath 'hindcast_gridspec.mat'],'geolat_t','geolon_t');
grid = csvread([cpath 'grid_csv.csv']);

%% Pick which time period mean
% 2050-2100
sp_smean=sp_mean5000;
sf_smean=sf_mean5000;
sd_smean=sd_mean5000;
mp_smean=mp_mean5000;
mf_smean=mf_mean5000;
md_smean=md_mean5000;
lp_smean=lp_mean5000;
ld_smean=ld_mean5000;
b_smean=b_mean5000;

%% Plots in space
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

%Biomass
Zsp=griddata(grid(:,2),grid(:,3),sp_smean(:,1),X,Y);
Zsf=griddata(grid(:,2),grid(:,3),sf_smean(:,1),X,Y);
Zsd=griddata(grid(:,2),grid(:,3),sd_smean(:,1),X,Y);
Zmp=griddata(grid(:,2),grid(:,3),mp_smean(:,1),X,Y);
Zmf=griddata(grid(:,2),grid(:,3),mf_smean(:,1),X,Y);
Zmd=griddata(grid(:,2),grid(:,3),md_smean(:,1),X,Y);
Zlp=griddata(grid(:,2),grid(:,3),lp_smean(:,1),X,Y);
Zld=griddata(grid(:,2),grid(:,3),ld_smean(:,1),X,Y);
Zb=griddata(grid(:,2),grid(:,3),b_smean(:,1),X,Y);

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

%% Catch

mf_mcatch(:,12) = mf_mcatch5000;
mp_mcatch(:,12) = mp_mcatch5000;
md_mcatch(:,12) = md_mcatch5000;
lp_mcatch(:,12) = lp_mcatch5000;
ld_mcatch(:,12) = ld_mcatch5000;

mf_tcatch(:,12) = mf_tcatch5000;
mp_tcatch(:,12) = mp_tcatch5000;
md_tcatch(:,12) = md_tcatch5000;
lp_tcatch(:,12) = lp_tcatch5000;
ld_tcatch(:,12) = ld_tcatch5000;

% g/m2 --> total g
area = grid(:,5);
Amf_mcatch = mf_mcatch .* repmat(area,1,12) * 365; %mean fish catch per yr
Amp_mcatch = mp_mcatch .* repmat(area,1,12) * 365;
Amd_mcatch = md_mcatch .* repmat(area,1,12) * 365;
Alp_mcatch = lp_mcatch .* repmat(area,1,12) * 365;
Ald_mcatch = ld_mcatch .* repmat(area,1,12) * 365;

Amf_tcatch = (mf_tcatch .* repmat(area,1,12)) /50; %total fish catch per yr
Amp_tcatch = (mp_tcatch .* repmat(area,1,12)) /50;
Amd_tcatch = (md_tcatch .* repmat(area,1,12)) /50;
Alp_tcatch = (lp_tcatch .* repmat(area,1,12)) /50;
Ald_tcatch = (ld_tcatch .* repmat(area,1,12)) /50;

%
Cmp=griddata(grid(:,2),grid(:,3),Amp_mcatch(:,12),X,Y);
Cmf=griddata(grid(:,2),grid(:,3),Amf_mcatch(:,12),X,Y);
Cmd=griddata(grid(:,2),grid(:,3),Amd_mcatch(:,12),X,Y);
Clp=griddata(grid(:,2),grid(:,3),Alp_mcatch(:,12),X,Y);
Cld=griddata(grid(:,2),grid(:,3),Ald_mcatch(:,12),X,Y);

CAll = Cmp+Cmf+Cmd+Clp+Cld;
CAllF = Cmf;
CAllP = Cmp+Clp;
CAllD = Cmd+Cld;
CAllM = Cmp+Cmf+Cmd;
CAllL = Clp+Cld;

%% bent
figure(50)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zb))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean benthic biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_BENT.png'])

% sp
figure(1)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zsp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean Larval P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_SP.png'])

% sf
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zsf))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean Larval F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_SF.png'])

% sd
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zsd))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean Larval D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_SD.png'])

% mp
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zmp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean Juvenile P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_MP.png'])

% mf
figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zmf))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_plot(-111,5,'o','w','MarkerSize',10);
%m_text(-111,5,'EEP','Color','red','HorizontalAlignment','center');
title('2050-2100 log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_MF.png'])

% md
figure(6)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zmd))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean Juvenile D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_MD.png'])

% lp
figure(7)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zlp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('2050-2100 log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_LP.png'])

% ld
figure(8)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zld))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('2050-2100 log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_LD.png'])

%% Diff maps of all fish

% ALL
figure(21)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(All))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean biomass All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-1 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_All.png'])

% all F
figure(22)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(AllF))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean biomass All F (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_AllF.png'])

% all D
figure(23)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(AllD))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean biomass All D (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_AllD.png'])

% All P
figure(24)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(AllP))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean biomass All P (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_AllP.png'])

% FracPD
figure(25)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPD)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 P:D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_FracPD.png'])

% FracPF
figure(26)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPF)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 P:F mean biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_FracPF.png'])

% FracPFvD
figure(27)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(FracPFvD)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 (P+F):D mean biomass(g m^-^2)')
colormap('jet')
colorbar('h')
caxis([0 1])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_FracPFvD.png'])

%% Fish catch
% mp
figure(51)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Cmp*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean Juvenile P annual catch (MT)')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_MP.png'])

% mf
figure(52)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Cmf*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_plot(-111,5,'o','w','MarkerSize',10);
%m_text(-111,5,'EEP','Color','red','HorizontalAlignment','center');
title('2050-2100 log10 mean Adult F annual catch (MT)')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_MF.png'])

% md
figure(53)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Cmd*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean Juvenile D annual catch (MT)')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_MD.png'])

% lp
figure(54)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Clp*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('2050-2100 log10 mean Adult P annual catch (MT)')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_LP.png'])

% ld
figure(55)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Cld*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('2050-2100 log10 mean Adult D annual catch (MT)')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_LD.png'])

% ALL
%sum All = 1.94e8 MT
figure(56)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(CAll*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean annual catch (MT) All Fishes')
colormap('jet')
colorbar('h')
caxis([1 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_All.png'])

% all F
figure(57)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(CAllF*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean annual catch (MT) All F')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_AllF.png'])

% all D
figure(58)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(CAllD*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean annual catch (MT) All D')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_AllD.png'])

% All P
figure(59)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(CAllP*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean annual catch (MT) All P')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_AllP.png'])

% all M
figure(60)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(CAllM*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean annual catch (MT) All M')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_AllM.png'])

% All L
figure(61)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(CAllL*1e-6))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('2050-2100 log10 mean annual catch (MT) All L')
colormap('jet')
colorbar('h')
caxis([0 5])
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_global_catch_AllL.png'])

%% N Pac

% regrid
gid = grid(:,1);

Zsf2 = NaN*ones(360,200);
Zsp2 = Zsf2;
Zsd2 = Zsf2;
Zmf2 = Zsf2;
Zmp2 = Zsf2;
Zmd2 = Zsf2;
Zlp2 = Zsf2;
Zld2 = Zsf2;
Zb2 = Zsf2;
Cmf2 = Zsf2;
Cmp2 = Zsf2;
Cmd2 = Zsf2;
Clp2 = Zsf2;
Cld2 = Zsf2;

Zsp2(gid)=sp_smean(:,1);
Zsf2(gid)=sf_smean(:,1);
Zsd2(gid)=sd_smean(:,1);
Zmp2(gid)=mp_smean(:,1);
Zmf2(gid)=mf_smean(:,1);
Zmd2(gid)=md_smean(:,1);
Zlp2(gid)=lp_smean(:,1);
Zld2(gid)=ld_smean(:,1);
Zb2(gid)=b_smean(:,1);

All2 = Zsp2+Zsf2+Zsd2+Zmp2+Zmf2+Zmd2+Zlp2+Zld2;
AllF2 = Zsf2+Zmf2;
AllP2 = Zsp2+Zmp2+Zlp2;
AllD2 = Zsd2+Zmd2+Zld2;

Cmp2(gid)=Amp_mcatch(:,12);
Cmf2(gid)=Amf_mcatch(:,12);
Cmd2(gid)=Amd_mcatch(:,12);
Clp2(gid)=Alp_mcatch(:,12);
Cld2(gid)=Ald_mcatch(:,12);

CAll2 = Cmp2+Cmf2+Cmd2+Clp2+Cld2;
CAllF2 = Cmf2;
CAllP2 = Cmp2+Clp2;
CAllD2 = Cmd2+Cld2;
CAllM2 = Cmp2+Cmf2+Cmd2;
CAllL2 = Clp2+Cld2;

%% plot info
GEOLON_T = double(geolon_t);
GEOLAT_T = double(geolat_t);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-180;
plotmaxlon=180;
latlim=[0 80];
lonlim=[plotminlon plotmaxlon];
% ENTER -100 TO MAP ORIGIN LONG

%% ALL
figure(31)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(All2)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
hcb = colorbar('h');
ylim(hcb,[-1 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean biomass All Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_All.png'])

% all F
figure(32)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(AllF2)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean biomass Forage Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_AllF.png'])

% all P
figure(33)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(AllP2)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean biomass Pelagic Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_AllP.png'])

% All D
figure(34)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(AllD2)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-2 1]);
hcb = colorbar('h');
ylim(hcb,[-2 1])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean biomass Demersal Fishes (g m^-^2)')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_AllD.png'])

%% N Pac catch
% All
figure(41)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(CAll2*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 5]);
hcb = colorbar('h');
ylim(hcb,[1 5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean annual catch (MT) All Fishes')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_catch_All.png'])

%% all F
figure(42)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(CAllF2*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 5]);
hcb = colorbar('h');
ylim(hcb,[0 5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean annual catch (MT) Forage Fishes')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_catch_AllF.png'])

% all P
figure(43)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(CAllP2*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 5]);
hcb = colorbar('h');
ylim(hcb,[0 5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean annual catch (MT) Pelagic Fishes')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_catch_AllP.png'])

% all D
figure(44)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(CAllD2*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 5]);
hcb = colorbar('h');
ylim(hcb,[0 5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean annual catch (MT) Demersal Fishes')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_catch_AllD.png'])

% all M
figure(45)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(CAllM2*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 5]);
hcb = colorbar('h');
ylim(hcb,[0 5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean annual catch (MT) Medium Fishes')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_catch_AllM.png'])

% all L
figure(46)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',[-255 -60],'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(geolat_t,geolon_t,real(log10(CAllL2*1e-6)))
colormap('jet')              
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 5]);
hcb = colorbar('h');
ylim(hcb,[0 5])                   %Set color axis if needed
set(gcf,'renderer','painters')
title('2050-2100 log10 mean annual catch (MT) Large Fishes')
stamp(cfile)
print('-dpng',[ppath 'Fore_fished_Pac_catch_AllL.png'])

