% Visualize output of POEM
% Pre-industrial spinup at all locations
% 100 years
% Saved as mat files

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Data_spinup_pristine_' cfile '.mat']);


%%
[loc,days]=size(SP.bio);
x=1:days;

lyr=x((end-12+1):end);

%temporal means of last year
sp_smean=nanmean(SP.bio(:,lyr),2);
sf_smean=nanmean(SF.bio(:,lyr),2);
sd_smean=nanmean(SD.bio(:,lyr),2);
mp_smean=nanmean(MP.bio(:,lyr),2);
mf_smean=nanmean(MF.bio(:,lyr),2);
md_smean=nanmean(MD.bio(:,lyr),2);
lp_smean=nanmean(LP.bio(:,lyr),2);
ld_smean=nanmean(LD.bio(:,lyr),2);
b_smean=nanmean(BENT.bio(:,lyr),2);

%spatial means over time
sp_tmean=nanmean(SP.bio,1);
sf_tmean=nanmean(SF.bio,1);
sd_tmean=nanmean(SD.bio,1);
mp_tmean=nanmean(MP.bio,1);
mf_tmean=nanmean(MF.bio,1);
md_tmean=nanmean(MD.bio,1);
lp_tmean=nanmean(LP.bio,1);
ld_tmean=nanmean(LD.bio,1);
b_tmean=nanmean(BENT.bio,1);

%% Plots in time
y = time;

% Piscivore
figure(1)
subplot(4,1,1)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Pre-Industrial Pelagic Piscivores')
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
print('-dpng',[ppath 'Preindust_pisc_time.png'])

% Planktivore
figure(2)
subplot(3,1,1)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Pre-Industrial Forage Fishes')
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

print('-dpng',[ppath 'Preindust_plan_time.png'])

% Detritivore
figure(3)
subplot(4,1,1)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('Pre-Industrial Demersal Piscivores')
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
print('-dpng',[ppath 'Preindust_detr_time.png'])

% All biomass in subplots
%SP
figure(4)
subplot(3,3,2)
plot(y,log10(sp_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('SP')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
stamp(cfile)

subplot(3,3,5)
plot(y,log10(mp_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('MP')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')

subplot(3,3,8)
plot(y,log10(lp_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('LP')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')

%FF
subplot(3,3,1)
plot(y,log10(sf_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('SF')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')

subplot(3,3,4)
plot(y,log10(mf_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('MF')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')

%Detritivore
subplot(3,3,3)
plot(y,log10(sd_tmean),'b','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('SD')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
stamp(cfile)

subplot(3,3,6)
plot(y,log10(md_tmean),'r','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('MD')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')

subplot(3,3,9)
plot(y,log10(ld_tmean),'k','Linewidth',1); hold on;
xlim([y(1) y(end)])
title('LD')
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
print('-dpng',[ppath 'Preindust_all_sizes_sub.png'])

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
xlim([y(1) y(end)])
xlabel('Time (mo)')
ylabel('log10 Biomass (g m^-^2)')
title('Pre-Industrial')
stamp(cfile)
print('-dpng',[ppath 'Preindust_all_sizes.png'])

%% Recruitment
sfrec_tmean=nanmean(SF.rec,1);
mprec_tmean=nanmean(MP.rec,1);
mdrec_tmean=nanmean(MD.rec,1);

st=1:12:length(y);
en=12:12:length(y);
MPy = NaN*ones(100,1);
SFy = MPy;
MDy = MPy;

for n=1:100
    MPy(n) = nansum(mprec_tmean(st(n):en(n)));
    SFy(n) = nansum(sfrec_tmean(st(n):en(n)));
    MDy(n) = nansum(mdrec_tmean(st(n):en(n)));
end

%%
figure(11)
subplot(3,1,2)
plot(1:100,log10(MPy),'b','Linewidth',2); hold on;
xlim([1 100])
ylabel('log10 annual Recruitment (g m^-^2)')
title('Pelagic piscivores')
stamp(cfile)

subplot(3,1,1)
plot(1:100,log10(SFy),'r','Linewidth',2); hold on;
xlim([1 100])
title('Forage fishes')

subplot(3,1,3)
plot(1:100,log10(MDy),'k','Linewidth',2); hold on;
xlim([1 100])
title('Demersal piscivores')
print('-dpng',[ppath 'Preindust_recruitment.png'])

%% Plots in space
grid = csvread([cpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

Zsp=griddata(grid(:,2),grid(:,3),sp_smean,X,Y);
Zsf=griddata(grid(:,2),grid(:,3),sf_smean,X,Y);
Zsd=griddata(grid(:,2),grid(:,3),sd_smean,X,Y);
Zmp=griddata(grid(:,2),grid(:,3),mp_smean,X,Y);
Zmf=griddata(grid(:,2),grid(:,3),mf_smean,X,Y);
Zmd=griddata(grid(:,2),grid(:,3),md_smean,X,Y);
Zlp=griddata(grid(:,2),grid(:,3),lp_smean,X,Y);
Zld=griddata(grid(:,2),grid(:,3),ld_smean,X,Y);
Zb=griddata(grid(:,2),grid(:,3),b_smean,X,Y);

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
print('-dpng',[ppath 'Preindust_global_BENT.png'])

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
print('-dpng',[ppath 'Preindust_global_SP.png'])

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
print('-dpng',[ppath 'Preindust_global_SF.png'])

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
print('-dpng',[ppath 'Preindust_global_SD.png'])

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
print('-dpng',[ppath 'Preindust_global_MP.png'])

% mf
figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zmf))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_plot(-111,5,'o','w','MarkerSize',10);
%m_text(-111,5,'EEP','Color','red','HorizontalAlignment','center');
title('log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_MF.png'])

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
print('-dpng',[ppath 'Preindust_global_MD.png'])

% lp
figure(7)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zlp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_LP.png'])

% ld
figure(8)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zld))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
%m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
title('log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_LD.png'])

%% Production
SP_prod=nanmean(SP.prod(:,lyr),2);
SF_prod=nanmean(SF.prod(:,lyr),2);
SD_prod=nanmean(SD.prod(:,lyr),2);
MP_prod=nanmean(MP.prod(:,lyr),2);
MF_prod=nanmean(MF.prod(:,lyr),2);
MD_prod=nanmean(MD.prod(:,lyr),2);
LP_prod=nanmean(LP.prod(:,lyr),2);
LD_prod=nanmean(LD.prod(:,lyr),2);

%
Psp=griddata(grid(:,2),grid(:,3),SP_prod,X,Y);
Psf=griddata(grid(:,2),grid(:,3),SF_prod,X,Y);
Psd=griddata(grid(:,2),grid(:,3),SD_prod,X,Y);
Pmp=griddata(grid(:,2),grid(:,3),MP_prod,X,Y);
Pmf=griddata(grid(:,2),grid(:,3),MF_prod,X,Y);
Pmd=griddata(grid(:,2),grid(:,3),MD_prod,X,Y);
Plp=griddata(grid(:,2),grid(:,3),LP_prod,X,Y);
Pld=griddata(grid(:,2),grid(:,3),LD_prod,X,Y);

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
figure(11)
m_proj('miller','lat',82);
m_pcolor(X,Y,log10(Psp)); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Larval P production (g g^-^1 m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_SP.png'])

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
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_SF.png'])

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
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_SD.png'])

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
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_MP.png'])

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
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_MF.png'])

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
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_prod_MD.png'])

% lp
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
caxis([-4 2])
print('-dpng',[ppath 'Preindust_global_prod_LP.png'])

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
caxis([-4 2])
print('-dpng',[ppath 'Preindust_global_prod_LD.png'])

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
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(All))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean biomass All Fishes (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_All.png'])

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
caxis([-2 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_AllF.png'])

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
caxis([-2 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_AllD.png'])

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
caxis([-2 2])
stamp(cfile)
print('-dpng',[ppath 'Preindust_global_AllP.png'])

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
print('-dpng',[ppath 'Preindust_global_FracPD.png'])

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
print('-dpng',[ppath 'Preindust_global_FracPF.png'])

% FracPFvD
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
print('-dpng',[ppath 'Preindust_global_FracPFvD.png'])

