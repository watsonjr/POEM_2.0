% Visualize output of POEM
% Spinup at 100 locations
% 50 years
% Saved as mat files

% clear all
% close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Mat_runs/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC275_RE0500';

fpath=['/Volumes/GFDL/NC/Matlab_runs/' cfile '/'];
ppath = [pp cfile '/'];

% load([fpath 'Means_spinup_' cfile '.mat']);

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

%% Planktivore
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

%% DEBUGGING -------------------------------------------------------------
% Bio
for n=46:2:52
    close all
    
    Zsp=griddata(grid(:,2),grid(:,3),SP.bio(:,n),X,Y);
    Zsf=griddata(grid(:,2),grid(:,3),SF.bio(:,n),X,Y);
    Zsd=griddata(grid(:,2),grid(:,3),SD.bio(:,n),X,Y);
    Zmp=griddata(grid(:,2),grid(:,3),MP.bio(:,n),X,Y);
    Zmf=griddata(grid(:,2),grid(:,3),MF.bio(:,n),X,Y);
    Zmd=griddata(grid(:,2),grid(:,3),MD.bio(:,n),X,Y);
    Zlp=griddata(grid(:,2),grid(:,3),LP.bio(:,n),X,Y);
    Zld=griddata(grid(:,2),grid(:,3),LD.bio(:,n),X,Y);
    Zb=griddata(grid(:,2),grid(:,3),BENT.bio(:,n),X,Y);
    
    figure(50)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zb))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean benthic biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_BENT' num2str(n) '.png'])
    
    % sp
    figure(1)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zsp))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval P biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SP' num2str(n) '.png'])
    
    % sf
    figure(2)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zsf))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval F biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SF' num2str(n) '.png'])
    
    % sd
    figure(3)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zsd))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval D biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SD' num2str(n) '.png'])
    
    % mp
    figure(4)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zmp))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile P biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MP' num2str(n) '.png'])
    
    % mf
    figure(5)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zmf))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    %m_plot(-111,5,'o','w','MarkerSize',10);
    m_text(-111,5,'EEP','Color','red','HorizontalAlignment','center');
    title(['log10 mean Adult F biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MF' num2str(n) '.png'])
    
    % md
    figure(6)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zmd))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile D biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MD' num2str(n) '.png'])
    
    % lp
    figure(7)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zlp))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
    title(['log10 mean Adult P biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LP' num2str(n) '.png'])
    
    % ld
    figure(8)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zld))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    m_text(-111,5,'EEP','Color','black','HorizontalAlignment','center');
    title(['log10 mean Adult D biomass (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LD' num2str(n) '.png'])
end

%% Con
L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

mass = [M_s;M_m;M_l];

A = 4.39;
    fc = 0.2;
    f0 = 0.6;
    epsassim = 0.6;
    n = 3/4;
    
    w = logspace(-3, 5);
    
    AvailEnergy = A*w.^n;
    Consumption = A / (epsassim*(f0-fc)) * w.^n;
    Consumption2 = A / (epsassim*(f0-fc)) * mass.^n;
    
    % Consump g/g/d --> g/d --> g/y
    
for n=48:2:54
    close all
    
    Zsp=griddata(grid(:,2),grid(:,3),SP.con(:,n),X,Y) .* M_s .* 365;
    Zsf=griddata(grid(:,2),grid(:,3),SF.con(:,n),X,Y) .* M_s .* 365;
    Zsd=griddata(grid(:,2),grid(:,3),SD.con(:,n),X,Y) .* M_s .* 365;
    Zmp=griddata(grid(:,2),grid(:,3),MP.con(:,n),X,Y) .* M_m .* 365;
    Zmf=griddata(grid(:,2),grid(:,3),MF.con(:,n),X,Y) .* M_m .* 365;
    Zmd=griddata(grid(:,2),grid(:,3),MD.con(:,n),X,Y) .* M_m .* 365;
    Zlp=griddata(grid(:,2),grid(:,3),LP.con(:,n),X,Y) .* M_l .* 365;
    Zld=griddata(grid(:,2),grid(:,3),LD.con(:,n),X,Y) .* M_l .* 365;
    
    % sp
    figure(1)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zsp))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval P con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SPcon' num2str(n) '.png'])
    
    % sf
    figure(2)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zsf))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval F con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SFcon' num2str(n) '.png'])
    
    % sd
    figure(3)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zsd))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval D con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SDcon' num2str(n) '.png'])
    
    % mp
    figure(4)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zmp))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile P con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MPcon' num2str(n) '.png'])
    
    % mf
    figure(5)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zmf))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    %m_plot(-111,5,'o','w','MarkerSize',10);
    title(['log10 mean Adult F con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MFcon' num2str(n) '.png'])
    
    % md
    figure(6)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zmd))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile D con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MDcon' num2str(n) '.png'])
    
    % lp
    figure(7)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zlp))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult P con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LPcon' num2str(n) '.png'])
    
    % ld
    figure(8)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zld))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult D con (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LDcon' num2str(n) '.png'])
end

%Time
clev_tmean(1,:)=mean(SF.clev(1:600),1);
clev_tmean(2,:)=mean(SP.clev,1);
clev_tmean(3,:)=mean(SD.clev,1);
clev_tmean(4,:)=mean(MF.clev,1);
clev_tmean(5,:)=mean(MP.clev,1);
clev_tmean(6,:)=mean(MD.clev,1);
clev_tmean(7,:)=mean(LP.clev,1);
clev_tmean(8,:)=mean(LD.clev,1);

figure
plot(1:600,clev_tmean)
legend('SF','SP','SD','MF','MP','MD','LP','LD')
legend('location','eastoutside')
xlim([1 600])
xlabel('Time (mo)')
ylabel('C/Cmax')
title('Spinup')
stamp(cfile)
print('-dpng',[ppath 'Spinup_all_sizes_clev.png'])

% Last year
lyr=time((end-12+1):end);
sp_mean=mean(SP.bio(:,lyr),2);
sf_mean=mean(SF.bio(:,lyr),2);
sd_mean=mean(SD.bio(:,lyr),2);
mp_mean=mean(MP.bio(:,lyr),2);
mf_mean=mean(MF.bio(:,lyr),2);
md_mean=mean(MD.bio(:,lyr),2);
lp_mean=mean(LP.bio(:,lyr),2);
ld_mean=mean(LD.bio(:,lyr),2);
b_mean=mean(BENT.bio(:,lyr),2);

%% Gamma/rep
for n=12:18
    close all
    
    Zsp=griddata(grid(:,2),grid(:,3),SP.gamma(:,n),X,Y);
    Zsf=griddata(grid(:,2),grid(:,3),SF.gamma(:,n),X,Y);
    Zsd=griddata(grid(:,2),grid(:,3),SD.gamma(:,n),X,Y);
    Zmp=griddata(grid(:,2),grid(:,3),MP.gamma(:,n),X,Y);
    Zmf=griddata(grid(:,2),grid(:,3),MF.rep(:,n),X,Y);
    Zmd=griddata(grid(:,2),grid(:,3),MD.gamma(:,n),X,Y);
    Zlp=griddata(grid(:,2),grid(:,3),LP.rep(:,n),X,Y);
    Zld=griddata(grid(:,2),grid(:,3),LD.rep(:,n),X,Y);
    
    % sp
    figure(1)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zsp))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval P gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SPgamma' num2str(n) '.png'])
    
    % sf
    figure(2)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zsf))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval F gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SFgamma' num2str(n) '.png'])
    
    % sd
    figure(3)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zsd))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval D gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SDgamma' num2str(n) '.png'])
    
    % mp
    figure(4)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zmp))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile P gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MPgamma' num2str(n) '.png'])
    
    % mf
    figure(5)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zmf))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    %m_plot(-111,5,'o','w','MarkerSize',10);
    title(['log10 mean Adult F rep (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MFrep' num2str(n) '.png'])
    
    % md
    figure(6)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zmd))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile D gamma (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MDgamma' num2str(n) '.png'])
    
    % lp
    figure(7)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zlp))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult P rep (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LPrep' num2str(n) '.png'])
    
    % ld
    figure(8)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,real(log10(Zld))); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult D rep (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    %caxis([-4 2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LDrep' num2str(n) '.png'])
end

%% Nu
for n=12:18
    close all
    
    Zsp=griddata(grid(:,2),grid(:,3),SP.nu(:,n),X,Y);
    Zsf=griddata(grid(:,2),grid(:,3),SF.nu(:,n),X,Y);
    Zsd=griddata(grid(:,2),grid(:,3),SD.nu(:,n),X,Y);
    Zmp=griddata(grid(:,2),grid(:,3),MP.nu(:,n),X,Y);
    Zmf=griddata(grid(:,2),grid(:,3),MF.nu(:,n),X,Y);
    Zmd=griddata(grid(:,2),grid(:,3),MD.nu(:,n),X,Y);
    Zlp=griddata(grid(:,2),grid(:,3),LP.nu(:,n),X,Y);
    Zld=griddata(grid(:,2),grid(:,3),LD.nu(:,n),X,Y);
    
    % sp
    figure(1)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zsp); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval P nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SPnu' num2str(n) '.png'])
    
    % sf
    figure(2)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zsf); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval F nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SFnu' num2str(n) '.png'])
    
    % sd
    figure(3)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zsd); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Larval D nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_SDnu' num2str(n) '.png'])
    
    % mp
    figure(4)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zmp); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile P nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MPnu' num2str(n) '.png'])
    
    % mf
    figure(5)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zmf); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    %m_plot(-111,5,'o','w','MarkerSize',10);
    title(['log10 mean Adult F nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MFnu' num2str(n) '.png'])
    
    % md
    figure(6)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zmd); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Juvenile D nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_MDnu' num2str(n) '.png'])
    
    % lp
    figure(7)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zlp); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult P nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LPnu' num2str(n) '.png'])
    
    % ld
    figure(8)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Zld); hold on;
    shading flat
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['log10 mean Adult D nu (g m^-^2) mo ' num2str(n)])
    colormap('jet')
    colorbar('h')
    caxis([-0.2 0.2])
    stamp(cfile)
    print('-dpng',[ppath 'Spinup_global_LDnu' num2str(n) '.png'])
end

%% Find negative values
tend=12;
bio=NaN*ones(9,tend);
con=NaN*ones(8,tend);
die=NaN*ones(8,tend);
gamma=NaN*ones(8,tend);
rec=NaN*ones(8,tend);

for n=1:tend
    %bio
    id = (SP.bio(:,n) < 0);
    bio(1,n) = sum(id);
    id = (SF.bio(:,n) < 0);
    bio(2,n) = sum(id);
    id = (SD.bio(:,n) < 0);
    bio(3,n) = sum(id);
    id = (MP.bio(:,n) < 0);
    bio(4,n) = sum(id);
    id = (MF.bio(:,n) < 0);
    bio(5,n) = sum(id);
    id = (MD.bio(:,n) < 0);
    bio(6,n) = sum(id);
    id = (LP.bio(:,n) < 0);
    bio(7,n) = sum(id);
    id = (LD.bio(:,n) < 0);
    bio(8,n) = sum(id);
    id = (BENT.bio(:,n) < 0);
    bio(9,n) = sum(id);
    %con
    id = (SP.con(:,n) < 0);
    con(1,n) = sum(id);
    id = (SF.con(:,n) < 0);
    con(2,n) = sum(id);
    id = (SD.con(:,n) < 0);
    con(3,n) = sum(id);
    id = (MP.con(:,n) < 0);
    con(4,n) = sum(id);
    id = (MF.con(:,n) < 0);
    con(5,n) = sum(id);
    id = (MD.con(:,n) < 0);
    con(6,n) = sum(id);
    id = (LP.con(:,n) < 0);
    con(7,n) = sum(id);
    id = (LD.con(:,n) < 0);
    con(8,n) = sum(id);
    %die
    id = (SP.die(:,n) < 0);
    die(1,n) = sum(id);
    id = (SF.die(:,n) < 0);
    die(2,n) = sum(id);
    id = (SD.die(:,n) < 0);
    die(3,n) = sum(id);
    id = (MP.die(:,n) < 0);
    die(4,n) = sum(id);
    id = (MF.die(:,n) < 0);
    die(5,n) = sum(id);
    id = (MD.die(:,n) < 0);
    die(6,n) = sum(id);
    id = (LP.die(:,n) < 0);
    die(7,n) = sum(id);
    id = (LD.die(:,n) < 0);
    die(8,n) = sum(id);
    %gamma/rep
    id = (SP.gamma(:,n) < 0);
    gamma(1,n) = sum(id);
    id = (SF.gamma(:,n) < 0);
    gamma(2,n) = sum(id);
    id = (SD.gamma(:,n) < 0);
    gamma(3,n) = sum(id);
    id = (MP.gamma(:,n) < 0);
    gamma(4,n) = sum(id);
    id = (MF.rep(:,n) < 0);
    gamma(5,n) = sum(id);
    id = (MD.gamma(:,n) < 0);
    gamma(6,n) = sum(id);
    id = (LP.rep(:,n) < 0);
    gamma(7,n) = sum(id);
    id = (LD.rep(:,n) < 0);
    gamma(8,n) = sum(id);
    %rec
    id = (SP.rec(:,n) < 0);
    rec(1,n) = sum(id);
    id = (SF.rec(:,n) < 0);
    rec(2,n) = sum(id);
    id = (SD.rec(:,n) < 0);
    rec(3,n) = sum(id);
    id = (MP.rec(:,n) < 0);
    rec(4,n) = sum(id);
    id = (MF.rec(:,n) < 0);
    rec(5,n) = sum(id);
    id = (MD.rec(:,n) < 0);
    rec(6,n) = sum(id);
    id = (LP.rec(:,n) < 0);
    rec(7,n) = sum(id);
    id = (LD.rec(:,n) < 0);
    rec(8,n) = sum(id);
    
end
%%
figure
plot(1:tend,bio)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD','B')
legend('location','northwest')
print('-dpng',[ppath 'biotest.png'])

figure
plot(1:tend,con)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'contest.png'])

figure
plot(1:tend,die)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'dietest.png'])

figure
plot(1:tend,gamma)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'gammatest.png'])

figure
plot(1:tend,rec)
ylim([0 10])
legend('SP','MP','LP','SF','MF','SD','MD','LD')
legend('location','northwest')
print('-dpng',[ppath 'rectest.png'])

% END DEBUGGING -------------------------------------------------------------

%%
Zsp=griddata(grid(:,2),grid(:,3),sp_mean,X,Y);
Zsf=griddata(grid(:,2),grid(:,3),sf_mean,X,Y);
Zsd=griddata(grid(:,2),grid(:,3),sd_mean,X,Y);
Zmp=griddata(grid(:,2),grid(:,3),mp_mean,X,Y);
Zmf=griddata(grid(:,2),grid(:,3),mf_mean,X,Y);
Zmd=griddata(grid(:,2),grid(:,3),md_mean,X,Y);
Zlp=griddata(grid(:,2),grid(:,3),lp_mean,X,Y);
Zld=griddata(grid(:,2),grid(:,3),ld_mean,X,Y);
Zb=griddata(grid(:,2),grid(:,3),b_mean,X,Y);

% bent
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
print('-dpng',[ppath 'Spinup_global_BENT.png'])

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
print('-dpng',[ppath 'Spinup_global_BENT_mgC.png'])

% sp
figure(11)
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
print('-dpng',[ppath 'Spinup_global_SP.png'])

% sf
figure(12)
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
print('-dpng',[ppath 'Spinup_global_SF.png'])

% sd
figure(13)
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
print('-dpng',[ppath 'Spinup_global_SD.png'])

% mp
figure(14)
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
print('-dpng',[ppath 'Spinup_global_MP.png'])

% mf
figure(15)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zmf))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Adult F biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_MF.png'])

% md
figure(16)
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
print('-dpng',[ppath 'Spinup_global_MD.png'])

% lp
figure(17)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zlp))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Adult P biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_LP.png'])

% ld
figure(18)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zld))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Adult D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-4 2])
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

