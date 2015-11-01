% Visualize output of POEM
% Spinup at all locations
% 30 years
% Saved as netcdf files

clear all
close all

dpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Data/NC/';
fpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

sname = 'Spinup_pristine_';

%% Load netcdf data
nname = [dpath 'Spinup_pristine_pisc.nc'];
ncid = netcdf.open(nname,'NC_NOWRITE');

% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
%biomass - locations x size classes x time(=months*years)
%S - size classes
%time - number of months simulated
%X - locations

pi_b=ncread(nname,'biomass');
pi_s=ncread(nname,'S');
time=ncread(nname,'time');
X=ncread(nname,'X');
netcdf.close(ncid);

pname = [dpath 'Spinup_pristine_plan.nc'];
ncid = netcdf.open(pname,'NC_NOWRITE');
pl_b=ncread(pname,'biomass');
pl_s=ncread(pname,'S');
netcdf.close(ncid);

dname = [dpath 'Spinup_pristine_detr.nc'];
ncid = netcdf.open(dname,'NC_NOWRITE');
de_b=ncread(dname,'biomass');
de_s=ncread(dname,'S');
netcdf.close(ncid);


%% Plots over time
% Sum over locations
pi_bsum = squeeze(sum(pi_b,1));
pl_bsum = squeeze(sum(pl_b,1));
de_bsum = squeeze(sum(de_b,1));

%Piscivore
figure(1)
subplot(1,2,1)
plot(time,pi_bsum);
xlim([time(1) time(end)])
title('Piscivore')
xlabel('Time (mo)')
ylabel('Biomass (g km^-^2)')
legend(num2str(pi_s))
subplot(2,2,2)
plot(time(1:24),pi_bsum(:,1:24),'Linewidth',2)
xlim([1 24])
subplot(2,2,4)
plot(time(336:360),pi_bsum(:,336:360),'Linewidth',2)
xlim([336 360])
print('-dpng',[fpath sname 'pisc_time.png'])

%Planktivore
figure(2)
subplot(1,2,1)
plot(time,pl_bsum)
xlim([time(1) time(end)])
title('Planktivore')
xlabel('Time (mo)')
ylabel('Biomass (g km^-^2)')
legend(num2str(pl_s))
subplot(2,2,2)
plot(time(1:24),pl_bsum(:,1:24),'Linewidth',2)
xlim([1 24])
subplot(2,2,4)
plot(time(336:360),pl_bsum(:,336:360),'Linewidth',2)
xlim([336 360])
print('-dpng',[fpath sname 'plan_time.png'])

%Detritivore
figure(3)
subplot(1,2,1)
plot(time,de_bsum)
xlim([time(1) time(end)])
title('Detritivore')
xlabel('Time (mo)')
ylabel('Biomass (g km^-^2)')
legend(num2str(de_s))
subplot(2,2,2)
plot(time(1:24),de_bsum(:,1:24),'Linewidth',2)
xlim([1 24])
subplot(2,2,4)
plot(time(336:360),de_bsum(:,336:360),'Linewidth',2)
xlim([336 360])
print('-dpng',[fpath sname 'detr_time.png'])

%% Plots in space
% Sum or average over time
%What amount of time?
t=1:360;
lyr=t((end-11):end);
%each size class
pi_mean=mean(pi_b(:,:,lyr),3);
pl_mean=mean(pl_b(:,:,lyr),3);
de_mean=mean(de_b(:,:,lyr),3);
%sum size classes then take mean
pi_ssum=sum(pi_b(:,:,lyr),2);
pl_ssum=sum(pl_b(:,:,lyr),2);
de_ssum=sum(de_b(:,:,lyr),2);
pi_smean=mean(pi_ssum,3);
pl_smean=mean(pl_ssum,3);
de_smean=mean(de_ssum,3);

load('MyColormaps.mat');
gpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
grid = csvread([gpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;
%WHY IS MIN LAT = -77.5?

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);
%Z=griddata(grid(:,2),grid(:,3),grid(:,4),X,Y);

% Mean biomass of all size classes together
mpi_smean = griddata(grid(:,2),grid(:,3),pi_smean,X,Y);
mpl_smean = griddata(grid(:,2),grid(:,3),pl_smean,X,Y);
mde_smean = griddata(grid(:,2),grid(:,3),de_smean,X,Y);
mpi_smean(mpi_smean<0) = 0;
mpl_smean(mpl_smean<0) = 0;
mde_smean(mde_smean<0) = 0;
lmpi_smean = log(mpi_smean);
lmpl_smean = log(mpl_smean);
lmde_smean = log(mde_smean);


%% Sum over size classes
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,lmpi_smean); hold on;
shading interp
colormap(jet)%(cmap_color_white0)
colorbar
caxis([-7 -3])
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('Piscivores (log biomass)')
print('-dpng',[fpath sname 'total_pisc_map.png'])

figure(5)
m_proj('miller','lat',82);
m_pcolor(X,Y,lmpl_smean); hold on;
shading interp
colormap(jet)%(cmap_color_white0)
colorbar
caxis([-125 0])
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('Planktivores (log biomass)')
print('-dpng',[fpath sname 'total_plan_map.png'])

figure(6)
m_proj('miller','lat',82);
m_pcolor(X,Y,lmde_smean); hold on;
shading interp
colormap(jet)%(cmap_color_white0)
colorbar
caxis([-200 -50])
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('Detritivores (log biomass)')
print('-dpng',[fpath sname 'total_detr_map.png'])

%% Each size classes

figure(7)
for n=1:length(pi_s);
    % Mean biomass of each size class
    mpi_mean = griddata(grid(:,2),grid(:,3),pi_mean(:,n),X,Y);
    subplot(4,3,n)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log(mpi_mean)); hold on;
    shading interp
    colormap(jet)%(cmap_color_white0)
    colorbar
    caxis([-10 -5])
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['Pisc ' num2str(n)])
end
print('-dpng',[fpath sname 'sizes__pisc_map.png'])
%
figure(8)
for n=1:length(pi_s);
    mpl_mean = griddata(grid(:,2),grid(:,3),pl_mean(:,n),X,Y);
    mpl_mean(mpl_mean<0) = 0;
    subplot(4,3,n)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log(mpl_mean)); hold on;
    shading interp
    colormap(jet)%(cmap_color_white0)
    colorbar
    caxis([-150 0])
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['Plan ' num2str(n)])
end
print('-dpng',[fpath sname 'sizes_plan_map.png'])
%
figure(9)
for n=1:length(pi_s);
    mde_mean = griddata(grid(:,2),grid(:,3),de_mean(:,n),X,Y);
    mde_mean(mde_mean<0) = 0;
    subplot(4,3,n)
    m_proj('miller','lat',82);
    m_pcolor(X,Y,log(mde_mean)); hold on;
    shading interp
    colormap(jet)%(cmap_color_white0)
    colorbar
    caxis([-200 -50])
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
    title(['Detr ' num2str(n)])
end
print('-dpng',[fpath sname 'sizes_detr_map.png'])







