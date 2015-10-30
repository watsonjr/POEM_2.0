% Visualize output of POEM
% Spinup at all locations
% 30 years
% Saved as netcdf files

clear all
close all

dpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Data/NC/';

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

%Planktivore
figure(2)
subplot(1,2,1)
plot(time,pl_bsum)
xlim([time(1) time(end)])
title('Planktivore')
xlabel('Time (d)')
ylabel('Biomass (g km^-^2)')
legend(num2str(pl_s))
subplot(2,2,2)
plot(time(1:24),pl_bsum(:,1:24),'Linewidth',2)
xlim([1 24])
subplot(2,2,4)
plot(time(336:360),pl_bsum(:,336:360),'Linewidth',2)
xlim([336 360])

%Detritivore
figure(3)
subplot(1,2,1)
plot(time,de_bsum)
xlim([time(1) time(end)])
title('Detritivore')
xlabel('Time (d)')
ylabel('Biomass (g km^-^2)')
legend(num2str(de_s))
subplot(2,2,2)
plot(time(1:24),de_bsum(:,1:24),'Linewidth',2)
xlim([1 24])
subplot(2,2,4)
plot(time(336:360),de_bsum(:,336:360),'Linewidth',2)
xlim([336 360])

%% Plots in space
% Sum or average over time
%What amount of time?

gpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
grid = csvread([gpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;
x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);
Z=griddata(grid(:,2),grid(:,3),grid(:,4),X,Y);


