% Visualize Arctic grid cells

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

grid = csvread('grid_csv.csv');
load('gridspec_forecast.mat');

[ni,nj] = size(geolon_t);
iindex = repmat([1:ni]',1,nj);
jindex = repmat(1:nj,ni,1);

%% 3. Stereographic projection of North Polar regions
% Note that coastline is drawn OVER the grid because of the order in which
% the two routines are called

figure(1)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,iindex);
colorbar
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','r');

figure(2)
m_proj('stereographic','lat',90,'long',30,'radius',30);
m_pcolor(geolon_t,geolat_t,jindex);
colorbar
m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
m_coast('patch',[.7 .7 .7],'edgecolor','r');

%% 2.  SSM/I Ice cover (data provided on a fixed grid)

figure(3)
m_proj('stereographic','latitude',90,'radius',55,'rotangle',45);

% Convert bottom and left corner points to screen coords. This
% is of course a kludge.
%[MAPX,dm]=m_ll2xy([279.26 350.03],[33.92 34.35],'clip','off');
%[dm,MAPY]=m_ll2xy([168.35 279.26],[30.98 33.92],'clip','off');

% Plot data as an image
m_pcolor(geolon_t,geolat_t,iindex);
set(gca,'ydir','normal');
%colormap([jet(100);0 0 0;1 1 1]);
m_coast('patch',[.6 .6 .6]);
m_grid('linewi',2,'tickdir','out');
h=colorbar('v');
set(get(h,'ylabel'),'string','i-index');

%%
pac=iindex(1:180,200);
atl=iindex(181:360,200);
atl=flipud(atl);
dpole=atl-pac;
npole(1:180)=dpole;
npole(181:360)=-1*flipud(dpole);
npole=npole';
csvwrite('npole_lon_shift.csv',npole)





