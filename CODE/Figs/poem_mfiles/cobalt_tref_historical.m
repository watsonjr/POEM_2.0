% FIND MIN T AT EACH ID DURING HISTORICAL PERIOD
% USE SAME TIME FRAME AS R ASCH METHODS 1901-1950
% Temps from hindcast netcdfs
% 'ocean.186101-200512.temp_100_avg.nc'

clear all
close all

hpath='/Volumes/GFDL/GCM_DATA/Hindcast/';

%% Hindcast
%Pel temp
ncid = netcdf.open([hpath 'ocean.186101-200512.temp_100_avg.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%% Bottom temp
ncid = netcdf.open([hpath 'ocean.186101-200512.bottom_temp.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
ptemp=TEMP_100;
nn=find(ptemp<=-1e33);
ptemp(nn)=nan(size(nn));
ptemp = ptemp - 273; %Kelvin to Celcius

btemp=bottom_temp;
nb=find(btemp<=-1e19);
btemp(nb)=nan(size(nb));

%%
nt=length(TIME);
t=1861:(1/12):(2006-(1/12));
%1901 through 1950
yr=481:1080;

ptemp = ptemp(:,:,yr);
btemp = btemp(:,:,yr);
ptempv = reshape(ptemp,360*200,600);
btempv = reshape(btemp,360*200,600);

% abs min
aminp = [nanmin(ptempv')]';
aminb = [nanmin(btempv')]';

%% min of monthly means
nyr=length(yr);
id1 = 1:12:(nyr-1);
id2 = 12:12:(nyr);
ID  = [id1; id2];
mptemp=NaN*ones(360*200,50);
mbtemp=NaN*ones(360*200,50);
for m=1:50
    mptemp(:,m) = nanmean(ptempv(:,id1(m):id2(m)),2);
    mbtemp(:,m) = nanmean(btempv(:,id1(m):id2(m)),2);
end
mminp = [nanmin(mptemp')]';
mminb = [nanmin(mbtemp')]';

%% Compare
aminpG = reshape(aminp,360,200);
aminbG = reshape(aminb,360,200);
mminpG = reshape(mminp,360,200);
mminbG = reshape(mminb,360,200);

%%
figure(1)
pcolor(aminpG')
colorbar
colormap(jet)
caxis([-2 25])
title('Abs min pel')

figure(2)
pcolor(mminpG')
colorbar
colormap(jet)
caxis([-2 25])
title('Clim mean min pel')

figure(3)
pcolor(aminbG')
colorbar
colormap(jet)
caxis([-2 10])
title('Abs min bot')

figure(4)
pcolor(mminbG')
colorbar
colormap(jet)
caxis([-2 10])
title('Clim mean min bot')

%% Add to phenol data
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/phenol/';
load([dpath 'data_for_colleen.mat']);
gpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
grid = csvread([gpath 'grid_csv.csv']);
gridID=csvread([gpath 'grid_200x300_phenolID.csv'],1,0);
load([gpath 'gridspec_forecast.mat'],'tmask');

% Map phenol params
Dt=NaN*ones(size(geolat_t));
Dt(:,43:200)=degree_days_mean';

T0=NaN*ones(size(geolat_t));
T0(:,43:200)=sst_min_baseline';

%% Land
surf_tmask = tmask(:,:,1);
lmask = surf_tmask;
lmask(lmask==0) = 999;
lmask(lmask==1) = NaN;

[nlon,nlat] = size(geolon_t);
t1 = [1:nlon*nlat]';
t1(:,2) = geolon_t(:);
t1(:,3) = geolat_t(:);
t1(:,4) = surf_tmask(:);
t1(:,5) = Dt(:);
t1(:,6) = T0(:);
t1(:,7) = mminp(:);
t1(:,8) = mminb(:);
t1(:,9) = aminp(:);
t1(:,10) = aminb(:);

ocean = find(t1(:,4)==1);
t2 = t1(ocean,:);

%% Georges Bank (NEP)
lon=find(t2(:,2)<=-66 & t2(:,2)>=-67);
lat=find(t2(:,3)<=42 & t2(:,3)>=41);
gid=intersect(lon,lat);

% Eastern Bering Sea (Pribs)
lon=find(t2(:,2)<=-169 & t2(:,2)>=-170);
lat=find(t2(:,3)<=57 & t2(:,3)>=56);
eid=intersect(lon,lat);

% Subarctic Pacific Gyre (Ocean Station Papa)
lon=find(t2(:,2)<=-145 & t2(:,2)>=-146);
lat=find(t2(:,3)<=51 & t2(:,3)>=50);
pid=intersect(lon,lat);

% HOT
lon=find(t2(:,2)<=-157 & t2(:,2)>=-158);
lat=find(t2(:,3)<=23 & t2(:,3)>=22);
hid=intersect(lon,lat);

% BATS
lon=find(t2(:,2)<=-64 & t2(:,2)>=-65);
lat=find(t2(:,3)<=32 & t2(:,3)>=31);
bid=intersect(lon,lat);

% North Sea
lon=find(t2(:,2)<=4 & t2(:,2)>=3);
lat=find(t2(:,3)<=57 & t2(:,3)>=56);
nid=intersect(lon,lat);

%% Save
names={'Georges Bank','Eastern Bering Sea','Ocean Station Papa',...
    'HOT','BATS','North Sea'};

ids(1,1)=gid;
ids(2,1)=eid;
ids(3,1)=pid;
ids(4,1)=hid;
ids(5,1)=bid;
ids(6,1)=nid;

lons = t2(ids,2);
lats = t2(ids,3);
DDs  = t2(ids,5);
T0s  = t2(ids,6);
cT0p  = t2(ids,7);
cT0b  = t2(ids,8);
aT0p  = t2(ids,9);
aT0b  = t2(ids,10);

%%
T=table(names',ids,lons,lats,DDs,T0s,cT0p,cT0b,aT0p,aT0b,...
    'VariableNames',{'Location','ID','Lon','Lat','Dthresh',...
    'Tref','cTminP','cTminB','aTminP','aTminB'});
writetable(T,[gpath 'grid_360x200_id_locs_clim_min.csv'],'Delimiter',',');
save([gpath 'grid_360x200_id_locs_clim_min.mat'],'T');

%%
G=table(t2(:,1),t2(:,2),t2(:,3),t2(:,5),t2(:,6),...
    t2(:,7),t2(:,8),'VariableNames',{'ID','Lon','Lat',...
    'Dthresh','Tref','TminP','TminB'});
writetable(G,[gpath 'grid_360x200_phenolID_clim_min.csv'],'Delimiter',',');

%%
Tp = t2(:,7);
Tb = t2(:,8);
csvwrite([gpath 'grid_phenol_T0p_clim_min_NOflip.csv'],Tp)
csvwrite([gpath 'grid_phenol_T0b_clim_min_NOflip.csv'],Tb)




