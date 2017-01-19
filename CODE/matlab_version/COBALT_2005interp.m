% Make a test mat file of interpolated time series from one year of COBALT
% Last year of ESM2M historical = 2005

clear all
close all

Cdir = '/Volumes/GFDL/GCM_DATA/ESM2M_hist/';
Time  = ncread([Cdir 'ocean.186101-200512.temp100.nc'],'TIME');
st_lyr = length(Time)-11;
en_lyr = length(Time);
Tp  = ncread([Cdir 'ocean.186101-200512.temp100.nc'],'TEMP100',[1 1 st_lyr],[Inf Inf 12],[1 1 1]);
Tb  = ncread([Cdir 'ocean_cobalt_btm.186101-200512.btm_temp.nc'],'btm_temp',[1 1 st_lyr],[Inf Inf 12],[1 1 1]);
Zm  = ncread([Cdir 'ocean_cobalt_biomass_100.186101-200512.nmdz_100.nc'],'nmdz_100',[1 1 st_lyr],[Inf Inf 12],[1 1 1]);
Zl  = ncread([Cdir 'ocean_cobalt_biomass_100.186101-200512.nlgz_100.nc'],'nlgz_100',[1 1 st_lyr],[Inf Inf 12],[1 1 1]);
det = ncread([Cdir 'ocean_cobalt_btm.186101-200512.fndet_btm.nc'],'fndet_btm',[1 1 st_lyr],[Inf Inf 12],[1 1 1]);
dZm = ncread([Cdir 'ocean_cobalt_miscflux_100.186101-200512.jhploss_nmdz_100.nc'],'jhploss_nmdz_100',[1 1 st_lyr],[Inf Inf 12],[1 1 1]);
dZl = ncread([Cdir 'ocean_cobalt_miscflux_100.186101-200512.jhploss_nlgz_100.nc'],'jhploss_nlgz_100',[1 1 st_lyr],[Inf Inf 12],[1 1 1]);

% index of water cells
WID = find(~isnan(Zm(:,:,1))); % spatial index of water cells
NID = length(WID); % number of water cells

% setup POEM data files
D_Tp  = zeros(NID,365);
D_Tb  = zeros(NID,365);
D_Zm  = zeros(NID,365);
D_Zl  = zeros(NID,365);
D_dZm = zeros(NID,365);
D_dZl = zeros(NID,365);
D_det = zeros(NID,365);
% D_u = zeros(NID,365);
% D_v = zeros(NID,365);

% NaN velocities
% u[find(u.==minimum(u))] = 0.0
% v[find(v.==minimum(v))] = 0.0

%% interpolate to daily resolution
for j = 1:NID
    % indexes
    [m,n] = ind2sub([360,200],WID(j)); % spatial index of water cell
    
%     % v currents in m/s 
%     Y = squeeze(V200(m,n,:));
%     yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
%     D_v(j,:) = yi;
% 
%     % u currents from m/s to m/d
%     Y = squeeze(U200(m,n,:));
%     yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
%     D_u(j,:) = yi;
    
    % pelagic temperature (from Kelvin to Celcius)
    Y = squeeze(Tp(m,n,:)) - 273;
    yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
    D_Tp(j,:) = yi;
    
    % bottom temperature (in Celcius)
    Y = squeeze(Tb(m,n,:));
    yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
    D_Tb(j,:) = yi;
    
    % medium zoo: from mol N m-2 to g(WW) m-2
    % 106/16 mol C in 1 mol N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    Y = squeeze(Zm(m,n,:));
    yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
    D_Zm(j,:) = yi * (106.0/16.0) * 12.01 * 9.0;
    
    % large zoo: from mol N m-2 to g(WW) m-2
    % 106/16 mol C in 1 mol N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    Y = squeeze(Zl(m,n,:));
    yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
    D_Zl(j,:) = yi * (106.0/16.0) * 12.01 * 9.0;
    
    % medium zoo mortality: from mol N m-2 s-1 to g(WW) m-2 d-1
    % 106/16 mol C in 1 mol N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    Y = squeeze(dZm(m,n,:));
    yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
    D_dZm(j,:) = yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
    
    % large zoo mortality: from mol N m-2 s-1 to g(WW) m-2 d-1
    % 106/16 mol C in 1 mol N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    Y = squeeze(dZl(m,n,:));
    yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
    D_dZl(j,:) = yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
    
    % detrital flux to benthos: from mol C m-2 s-1 to g(WW) m-2 d-1
    % 106/16 mol C in 1 mol N
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    Y = squeeze(det(m,n,:));
    yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
    D_det(j,:) = yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24;
    
end

% convert to single precision to save space
COBALT.Tp = D_Tp;
COBALT.Tb = D_Tb;
COBALT.Zm = D_Zm;
COBALT.Zl = D_Zl;
COBALT.dZm = D_dZm;
COBALT.dZl = D_dZl;
COBALT.det = D_det;
% COBALT.U = D_u;
% COBALT.V = D_v;

% save

di = '/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2005';
save([di '.mat'], 'COBALT');



