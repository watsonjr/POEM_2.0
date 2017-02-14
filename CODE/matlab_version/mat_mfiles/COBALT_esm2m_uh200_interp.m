%Interpolate velocities and transports to daily

clear all
close all

% 200-m weighted mean transports and velocities
fpath='/Volumes/GFDL/GCM_DATA/ESM2M_hist/';
fname='ocean.199601-200012';
load([fpath  fname '_uh200_vh200.mat']);

% Save to Julia folder
spath = '/Volumes/GFDL/POEM_JLD/esm2m_hist/';

%%
yrs= 1996:2000;
st = 1:12:length(TIME);
en = 12:12:length(TIME);
Time=TIME-TIME(1)+15;

[ni,nj,nk]=size(u200);
%%
U = reshape(u200,ni*nj,nk);
V = reshape(v200,ni*nj,nk);
UH = reshape(uh200,ni*nj,nk);
VH = reshape(vh200,ni*nj,nk);
Tmat = repmat(Time(1:12),1,ni*nj);
Tmat = Tmat';
Dmat = repmat(1:365,ni*nj,1);
%%
for y = 1:length(yrs)
    y
    
    D_u = zeros(ni*nj,365);
    D_v = zeros(ni*nj,365);
    D_uh = zeros(ni*nj,365);
    D_vh = zeros(ni*nj,365);
    
    %% interpolate to daily resolution
    for j = 1:ni*nj
        % v currents in m/s
%         Y = V(:,st(y):en(y));
%         yi = interp2(Tmat, Y, Dmat,'linear','extrap');
%         D_v = yi;
        Y = V(j,st(y):en(y));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_v(j,:) = yi;
        
        % u currents in m/s
        Y = U(j,st(y):en(y));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_u(j,:) = yi;
        
        % u transports in m2/s
        Y = UH(j,st(y):en(y));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_uh(j,:) = yi;
        
        % v transports in m2/s
        Y = VH(j,st(y):en(y));
        yi = interp1(Time(1:12), Y, 1:365,'linear','extrap');
        D_vh(j,:) = yi;
        
    end
    
    % Put back into grid
    u = reshape(D_u,ni,nj,365);
    v = reshape(D_v,ni,nj,365);
    uh = reshape(D_uh,ni,nj,365);
    vh = reshape(D_vh,ni,nj,365);
    
    % save
    save([spath 'Vel200_ESM2Mhist_' num2str(yrs(y)) '.mat'], 'u','v','uh','vh');
end


