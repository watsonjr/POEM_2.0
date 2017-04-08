% Input data and params needed in advection-diffusion scheme

clear all
close all

% Transports path
vpath = '/Volumes/GFDL/POEM_JLD/esm2m_hist/';

% Grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_hindcast_grid_cp2D.mat')

%% number of water cells
ID = find(GRD.mask==1);
NX = length(ID);

% grid size
[ni,nj] = size(GRD.mask);
isd = 1;
jsd = 1;
ied = ni;
jed = nj;

%% define a patch to advect
bio = zeros(ni,nj);
%Global
bio = 100*ones(ni,nj);   %Global
%bio(220:240,:) = 1.0; bio(121:141,195:200) = 1.0; %Atl-Arctic
%bio(:,84:109) = 1.0e2;     %seed equator
%bio(220:240,:) = 1.0e2;    %seed Atl
%bio(59:79,:) = 1.0e2;      %seed Pac
%bio(5:25,:) = 1.0e2;       %seed Indian W
%bio(340:360,:) = 1.0e2;    %seed Indian E
%bio(:,181:200) = 1.0e2;    %seed Arctic
%bio(:,12:32) = 1.0e2;      %seed Antarctic
%bio(206:295,150:177) = 1.0e2; %w/i Natasha Atl

bio = bio .* GRD.mask;

% define time
YEARS = 1;
DAYS = 365;
tstep = 1; %time step in hours

% Files to save
cname='Global_even_dt1hr_esm2m2000_vel_b100_area';
biov = zeros(NX,DAYS*YEARS);

%% do advec-diff
% % define diffusivity
% K = 600.0;
% 
% n=0;
% for Y=1:YEARS
%     yr = num2str(Y+1988-1);
%     % Velocities
%     %load([vpath 'Vel200_ESM2Mhist_' num2str(yr) '.mat'],'u','v');
%     load([vpath 'Vel200_ESM2Mhist_2000.mat'],'u','v');
%     for DAY = 1:DAYS
%         DAY
%         n=n+1;
%         U = u(:,:,DAY); 
%         V = v(:,:,DAY);
%         bio = sub_advec_diff_vel(GRD,bio,K,U,V,ni,nj,tstep);
%         biov(:,n) = bio(ID);
%     end
% end
% 
% % Save
% %csvwrite(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_diff_' cname '.csv'],biov);
% save(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_diff_' cname '.mat'],'biov');

%% do advec only

% bio = 100*ones(ni,nj);   %Global
% bio = bio .* GRD.mask;
% 
% n=0;
% for Y=1:YEARS
%     yr = num2str(Y+1988-1);
%     % Velocities
%     %load([vpath 'Vel200_ESM2Mhist_' num2str(yr) '.mat'],'u','v');
%     load([vpath 'Vel200_ESM2Mhist_2000.mat'],'u','v');
%     for DAY = 1:DAYS
%         DAY
%         n=n+1;
%         U = u(:,:,DAY); 
%         V = v(:,:,DAY);
%         bio = sub_advec_vel(GRD,bio,U,V,ni,nj,tstep);
%         biov(:,n) = bio(ID);
%     end
% end
% 
% % Save
% %csvwrite(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_' cname '.csv'],biov);
% save(['/Volumes/GFDL/CSV/advect_tests/Matlab_adv_' cname '.mat'],'biov');

%% do diff only
%Global
% bio = 100*rand(ni,nj);   
%Gradient highest at eq
g1 = 1:100;
g2 = 100:-1:1;
bio = zeros(ni,nj);
bio(:,1:100) = repmat(g1,ni,1);
bio(:,101:200) = repmat(g2,ni,1);
%Gradient highest at poles
% g1 = 1:100;
% g2 = 100:-1:1;
% bio = zeros(ni,nj);
% bio(:,1:100) = repmat(g2,ni,1);
% bio(:,101:200) = repmat(g1,ni,1);

bio = bio .* GRD.mask;

cname='Eq_grad_dt1hr_esm2m2000_vel_b100_area';

% define diffusivity
K = 600.0;
U = zeros(ni,nj);
V = zeros(ni,nj);
n=0;
for Y=1:YEARS
    yr = num2str(Y+1988-1);
    for DAY = 1:DAYS
        DAY
        n=n+1;
        bio = sub_advec_diff_vel(GRD,bio,K,U,V,ni,nj,tstep);
        biov(:,n) = bio(ID);
    end
end

% Save
%csvwrite(['/Volumes/GFDL/CSV/advect_tests/Matlab_diff_' cname '.csv'],biov);
save(['/Volumes/GFDL/CSV/advect_tests/Matlab_diff_' cname '.mat'],'biov');
