% Test function for 2D advection of biomass
clear all

% Create dummy fields *****************************************************
% Biomass distrib
B=zeros(362,52);
B(51:81,21:31)=10;

% Domain mask
mask=ones(360,50);
% mask(randi([1,360*50],360*10,1))=0
% mask(70:90,20:30)=0;
mask(:,1)=0;
mask(:,end)=0;

% Coordinates
lat0=[-25.5:1:25.5];
lon0=[-180.5:1:180.5]';
lon=repmat(lon0,[1,52]);
lat=repmat(lat0,[362,1]);

% Velocities
u=ones(362,52).*cosd(lon)*30;
v=ones(362,52).*cosd(lat)*10;


% GET GRID ****************************************************************
% Determines the features of the Arakawa C-grid used for the computational 
% domain 
[dxB, dyB, dxU, dyU, dxV, dyV]=Grid(lon,lat);

% ITERATE *****************************************************************
B(2:end-1,2:end-1)=B(2:end-1,2:end-1)./dxB./dyB.*mask;                         % Normalize at grid cell                                              
for t=1:100                                                                    % Loop on time
    % GET VELOCITIES ******************************************************
    % Velocity at cell boundary at each time step
    [Udemi, Vdemi]=Velocity(u,v);
    % ADVECTION ***********************************************************
    % Advection of biomass B in the domain (LOOP OVER SIZE CLASS)
    [B(2:end-1,2:end-1)]=AdvTVD(B(2:end-1,2:end-1),Udemi,Vdemi,dxB, dyB, dxU, dyU, dxV, dyV, 1, mask);   
    Btot(t)=sum(sum(B(2:end-1,2:end-1).*dxB.*dyB));                            % Compute Btot to test conservation
end

