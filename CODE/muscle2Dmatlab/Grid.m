% GRID DEFINITION FOR AdvTVD.m
% J Guiet
% 19-12-2016
%
% Arakawa C-grid
%
%           (xP_i+1,yP_j) _ _ _ (xV_i,yV_j+1) _ _ _ (xP_i+1,yP_j+1)
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%           (xU_i,yU_j)          (xB_i,yB_j)         (xU_i+1,yU_j)
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%                 |                                        |
%           (xP_i,yP_j) _ _ _ _  (xV_i,yV_j) _ _ _ _ (xP_i+1,yP_j)
%
%
%
% Decomposition of advection velocity
%
%                                  ^        
%                                V | 
%                                  |
%                                   -----> 
%                                     U  
%
%
% INPUTS
% x:        x_i,j       -> xB_i,j
% y:        y_i,j       -> yB_i,j
% Note that x and y have a size nx+2 * ny+2, a one node ghost
% layer all arround the domain will allow the definition of the mid points.
% The computation domain will then have a size nx*ny.
%
%
% OUTPUT
% dxB:      dxB_i,j     -> (xU_i+1,yU_j)    -   (xU_i,yU_j)
% dyB:      dyB_i,j     -> (xV_i,yV_j+1)    -   (xV_i,yV_j)
% dxU:      dxU_i,j     -> (xB_i,yB_j)      -   (xB_i-1,yB_j)
% dyU:      dyU_i,j     -> (xP_i+1,yP_j)    -   (xP_i,yP_j)
% dxV:      dxV_i,j     -> (xP_i+1,yP_j)    -   (xP_i,yP_j)
% dyV:      dyV_i,j     -> (xB_i,yB_j)      -   (xB_i,yB_j-1)


function [dxB, dyB, dxU, dyU, dxV, dyV]=Grid(x,y)

% Define grid coordinates *************************************************
% Define grid dimensions
xB=zeros(size(x)-2);                                                        % Size nx * ny
yB=zeros(size(x)-2);                                                        % Size nx * ny
xU=zeros(size(x,1)-1,size(x,2)-2);                                          % Size nx+1 * ny
yU=zeros(size(x,1)-1,size(x,2)-2);                                          % Size nx+1 * ny
xV=zeros(size(x,1)-2,size(x,2)-1);                                          % Size nx * ny+1
yV=zeros(size(x,1)-2,size(x,2)-1);                                          % Size nx * ny+1
xP=zeros(size(x,1)-1,size(x,2)-1);                                          % Size nx+1 * ny+1
yP=zeros(size(x,1)-1,size(x,2)-1);                                          % Size nx+1 * ny+1

% Set coordinates values 
for j=1:size(x,2)-2                                                         % B points 
    for i=1:size(x,1)-2
        xB(i,j)=x(i+1,j+1);
        yB(i,j)=y(i+1,j+1);
    end
end

for j=1:size(x,2)-2                                                         % rho2U points 
    for i=1:size(x,1)-1
        xU(i,j)=(x(i+1,j+1)+x(i,j+1))/2;
        yU(i,j)=y(i+1,j+1);
    end
end

for j=1:size(x,2)-1                                                         % rho2V points 
    for i=1:size(x,1)-2
        xV(i,j)=x(i+1,j+1);
        yV(i,j)=(y(i+1,j+1)+y(i+1,j))/2;
    end
end

for j=1:size(x,2)-1                                                         % rho2P points 
    for i=1:size(x,1)-1
        xP(i,j)=((x(i+1,j)+x(i,j))/2+(x(i+1,j+1)+x(i,j+1))/2)/2;
        yP(i,j)=((y(i,j+1)+y(i,j))/2+(y(i+1,j+1)+y(i+1,j))/2)/2;
    end
end


% Define spatial steps ****************************************************
% Define grid dimensions
dxB=zeros(size(x)-2);                                                       % Size nx * ny
dyB=zeros(size(x)-2);                                                       % Size nx * ny
dxU=zeros(size(x,1)-1,size(x,2)-2);                                         % Size nx+1 * ny
dyU=zeros(size(x,1)-1,size(x,2)-2);                                         % Size nx+1 * ny
dxV=zeros(size(x,1)-2,size(x,2)-1);                                         % Size nx * ny+1
dyV=zeros(size(x,1)-2,size(x,2)-1);                                         % Size nx * ny+1

% Set spatial step values in miles !
for j=1:size(x,2)-2                                                         % Steps around B points 
    for i=1:size(x,1)-2        
        dxB(i,j)=120*distance(yU(i,j),xU(i,j),yU(i+1,j),xU(i+1,j));
        dyB(i,j)=120*distance(yV(i,j),xV(i,j),yV(i,j+1),xV(i,j+1));
    end
end

for j=1:size(x,2)-2                                                         % Steps around U points
    for i=1:size(x,1)-1
        if i==1                                                             % Periodicity in longitude
            dxU1=120*distance(yU(1,j),xU(1,j),yB(1,j),xB(1,j));
            dxU2=120*distance(yB(end,j),xB(end,j),yU(end,j),xU(end,j));
            dxU(i,j)=dxU1+dxU2;
        elseif i==size(x,1)-1                                               % Periodicity in longitude
            dxU(i,j)=dxU(1,j);
        else
            dxU(i,j)=120*distance(yB(i-1,j),xB(i-1,j),yB(i,j),xB(i,j));    
        end
        dyU(i,j)=120*distance(yP(i,j),xP(i,j),yP(i,j+1),xP(i,j+1));         
    end
end

for j=1:size(x,2)-1                                                         % Steps around V points
    for i=1:size(x,1)-2
        dxV(i,j)=120*distance(yP(i,j),xP(i,j),yP(i+1,j),xP(i+1,j));
        if j==1                                                             % Periodicity in latitude
            dyV1=120*distance(yV(i,1),xV(i,1),yB(i,1),xB(i,1));
            dyV2=120*distance(yB(i,end),xB(i,end),yV(i,end),xV(i,end));
            dyV(i,j)=dyV1+dyV2;
        elseif j==size(x,2)-1                                               % Periodicity in latitude
            dyV(i,j)=dyV(i,1);
        else
            dyV(i,j)=120*distance(yB(i,j-1),xB(i,j-1),yB(i,j),xB(i,j));
        end
    end
end


end
