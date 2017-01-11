% VELOCITY DEFINITION FOR AdvTVD.m
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
% u:        u_i,j       -> U velocity component at (xB_i,yB_j)
% v:        v_i,j       -> V velocity component at (xB_i,yB_j)
% Note that u, v have a size ny+2 * nx+2, a one node ghost
% layer all arround the domain will allow the definition of the mid points.
% The computation domain will then have a size ny*nx.
%
%
% OUTPUT
% Udemi:    Udemi_i,j   -> U velocity component at (xU_i,yU_j)
% Vdemi:    Vdemi_i,j   -> V velocity component at (xV_i,yV_j)


function [Udemi, Vdemi]=Velocity(u,v)

% Define velocities at cell boundary **************************************
% Define grid dimensions
Udemi=zeros(size(u,1)-1,size(u,2)-2);                                       % Size nx+1 * ny
Vdemi=zeros(size(u,1)-2,size(u,2)-1);                                       % Size nx * ny+1

% Set velocity values
for j=1:size(u,2)-2                                                         % Udemi at U points 
    for i=1:size(u,1)-1
        Udemi(i,j)=u(i,j+1)+(u(i+1,j+1)-u(i,j+1))/2;
    end
end
for j=1:size(u,2)-1                                                         % Vdemi at V points
    for i=1:size(u,1)-2
        Vdemi(i,j)=v(i+1,j)+(v(i+1,j+1)-v(i+1,j))/2;
    end
end

end
