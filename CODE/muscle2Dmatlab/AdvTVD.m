% ADVECTION SCHEME WITH FLUX LIMITER (TVD, MINMOD)
% Jérôme Guiet
% 19-12-2016
%
%
% This numerical scheme is detailed in the chapter 4 of the lectures of V. 
% Springel and C.P. Dullemond on numerical methods (see  http://www.ita.uni
% -heidelberg.de/~dullemond/lectures/num_fluid_2012/)
% It is also partly inspired by the advection scheme for ecosystem models
% detailed in Faugeras and Maury 2005, Mathematical Biosciences and
% engineering 2(4):1-23.
%
%
% Advection on a Arakawa C-grid
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
% Bin:      Bin_i,j     -> Biomass/Abundance/Energy density at (xB_i,yB_j)
% Udemi:    Udemi_i,j   -> U velocity component at (xU_i,yU_j)
% Vdemi:    Vdemi_i,j   -> V velocity component at (xV_i,yV_j)
% dxB:      dxB_i,j     -> (xU_i+1,yU_j)    -   (xU_i,yU_j)
% dyB:      dyB_i,j     -> (xV_i,yV_j+1)    -   (xV_i,yV_j)
% dxU:      dxU_i,j     -> (xB_i,yB_j)      -   (xB_i-1,yB_j)
% dyU:      dyU_i,j     -> (xP_i+1,yP_j)    -   (xP_i,yP_j)
% dxV:      dxV_i,j     -> (xP_i+1,yP_j)    -   (xP_i,yP_j)
% dyV:      dyV_i,j     -> (xB_i,yB_j)      -   (xB_i,yB_j-1)  
% dt:       Simulation time step
% mask:     mask for dry and wet grid points
%
%
% OUTPUTS
% Bout:     Bout_i,j     -> Biomass/Abundance/Energy density at (xB_i,yB_j)


function [Bout]=AdvTVD(Bin,Udemi,Vdemi,dxB,dyB,dxU,dyU,dxV,dyV,dt,mask)

% Computation matrix ******************************************************
% Two ghost layers are added on the computational domain in order simplify 
% the computation. The Bin, dxU and dyV, mask are extended.
% Define matrix dimension
Bcomp=zeros(size(Bin)+4);                                                   % size nx+4 * ny+4
dxUcomp=zeros(size(dxU)+2);                                                 % Size nx+3 * ny+2
dyVcomp=zeros(size(dyV)+2);                                                 % Size nx+2 * ny+3
maskcomp=zeros(size(Bin)+4);                                                % size nx+4 * ny+4

% Set extended matrix
for j=1:size(Bin,2)                                                         % Bcomp for Bin 
    for i=1:size(Bin,1)  
        Bcomp(i+2,j+2)=Bin(i,j);
    end
end
for i=1:size(Bin,1)                                                         % Bcomp periodicity in longitude
    Bcomp(i,1)=Bcomp(i,end-3);
    Bcomp(i,2)=Bcomp(i,end-2);
    Bcomp(i,end-1)=Bcomp(i,3);
    Bcomp(i,end)=Bcomp(i,4);    
end
for j=1:size(Bin,2)                                                         % Bcomp periodicity in latitude
    Bcomp(1,j)=Bcomp(end-3,j);                                              
    Bcomp(2,j)=Bcomp(end-2,j);
    Bcomp(end-1,j)=Bcomp(3,j);
    Bcomp(end,j)=Bcomp(4,j);
end

for j=1:size(Bin,2)                                                         % dxUcomp for dxU 
    for i=1:size(Bin,1)+1  
        dxUcomp(i+1,j+1)=dxU(i,j);
    end
end
for i=1:size(Bin,1)+1                                                       % dxUcomp periodicity in longitude
    dxUcomp(i,1)=dxUcomp(i,end-1);
    dxUcomp(i,end)=dxUcomp(i,2);
end
for j=1:size(Bin,2)                                                         % dxUcomp periodicity in latitude
    dxUcomp(1,j)=dxUcomp(end-1,j);
    dxUcomp(end,j)=dxUcomp(2,j);
end

for j=1:size(Bin,2)+1                                                       % dyVcomp for dyV 
    for i=1:size(Bin,1)  
        dyVcomp(i+1,j+1)=dyV(i,j);
    end
end
for i=1:size(Bin,1)                                                         % dyVcomp periodicity in longitude
    dyVcomp(i,1)=dyVcomp(i,end-1);
    dyVcomp(i,end)=dyVcomp(i,2);
end
for j=1:size(Bin,2)+1                                                       % dyVcomp periodicity in latitude
    dyVcomp(1,j)=dyVcomp(end-1,j);
    dyVcomp(end,j)=dyVcomp(2,j);
end

for j=1:size(Bin,2)                                                         % maskcomp for mask 
    for i=1:size(Bin,1)  
        maskcomp(i+2,j+2)=mask(i,j);
    end
end
for i=1:size(Bin,1)                                                         % maskcomp periodicity in longitude
    maskcomp(i,1)=maskcomp(i,end-3);
    maskcomp(i,2)=maskcomp(i,end-2);
    maskcomp(i,end-1)=maskcomp(i,3);
    maskcomp(i,end)=maskcomp(i,4);    
end
for j=1:size(Bin,2)                                                         % maskcomp periodicity in latitude
    maskcomp(1,j)=maskcomp(end-3,j);                                              
    maskcomp(2,j)=maskcomp(end-2,j);
    maskcomp(end-1,j)=maskcomp(3,j);
    maskcomp(end,j)=maskcomp(4,j);
end


% Compute flux at boundaries of cells *************************************
% Define flux matrix
Fxdemi=zeros(size(Bin,1)+1,size(Bin,2));                                    % Size nx+1 *ny
Fydemi=zeros(size(Bin,1),size(Bin,2)+1);                                    % Size nx * ny+1

% Set flux
for j=1:size(Bin,2)                                                         % Fxdemi at U points
    for i=1:size(Bin,1)+1
        if (maskcomp(i+2,j+2)==1)                                           % Test for dry and wet points     
            if (Udemi(i,j)>0)                                               % Test for advection direction
                if (maskcomp(i+1,j+2)==0)                   
                    % Flux nul
                    Fxdemi(i,j)=0;
                elseif (maskcomp(i,j+2)==0)
                    % Flux first order
                    Fxdemi(i,j)=dyU(i,j)*Udemi(i,j)         *Bcomp(i+1,j+2);
                else                
                    % Flux correction (here minmod)
                    fluxcor= minmod(  (Bcomp(i+1,j+2)-Bcomp(i,j+2))/dxUcomp(i,j+1),        (Bcomp(i+2,j+2)-Bcomp(i+1,j+2))/dxUcomp(i+1,j+1)); 
                    % Flux computation      
                    Fxdemi(i,j)=dyU(i,j)*Udemi(i,j)         *Bcomp(i+1,j+2)  +...
                                dyU(i,j)*abs(Udemi(i,j))/2  *fluxcor         *(dxUcomp(i+1,j+1)-abs(Udemi(i,j)*dt));
                end
            else
                if (maskcomp(i+3,j+2)==0)
                    % Flux first order
                    Fxdemi(i,j)=dyU(i,j)*Udemi(i,j)         *Bcomp(i+2,j+2);
                else                
                    % Flux correction (here minmod)
                    fluxcor= minmod(  (Bcomp(i+2,j+2)-Bcomp(i+1,j+2))/dxUcomp(i+1,j+1),    (Bcomp(i+3,j+2)-Bcomp(i+2,j+2))/dxUcomp(i+2,j+1)); 
                    % Flux computation      
                    Fxdemi(i,j)=dyU(i,j)*Udemi(i,j)         *Bcomp(i+2,j+2)  +...
                                dyU(i,j)*abs(Udemi(i,j))/2  *fluxcor         *(dxUcomp(i+1,j+1)-abs(Udemi(i,j)*dt));
                end
            end
        else
            % Flux nul
            Fxdemi(i,j)=0;
        end
    end
end

for j=1:size(Bin,2)+1                                                       % Fydemi at V points
    for i=1:size(Bin,1)
        if (maskcomp(i+2,j+2)==1)                                           % Test for dry and wet points
            if (Vdemi(i,j)>0)                                               % Test for advection direction
                if (maskcomp(i+2,j+1)==0)
                    % Flux nul
                    Fydemi(i,j)=0;
                elseif (maskcomp(i+2,j)==0)
                    % Flux first order
                    Fydemi(i,j)=dxV(i,j)*Vdemi(i,j)         *Bcomp(i+2,j+1);
                else
                    % Flux correction (here minmod)
                    fluxcor= minmod(  (Bcomp(i+2,j+1)-Bcomp(i+2,j))/dyVcomp(i+1,j),      (Bcomp(i+2,j+2)-Bcomp(i+2,j+1))/dyVcomp(i+1,j+1));
                    % Flux computation
                    Fydemi(i,j)=dxV(i,j)*Vdemi(i,j)         *Bcomp(i+2,j+1)  +...
                                dxV(i,j)*abs(Vdemi(i,j))/2  *fluxcor         *(dyVcomp(i+1,j+1)-abs(Vdemi(i,j)*dt));
                end
            else
                if (maskcomp(i+2,j+3)==0)
                    % Flux first order
                    Fydemi(i,j)=dxV(i,j)*Vdemi(i,j)         *Bcomp(i+2,j+2);
                else
                    % Flux correction (here minmod)
                    fluxcor= minmod(  (Bcomp(i+2,j+2)-Bcomp(i+2,j+1))/dyVcomp(i+1,j+1),  (Bcomp(i+2,j+3)-Bcomp(i+2,j+2))/dyVcomp(i+1,j+2));
                    % Flux computation
                    Fydemi(i,j)=dxV(i,j)*Vdemi(i,j)         *Bcomp(i+2,j+2)  +...
                                dxV(i,j)*abs(Vdemi(i,j))/2  *fluxcor         *(dyVcomp(i+1,j+1)-abs(Vdemi(i,j)*dt));
                end
            end
        else
            % Flux nul
            Fydemi(i,j)=0;
        end
    end
end

% Update biomass **********************************************************
for j=1:size(Bin,2)                                                     
    for i=1:size(Bin,1)
        Bout(i,j)=Bin(i,j)  -   dt/dxB(i,j)/dyB(i,j)    *((Fxdemi(i+1,j)-Fxdemi(i,j))+(Fydemi(i,j+1)-Fydemi(i,j)));
    end
end

end




