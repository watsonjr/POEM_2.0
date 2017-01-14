% Advection and diffusion using upwind scheme adapted from MOM code

clear all
close all

load('/Volumes/GFDL/GCM_DATA/CORE-forced/feb152013_run25_ocean.198801-200712_uh200_vh200.mat',...
    'uh200','vh200');
Uth_200 = uh200(:,:,1);
Vth_200 = vh200(:,:,1);

load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'AREA_OCN','dat','dxtn','dyte','ht','geolon_t','geolat_t');

%%
area = AREA_OCN;
area = (area)*510072000*1e6;
area = max(area,1);

% depth of the surface layer, 200m or less
eps = 1;
dep = min(ht,200);
dep = max(dep,eps);

% grid size
ni = 360;
nj = 200;
isd = 1;
jsd = 1;
ied = ni;
jed = nj;

% land mask
mask = zeros(ni,nj);
aa = find(ht > 0);
mask(aa) = 1;

%define a patch to advect
TF = zeros(360,200);
TF2 = zeros(360,200);
lonmin = -280;
lonmax = 80;
latmin = -90;
latmax = 90;
aa = find( (geolon_t > lonmin) & (geolon_t < lonmax) & (geolat_t > latmin) & ...
    (geolat_t < latmax) & (ht > 0) );
TF(aa) = 100*rand(size(aa));

% TF(220:240,:) = 100.0;
% TF = TF .* mask;

total_mass(1) = sum(TF(:).*area(:));

%% Following Advect_upwind_2D
dt = 3600;
ntime = 365*24;
% uvel = Uth_200;
% vvel = Vth_200;
uvel = zeros(ni,nj);
vvel = zeros(ni,nj);
K = 600.0;

fe = zeros(ni,nj);
fn = zeros(ni,nj);
dfe = zeros(ni,nj);
dfn = zeros(ni,nj);
gradTi = zeros(ni,nj);
gradTj = zeros(ni,nj);
upwind = zeros(ni,nj);
dupwind = zeros(ni,nj);
gT = NaN(ni,nj,ntime);

%% Advection loop
for n = 1:ntime
    n
    
    % Calculate biomass gradient
    %Gradient i
    for j=jsd:jed
        for i=isd:ied
            if (i == ied)
                gradTi(i,j) = (TF(isd,j) - TF(i,j)) ./ dyte(i,j) *mask(i,j)*mask(isd,j);
            else
                gradTi(i,j) = (TF(i+1,j) - TF(i,j)) ./ dyte(i,j) *mask(i,j)*mask(i+1,j);
            end
        end
    end
    %Gradient j
    for j=jsd:jed
        for i=isd:ied
            if (j < jed)
                gradTj(i,j) = (TF(i,j+1) - TF(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(i,j+1);
            else
                gradTj(i,j) = (TF(ni-i+1,j) - TF(i,j)) ./ dxtn(i,j) *mask(i,j)*mask(ni-i+1,j);
            end
        end
    end
    gradT = (gradTi + gradTj); %.* mask;
    gT(:,:,n) = gradT;
    
    diffusiv = 0.5*K;
    kpos     = diffusiv + abs(diffusiv);
    kneg     = diffusiv - abs(diffusiv);
    
    % Westward flux
    for j = jsd:jed
        for i = isd:ied
            velocity = 0.5*uvel(i,j);
            upos = velocity + abs(velocity);
            uneg = velocity - abs(velocity);
            
            % define only for ocean cells
            if (mask(i,j) > 0)
                
                if (i == ied)
                    fe(i,j) = dyte(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(isd,j))/dep(isd,j)* ...
                        mask(i,j)*mask(isd,j);
                    dfe(i,j)  = dyte(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(isd,j)) ...
                        .* mask(i,j) .* mask(isd,j);
                else
                    fe(i,j) = dyte(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i+1,j)/dep(i+1,j))* ...
                        mask(i,j)*mask(i+1,j);
                    dfe(i,j)  = dyte(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(i+1,j)) ...
                        .* mask(i,j) .* mask(i+1,j);
                end
                
            end
        end
    end
    
    % Northward flux
    for j = jsd:jed
        for i = isd:ied
            velocity = 0.5*vvel(i,j);
            upos = velocity + abs(velocity);
            uneg = velocity - abs(velocity);
            
            % define only for ocean cells
            if (mask(i,j) > 0)
                
                if (j < jed)
                    fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(i,j+1)/dep(i,j+1))* ...
                        mask(i,j)*mask(i,j+1);
                    dfn(i,j)  = dxtn(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(i,j+1)) .* ...
                        mask(i,j) .* mask(i,j+1);
                else
                    fn(i,j) = dxtn(i,j)*(upos.*TF(i,j)/dep(i,j) + uneg.*TF(ni-i+1,j)/dep(ni-i+1,j))* ...
                        mask(i,j)*mask(ni-i+1,j);
                    dfn(i,j)  = dxtn(i,j).*(kpos.*gradT(i,j) + kneg.*gradT(ni-i+1,j)) .* ...
                        mask(i,j) .* mask(ni-i+1,j);
                end
                
            end
        end
    end
    
    % Combine fluxes
    for j = jsd:jed
        for i = isd:ied
            if (j > 1)
                if (i > 1)
                    upwind(i,j) = mask(i,j).*(fe(i-1,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
                    dupwind(i,j) = mask(i,j).*(dfe(i-1,j)-dfe(i,j)+dfn(i,j-1)-dfn(i,j));
                else
                    upwind(i,j) = mask(i,j).*(fe(ied,j)-fe(i,j)+fn(i,j-1)-fn(i,j));
                    dupwind(i,j) = mask(i,j).*(dfe(ied,j)-dfe(i,j)+dfn(i,j-1)-dfn(i,j));
                end
            end
        end
    end
    
    % Update tracers
    for j = jsd:jed
        for i = isd:ied
            TF2(i,j) = TF(i,j) + (dt*upwind(i,j))/area(i,j) - (dt*dupwind(i,j))/area(i,j);
        end
    end
    total_mass(n+1) = sum(TF2(:).*area(:));
    
    % Plot, do mass balance, reset tracer fields
    aa = find(ht == 0);
    TF(aa) = -999;
    if n == 1
        figure(1)
        surf(geolon_t,geolat_t,TF); view(2); shading interp; caxis([0 100]);
        %pause
    end
    
    TF2(aa) = -999;
    
    if n == 8760
        figure(2)
        clf
        surf(geolon_t,geolat_t,TF2); view(2); shading interp;
        caxis([0 100]); %caxis([0 2e3]) to see how big instabilities are
        pdiff = 100*(total_mass(n+1) - total_mass(n))/total_mass(n);
        title(['%diff = ', num2str(pdiff,'%10.3e')]);
        %pause
    end
    
    TF = TF2;
    
end

%% Plot gradT
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

t = 1:2189:ntime;
d = t/24;
d = round(d);

for i=1:length(t)
    figure
    surf(geolon_t,geolat_t,gT(:,:,t(i)));
    view(2); shading interp;
    caxis([-4e-3 4e-3])
    colorbar
    title(['GradT day ', num2str(d(i))]);
    print('-dpng',[fpath 'POEM_diff_test_gradT_' num2str(d(i)) '.png'])
end

%% Arctic projection
for i=1:length(t)
    figure
    m_proj('stereographic','lat',90,'long',30,'radius',30);
    m_pcolor(geolon_t,geolat_t,gT(:,:,t(i)));
    shading interp
    colorbar
    colormap('jet')
    caxis([-4e-3 4e-3])
    m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['GradT day ' num2str(d(i))])
    print('-dpng',[fpath 'POEM_diff_test_gradT_arcticproj_' num2str(d(i)) '.png'])
end




