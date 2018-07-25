% Make GRD file for FEISTY input

clear all
close all

Cdir = '/Volumes/GFDL/NEMURO/10km/';

ncid = netcdf.open([Cdir 'wc12_avgmo_1988.nc'],'NC_NOWRITE');

[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    %eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%% Retain only water cells
ID = find(MASK(:)>0);
GRD.ID = ID;
GRD.N = length(ID);
GRD.LON = LON(ID);
GRD.LAT = LAT(ID);
GRD.Z   = H(ID);
GRD.AREA  = AREA(ID);
GRD.lmask = MASK(ID);

%% Save needed variables
save([Cdir 'gridspec_10km.mat'],'AREA','H','LAT','LON','MASK');
save([Cdir 'Data_grid_10km_hist.mat'],'GRD');
          
        