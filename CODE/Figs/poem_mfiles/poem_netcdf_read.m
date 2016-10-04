% POEM output at all locations

clear all
close all

fpath='/Volumes/GFDL/NC/Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10/';

%% SP
ncid = netcdf.open([fpath 'Data_spinup_pristine_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
SP.clev = clev;
SP.con = con;
SP.DD = DD;
SP.die = die;
SP.egg = egg;
SP.gamma = gamma;
SP.nu = nu;
SP.prod = prod;
SP.rec = rec;
SP.rep = rep;
SP.S = S;
SP.X = X;
SP.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

%% SF
ncid = netcdf.open([fpath 'Data_spinup_pristine_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass;
SF.clev = clev;
SF.con = con;
SF.DD = DD;
SF.die = die;
SF.egg = egg;
SF.gamma = gamma;
SF.nu = nu;
SF.prod = prod;
SF.rec = rec;
SF.rep = rep;
SF.S = S;
SF.X = X;
SF.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% SD
ncid = netcdf.open([fpath 'Data_spinup_pristine_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.clev = clev;
SD.con = con;
SD.DD = DD;
SD.die = die;
SD.egg = egg;
SD.gamma = gamma;
SD.nu = nu;
SD.prod = prod;
SD.rec = rec;
SD.rep = rep;
SD.S = S;
SD.X = X;
SD.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MP
ncid = netcdf.open([fpath 'Data_spinup_pristine_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.clev = clev;
MP.con = con;
MP.DD = DD;
MP.die = die;
MP.egg = egg;
MP.gamma = gamma;
MP.nu = nu;
MP.prod = prod;
MP.rec = rec;
MP.rep = rep;
MP.S = S;
MP.X = X;
MP.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MF
ncid = netcdf.open([fpath 'Data_spinup_pristine_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.clev = clev;
MF.con = con;
MF.DD = DD;
MF.die = die;
MF.egg = egg;
MF.gamma = gamma;
MF.nu = nu;
MF.prod = prod;
MF.rec = rec;
MF.rep = rep;
MF.S = S;
MF.X = X;
MF.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MD
ncid = netcdf.open([fpath 'Data_spinup_pristine_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.clev = clev;
MD.con = con;
MD.DD = DD;
MD.die = die;
MD.egg = egg;
MD.gamma = gamma;
MD.nu = nu;
MD.prod = prod;
MD.rec = rec;
MD.rep = rep;
MD.S = S;
MD.X = X;
MD.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LP
ncid = netcdf.open([fpath 'Data_spinup_pristine_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.clev = clev;
LP.con = con;
LP.DD = DD;
LP.die = die;
LP.egg = egg;
LP.gamma = gamma;
LP.nu = nu;
LP.prod = prod;
LP.rec = rec;
LP.rep = rep;
LP.S = S;
LP.X = X;
LP.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LD
ncid = netcdf.open([fpath 'Data_spinup_pristine_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.clev = clev;
LD.con = con;
LD.DD = DD;
LD.die = die;
LD.egg = egg;
LD.gamma = gamma;
LD.nu = nu;
LD.prod = prod;
LD.rec = rec;
LD.rep = rep;
LD.S = S;
LD.X = X;
LD.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

%% Benthic material
ncid = netcdf.open([fpath 'Data_spinup_pristine_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

BENT.bio = biomass;
clear biomass 

%%
save([fpath 'Data_spinup_pristine_Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10.mat'])


