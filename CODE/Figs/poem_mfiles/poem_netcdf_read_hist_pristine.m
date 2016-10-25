% POEM output at all locations

clear all
close all

fpath='/Volumes/GFDL/NC/Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/';

%% SP
ncid = netcdf.open([fpath 'Data_hist_pristine_sml_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([fpath 'Data_hist_pristine_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass;
% SF.clev = clev;
% SF.con = con;
% SF.DD = DD;
% SF.die = die;
% SF.egg = egg;
% SF.gamma = gamma;
% SF.nu = nu;
SF.prod = prod;
SF.rec = rec;
% SF.rep = rep;
% SF.S = S;
% SF.X = X;
% SF.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% SD
ncid = netcdf.open([fpath 'Data_hist_pristine_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
% SD.clev = clev;
% SD.con = con;
% SD.DD = DD;
% SD.die = die;
% SD.egg = egg;
% SD.gamma = gamma;
% SD.nu = nu;
SD.prod = prod;
SD.rec = rec;
% SD.rep = rep;
% SD.S = S;
% SD.X = X;
% SD.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MP
ncid = netcdf.open([fpath 'Data_hist_pristine_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
% MP.clev = clev;
% MP.con = con;
% MP.DD = DD;
% MP.die = die;
% MP.egg = egg;
% MP.gamma = gamma;
% MP.nu = nu;
MP.prod = prod;
MP.rec = rec;
% MP.rep = rep;
% MP.S = S;
% MP.X = X;
% MP.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MF
ncid = netcdf.open([fpath 'Data_hist_pristine_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
% MF.clev = clev;
% MF.con = con;
% MF.DD = DD;
% MF.die = die;
% MF.egg = egg;
% MF.gamma = gamma;
% MF.nu = nu;
MF.prod = prod;
MF.rec = rec;
% MF.rep = rep;
% MF.S = S;
% MF.X = X;
% MF.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MD
ncid = netcdf.open([fpath 'Data_hist_pristine_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
% MD.clev = clev;
% MD.con = con;
% MD.DD = DD;
% MD.die = die;
% MD.egg = egg;
% MD.gamma = gamma;
% MD.nu = nu;
MD.prod = prod;
MD.rec = rec;
% MD.rep = rep;
% MD.S = S;
% MD.X = X;
% MD.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LP
ncid = netcdf.open([fpath 'Data_hist_pristine_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
% LP.clev = clev;
% LP.con = con;
% LP.DD = DD;
% LP.die = die;
% LP.egg = egg;
% LP.gamma = gamma;
% LP.nu = nu;
LP.prod = prod;
LP.rec = rec;
% LP.rep = rep;
% LP.S = S;
% LP.X = X;
% LP.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LD
ncid = netcdf.open([fpath 'Data_hist_pristine_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
% LD.clev = clev;
% LD.con = con;
% LD.DD = DD;
% LD.die = die;
% LD.egg = egg;
% LD.gamma = gamma;
% LD.nu = nu;
LD.prod = prod;
LD.rec = rec;
% LD.rep = rep;
% LD.S = S;
% LD.X = X;
% LD.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% Benthic material
ncid = netcdf.open([fpath 'Data_hist_pristine_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

BENT.bio = biomass;
clear biomass 

%% Take means

%WHOLE TIME 1861-2005
sp_mean=nanmean(SP.bio,2);
sf_mean=nanmean(SF.bio,2);
sd_mean=nanmean(SD.bio,2);
mp_mean=nanmean(MP.bio,2);
mf_mean=nanmean(MF.bio,2);
md_mean=nanmean(MD.bio,2);
lp_mean=nanmean(LP.bio,2);
ld_mean=nanmean(LD.bio,2);
b_mean=nanmean(BENT.bio,2);

SP_prod=nanmean(SP.prod,2);
SF_prod=nanmean(SF.prod,2);
SD_prod=nanmean(SD.prod,2);
MP_prod=nanmean(MP.prod,2);
MF_prod=nanmean(MF.prod,2);
MD_prod=nanmean(MD.prod,2);
LP_prod=nanmean(LP.prod,2);
LD_prod=nanmean(LD.prod,2);

%% DECADAL
st=1:120:length(time);
en=120:120:length(time);
en(15)=time(end);

for n=1:length(st)
    sp_mean(:,n+1)=nanmean(SP.bio(:,st(n):en(n)),2);
    sf_mean(:,n+1)=nanmean(SF.bio(:,st(n):en(n)),2);
    sd_mean(:,n+1)=nanmean(SD.bio(:,st(n):en(n)),2);
    mp_mean(:,n+1)=nanmean(MP.bio(:,st(n):en(n)),2);
    mf_mean(:,n+1)=nanmean(MF.bio(:,st(n):en(n)),2);
    md_mean(:,n+1)=nanmean(MD.bio(:,st(n):en(n)),2);
    lp_mean(:,n+1)=nanmean(LP.bio(:,st(n):en(n)),2);
    ld_mean(:,n+1)=nanmean(LD.bio(:,st(n):en(n)),2);
    b_mean(:,n+1)=nanmean(BENT.bio(:,st(n):en(n)),2);
    
    SP_prod(:,n+1)=nanmean(SP.prod(:,st(n):en(n)),2);
    SF_prod(:,n+1)=nanmean(SF.prod(:,st(n):en(n)),2);
    SD_prod(:,n+1)=nanmean(SD.prod(:,st(n):en(n)),2);
    MP_prod(:,n+1)=nanmean(MP.prod(:,st(n):en(n)),2);
    MF_prod(:,n+1)=nanmean(MF.prod(:,st(n):en(n)),2);
    MD_prod(:,n+1)=nanmean(MD.prod(:,st(n):en(n)),2);
    LP_prod(:,n+1)=nanmean(LP.prod(:,st(n):en(n)),2);
    LD_prod(:,n+1)=nanmean(LD.prod(:,st(n):en(n)),2);
end

%%
save([fpath 'Means_hist_pristine_Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05.mat'],...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean');


