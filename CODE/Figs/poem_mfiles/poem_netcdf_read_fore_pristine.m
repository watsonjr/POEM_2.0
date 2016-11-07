% POEM output at all locations

clear all
close all

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05';

fpath=['/Volumes/GFDL/NC/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Data_fore_pristine2006-2085_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    if (strcmp(varname,'catch'))
        varname = 'caught';
    end
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.bio = biomass;
% SP.catch = caught;
% SP.clev = clev;
% SP.con = con;
% SP.DD = DD;
% SP.die = die;
% SP.egg = egg;
% SP.gamma = gamma;
% SP.nu = nu;
SP.prod = prod;
SP.rec = rec;
% SP.rep = rep;
% SP.S = S;
% SP.X = X;
% SP.time = time;

clear biomass clev con DD die egg gamma nu rec rep S X time prod caught

%% SF
ncid = netcdf.open([fpath 'Data_fore_pristine2006-2085_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    if (strcmp(varname,'catch'))
        varname = 'caught';
    end
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass;
% SF.catch = caught;
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

clear biomass clev con DD die egg gamma nu rec rep S X time prod caught

% SD
ncid = netcdf.open([fpath 'Data_fore_pristine2006-2085_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    if (strcmp(varname,'catch'))
        varname = 'caught';
    end
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
% SD.catch = caught;
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

clear biomass clev con DD die egg gamma nu rec rep S X time prod caught

% MP
ncid = netcdf.open([fpath 'Data_fore_pristine2006-2085_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    if (strcmp(varname,'catch'))
        varname = 'caught';
    end
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
% MP.catch = caught;
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

clear biomass clev con DD die egg gamma nu rec rep S X time prod caught

% MF
ncid = netcdf.open([fpath 'Data_fore_pristine2006-2085_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    if (strcmp(varname,'catch'))
        varname = 'caught';
    end
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
% MF.catch = caught;
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

clear biomass clev con DD die egg gamma nu rec rep S X time prod caught

% MD
ncid = netcdf.open([fpath 'Data_fore_pristine2006-2085_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    if (strcmp(varname,'catch'))
        varname = 'caught';
    end
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
% MD.catch = caught;
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

clear biomass clev con DD die egg gamma nu rec rep S X time prod caught

% LP
ncid = netcdf.open([fpath 'Data_fore_pristine2006-2085_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    if (strcmp(varname,'catch'))
        varname = 'caught';
    end
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
% LP.catch = caught;
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

clear biomass clev con DD die egg gamma nu rec rep S X time prod caught

% LD
ncid = netcdf.open([fpath 'Data_fore_pristine2006-2085_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    if (strcmp(varname,'catch'))
        varname = 'caught';
    end
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
% LD.catch = caught;
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

clear biomass clev con DD die egg gamma nu rec rep S X time prod caught

%% Benthic material
ncid = netcdf.open([fpath 'Data_fore_pristine2006-2085_bent.nc'],'NC_NOWRITE');
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

%WHOLE TIME 2006-2100
sp_mean=nanmean(SP.bio,2);
sf_mean=nanmean(SF.bio,2);
sd_mean=nanmean(SD.bio,2);
mp_mean=nanmean(MP.bio,2);
mf_mean=nanmean(MF.bio,2);
md_mean=nanmean(MD.bio,2);
lp_mean=nanmean(LP.bio,2);
ld_mean=nanmean(LD.bio,2);
b_mean=nanmean(BENT.bio,2);

sp_prod=nanmean(SP.prod,2);
sf_prod=nanmean(SF.prod,2);
sd_prod=nanmean(SD.prod,2);
mp_prod=nanmean(MP.prod,2);
mf_prod=nanmean(MF.prod,2);
md_prod=nanmean(MD.prod,2);
lp_prod=nanmean(LP.prod,2);
ld_prod=nanmean(LD.prod,2);

%% DECADAL
time=[1:size(SP.bio,2)]';
st=1:120:length(time);
en=120:120:length(time);
en(10)=time(end);

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
    
    sp_prod(:,n+1)=nanmean(SP.prod(:,st(n):en(n)),2);
    sf_prod(:,n+1)=nanmean(SF.prod(:,st(n):en(n)),2);
    sd_prod(:,n+1)=nanmean(SD.prod(:,st(n):en(n)),2);
    mp_prod(:,n+1)=nanmean(MP.prod(:,st(n):en(n)),2);
    mf_prod(:,n+1)=nanmean(MF.prod(:,st(n):en(n)),2);
    md_prod(:,n+1)=nanmean(MD.prod(:,st(n):en(n)),2);
    lp_prod(:,n+1)=nanmean(LP.prod(:,st(n):en(n)),2);
    ld_prod(:,n+1)=nanmean(LD.prod(:,st(n):en(n)),2);
    
end

%% last 50 yrs 2050-2100
mo=time/12;
mo=mo+2005;
yr50=find(mo>2050 & mo<=2100);

sp_mean5000=nanmean(SP.bio(:,yr50),2);
sf_mean5000=nanmean(SF.bio(:,yr50),2);
sd_mean5000=nanmean(SD.bio(:,yr50),2);
mp_mean5000=nanmean(MP.bio(:,yr50),2);
mf_mean5000=nanmean(MF.bio(:,yr50),2);
md_mean5000=nanmean(MD.bio(:,yr50),2);
lp_mean5000=nanmean(LP.bio(:,yr50),2);
ld_mean5000=nanmean(LD.bio(:,yr50),2);
b_mean5000=nanmean(BENT.bio(:,yr50),2);

sp_prod5000=nanmean(SP.prod(:,yr50),2);
sf_prod5000=nanmean(SF.prod(:,yr50),2);
sd_prod5000=nanmean(SD.prod(:,yr50),2);
mp_prod5000=nanmean(MP.prod(:,yr50),2);
mf_prod5000=nanmean(MF.prod(:,yr50),2);
md_prod5000=nanmean(MD.prod(:,yr50),2);
lp_prod5000=nanmean(LP.prod(:,yr50),2);
ld_prod5000=nanmean(LD.prod(:,yr50),2);


%%
save([fpath 'Means_fore_pristine2006-2085_' cfile '.mat'],...
    'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
    'lp_mean','ld_mean','sf_prod','sp_prod','sd_prod','mf_prod','mp_prod',...
    'md_prod','lp_prod','ld_prod','sf_mean5000','sp_mean5000','sd_mean5000',...
    'mf_mean5000',...
    'mp_mean5000','md_mean5000','b_mean5000','lp_mean5000','ld_mean5000',...
    'sf_prod5000','sp_prod5000','sd_prod5000','mf_prod5000','mp_prod5000',...
    'md_prod5000','lp_prod5000','ld_prod5000');

