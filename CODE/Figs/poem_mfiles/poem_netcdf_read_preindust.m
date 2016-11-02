% POEM output at all locations

clear all
close all

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05';

fpath=['/Volumes/GFDL/NC/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Data_preindust_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
SP.bio = biomass;
SP.prod = prod;
SP.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

%% SF
ncid = netcdf.open([fpath 'Data_preindust_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass;
SF.prod = prod;
SF.rec = rec;
SF.time = time;
SF.X = X;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% SD
ncid = netcdf.open([fpath 'Data_preindust_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
SD.prod = prod;
SD.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MP
ncid = netcdf.open([fpath 'Data_preindust_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
MP.prod = prod;
MP.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MF
ncid = netcdf.open([fpath 'Data_preindust_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
MF.prod = prod;
MF.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% MD
ncid = netcdf.open([fpath 'Data_preindust_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
MD.prod = prod;
MD.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LP
ncid = netcdf.open([fpath 'Data_preindust_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
LP.prod = prod;
LP.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% LD
ncid = netcdf.open([fpath 'Data_preindust_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
LD.prod = prod;
LD.rec = rec;

clear biomass clev con DD die egg gamma nu rec rep S X time prod

% Benthic material
ncid = netcdf.open([fpath 'Data_preindust_bent.nc'],'NC_NOWRITE');
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
%save([fpath 'Data_preindust_' cfile '.mat'])

%% 1800-1850
mo=(time-1)/12;
mo=mo+1760;
yr50=find(mo>=1800 & mo<1850);

sp_smean=nanmean(SP.bio(:,yr50),2);
sf_smean=nanmean(SF.bio(:,yr50),2);
sd_smean=nanmean(SD.bio(:,yr50),2);
mp_smean=nanmean(MP.bio(:,yr50),2);
mf_smean=nanmean(MF.bio(:,yr50),2);
md_smean=nanmean(MD.bio(:,yr50),2);
lp_smean=nanmean(LP.bio(:,yr50),2);
ld_smean=nanmean(LD.bio(:,yr50),2);
b_smean=nanmean(BENT.bio(:,yr50),2);

sp_prod=nanmean(SP.prod(:,yr50),2);
sf_prod=nanmean(SF.prod(:,yr50),2);
sd_prod=nanmean(SD.prod(:,yr50),2);
mp_prod=nanmean(MP.prod(:,yr50),2);
mf_prod=nanmean(MF.prod(:,yr50),2);
md_prod=nanmean(MD.prod(:,yr50),2);
lp_prod=nanmean(LP.prod(:,yr50),2);
ld_prod=nanmean(LD.prod(:,yr50),2);

sp_rec=nanmean(SP.rec(:,yr50),2);
sf_rec=nanmean(SF.rec(:,yr50),2);
sd_rec=nanmean(SD.rec(:,yr50),2);
mp_rec=nanmean(MP.rec(:,yr50),2);
mf_rec=nanmean(MF.rec(:,yr50),2);
md_rec=nanmean(MD.rec(:,yr50),2);
lp_rec=nanmean(LP.rec(:,yr50),2);
ld_rec=nanmean(LD.rec(:,yr50),2);

%%
save([fpath 'Means_preindust_' cfile '.mat'],...
    'sf_smean','sp_smean','sd_smean','mf_smean','mp_smean','md_smean','b_smean',...
    'lp_smean','ld_smean','sf_prod','sp_prod','sd_prod','mf_prod','mp_prod',...
    'md_prod','lp_prod','ld_prod','sf_rec','sp_rec','sd_rec','mf_rec',...
    'mp_rec','md_rec','lp_rec','ld_rec');

