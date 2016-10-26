%POEM catch vs. SAUP catch by LME
%Assume each location represents its LME

clear all
close all

datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

load([spath 'LME_Catch_annual.mat']);
load([cpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([cpath 'lme_mask_esm2m.mat']);
grid = csvread([cpath 'grid_csv.csv']);

npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish05/'];
npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish10/'];
npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish15/'];
npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish20/'];
npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish30/'];
npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish40/'];
npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish50/'];
npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fishM43L12/'];
npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish60/'];
npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish70/'];
npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish80/'];
npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish90/'];
npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish100/'];
npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish110/'];
npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish120/'];
npath16 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish130/'];
npath17 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish140/'];
npath18 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish150/'];
npath19 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish160/'];
npath20 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish170/'];
npath21 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish180/'];
npath22 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish190/'];

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

dp = {npath2;npath4;npath5;npath6;npath7;npath9;npath10;npath11;...
    npath12;npath13;npath14;npath15;npath16;npath17;npath18;npath19;npath20;...
    npath21;npath22};
sims = [0.1:0.1:1.9];
cfile = 'Dc_MFeqMP_MZ01_BE05_fishing_catch';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

%%

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met'};
cols=cols';

load('cmap_ppt_angles.mat')

%% Total catch
Tot = NaN*ones(length(dp),length(spots));
for d=1:length(dp)
    
    dpath = char(dp(d));
    load([dpath sname sname2 'lastyr_sum_mean_biom']);
    
    Ftot = sum(Ftotcatch);
    Ptot = sum(Ptotcatch);
    Dtot = sum(Dtotcatch);
    Tot(d,:) = Ftot+Ptot+Dtot;
end

%% Which LME
Lat = geolat_t;
Lon = geolon_t;
%fix lon shift
id=find(Lon(:)<-180);
Lon(id)=Lon(id)+360;

% Georges Bank (NEP)
lon=find(Lon<=-66 & Lon>=-67);
lat=find(Lat<=42 & Lat>=41);
gid=intersect(lon,lat);

% Eastern Bering Sea (Pribs)
lon=find(Lon<=-169 & Lon>=-170);
lat=find(Lat<=57 & Lat>=56);
eid=intersect(lon,lat);

% Subarctic Pacific Gyre (Ocean Station Papa)
lon=find(Lon<=-145 & Lon>=-146);
lat=find(Lat<=51 & Lat>=50);
pid=intersect(lon,lat);

% HOT
lon=find(Lon<=-157 & Lon>=-158);
lat=find(Lat<=23 & Lat>=22);
hid=intersect(lon,lat);

% BATS
lon=find(Lon<=-64 & Lon>=-65);
lat=find(Lat<=32 & Lat>=31);
bid=intersect(lon,lat);

% North Sea
lon=find(Lon<=4 & Lon>=3);
lat=find(Lat<=57 & Lat>=56);
nid=intersect(lon,lat);

% Eastern Equatorial Pacific
lon=find(Lon<=-110 & Lon>=-111);
lat=find(Lat<=5.5 & Lat>=5);
qid=intersect(lon,lat);

%Subpolar W Pac station K2: 47oN, 160oE
lon=find(Lon<=161 & Lon>=160);
lat=find(Lat<=48 & Lat>=47);
kid=intersect(lon,lat);

%Subtropical W Pac station S1: 30oN, 145oE
lon=find(Lon<=146 & Lon>=145);
lat=find(Lat<=31 & Lat>=30);
sid=intersect(lon,lat);

tlme = lme_mask_esm2m';
lme_ids(1,1)=tlme(gid);
lme_ids(2,1)=tlme(eid);
lme_ids(3,1)=tlme(pid);
lme_ids(4,1)=tlme(hid);
lme_ids(5,1)=tlme(bid);
lme_ids(6,1)=tlme(nid);
lme_ids(7,1)=tlme(qid);
lme_ids(8,1)=tlme(kid);
lme_ids(9,1)=tlme(sid);

lme_ids(4,1)=10; %HOT
lme_ids(8,1)=51; %K2
lme_ids(9,1)=49; %S1

%% Areal totals
lmes = find(~isnan(lme_ids));
lme_area = NaN*ones(size(lme_ids));
loc_catch = NaN*ones(size(lme_ids));
for n=1:length(lmes)
    L = lme_ids(lmes(n));
    r = find(lme_ids==L);
    lid = find(tlme==L);
    lme_area(r,1) = nansum(dat(lid));
    
    loc_catch(r,1) = lme_catch(31,L);
end

lme_Tot = Tot';
%LME biomass in MT
lme_aTot = lme_Tot .* repmat(lme_area,1,length(dp)) * 1e-6;

%Difference
diff_aTot = lme_aTot - repmat(loc_catch,1,length(dp));
sum_diff = nansum(diff_aTot);



