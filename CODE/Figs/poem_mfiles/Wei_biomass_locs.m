% Save Wei biomasses at locations

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

seafl = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei2010_Global_seafloor_biomass.csv',1,0);

%% Convert Wei units
Wcol = {'latitude','longitude','depth','bact.biom.mean','meio.biom.mean',...
    'macro.biom.mean','mega.biom.mean','inv.biom.mean','fis.biom.mean'};
Wcol = Wcol';
% all mean biomasses in log10 mg C/m2
invert = seafl(:,8);
fish = seafl(:,9);
% convert to g WW/m2
invert = 10.^(invert) * 1e-3 * 9.0;
fish = 10.^(fish) * 1e-3 * 9.0;

%% Plots in space
grid = csvread([cpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

Zi=griddata(seafl(:,2),seafl(:,1),invert,X,Y);
Zf=griddata(seafl(:,2),seafl(:,1),fish,X,Y);

glon = seafl(:,2);
glat = seafl(:,1);

%% Locations  
% Georges Bank (NEP)
lon=find(glon<=-66 & glon>=-67);
lat=find(glat<=42 & glat>=41);
gid=intersect(lon,lat);

% Eastern Bering Sea (Pribs)
lon=find(glon<=-169 & glon>=-170);
lat=find(glat<=57 & glat>=56);
eid=intersect(lon,lat);

% Subarctic Pacific Gyre (Ocean Station Papa)
lon=find(glon<=-145 & glon>=-146);
lat=find(glat<=51 & glat>=50);
pid=intersect(lon,lat);

% HOT
lon=find(glon<=-157 & glon>=-158);
lat=find(glat<=23 & glat>=22);
hid=intersect(lon,lat);

% BATS
lon=find(glon<=-64 & glon>=-65);
lat=find(glat<=32 & glat>=31);
bid=intersect(lon,lat);

% North Sea
lon=find(glon<=4 & glon>=3);
lat=find(glat<=57 & glat>=56);
nid=intersect(lon,lat);

% E Eq Pac
lon=find(glon<=-110 & glon>=-111);
lat=find(glat<=5.5 & glat>=5);
qid=intersect(lon,lat);

% Subpolar W Pac station K2: 47oN, 160oE
lon=find(glon<=161 & glon>=160);
lat=find(glat<=48 & glat>=47);
kid=intersect(lon,lat);

% Subtropical W Pac station S1: 30oN, 145oE
lon=find(glon<=146 & glon>=145);
lat=find(glat<=31 & glat>=30);
sid=intersect(lon,lat);

% Australia
lon=find(glon<=-125 & glon>=-126);
lat=find(glat<=-13.5 & glat>=-13.75);
aid=intersect(lon,lat);

% Peru Upwelling
lon=find(glon<=-79 & glon>=-80);
lat=find(glat<=13 & glat>=12);
uid=intersect(lon,lat);

names={'Georges Bank','Eastern Bering Sea','Ocean Station Papa','HOT',...
    'BATS','North Sea','Eastern Equatorial Pacific','K2','S1','Aus','Peru Upwell'};

ids(1,1)=gid;
ids(2,1)=eid;
ids(3,1)=pid;
ids(4,1)=hid;
ids(5,1)=bid;
ids(6,1)=nid;
ids(7,1)=qid;
ids(8,1)=kid;
ids(9,1)=sid;
ids(10,1)=aid;
ids(11,1)=uid;

B = invert(ids);
F = fish(ids);

%
T=table(names',B,F,...
    'VariableNames',{'Location','Inverts','Fish'});
writetable(T,'/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei_inverts_fish_gWWm2_locs.csv','Delimiter',',');
save('/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei_inverts_fish_gWWm2_locs.mat','T');