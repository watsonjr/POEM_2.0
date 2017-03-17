clear all
close all

%% COBALT file
%! Setup spinup (loop year of COBALT)
load('/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2000.mat');

%! Where to run the model
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_grid_hindcast_NOTflipped.mat');
ids = [40319,42639,41782,36334,38309,42744,30051,41284,38003,19327,20045];
names = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};

%! Simname & directory
dpath = '/Volumes/GFDL/CSV/Matlab_new_size/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

ENVR = sub_init_env(ids);
DY=1:365;
ENVR = get_COBALT(COBALT,GRD,ids,DY);

%% Run Ken's model

%Set parameters
baseparameters

%Fishing rate
param.F = 0.0*[0 0 0 0 0 0 0.1 1]';
%Length of run (years)
param.tEnd = 200;
%Results
result = poem(param);
%Plot
plotPoem(param, result)