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

%% Compare
lyr = find(result.t<= 200 & result.t> 199);
R = mean(result.R(lyr,:));
Zm = mean(ENVR.Zm,2);
Zl = mean(ENVR.Zl,2);
dZm = mean(ENVR.dZm,2);
dZl = mean(ENVR.dZl,2);

%%
result0=result;
clear result
load([dpath,simname,'/baserun_log_S1.mat'],'result')
resultS=result;
clear result
load([dpath,simname,'/baserun_log_GB.mat'],'result')
resultG=result;
clear result
lyr = (resultS.t(end)-364):resultS.t(end);
RS = mean(resultS.R(lyr,:));
lyr = (resultG.t(end)-364):resultG.t(end);
RG = mean(resultG.R(lyr,:));

Rs = R./RS;
Rg = R./RG;

figure
subplot(2,2,1)
bar([R;RS;RG]')

figure
subplot(2,2,1)
bar(log10([Rs;Rg]'))

%%
figure
subplot(2,2,1)
plot(DY,ENVR.Zm,'LineWidth',1); hold on
plot(DY,ENVR.dZm,'--','LineWidth',1);

figure
subplot(2,2,1)
plot(DY,ENVR.Zl,'LineWidth',1); hold on
plot(DY,ENVR.dZl,'--','LineWidth',1);


