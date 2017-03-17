%
% run all 3 types at one COBALT grid cell
%
clear all
close all

tic

%! Set parameters
baseparams_1D

%! Setup spinup (loop year of COBALT)
load('/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2000.mat');

%! Where to run the model
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_grid_hindcast_NOTflipped.mat');
ids = [40319,42639,41782,36334,38309,42744,30051,41284,38003,19327,20045];
names = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};

%! Fishing rate
param.F = param.dt * 0.0 *[0 0 0 0 0 0 0.1 1];

%! Length of run (years)
param.tEnd = 100;

%! Simname & directory
dirname

%% ! Results
for L = 1:length(ids);
    ID = ids(L);
    loc = names{L};
    NX = length(ID);
    ENVR = sub_init_env(ID);
    t=0;
    for YR = 1:param.tEnd % years
        for DAY = 1:365
            t=t+1;
            ENVR = get_COBALT(COBALT,GRD,ID,DAY);
            if (YR==1 && DAY==1)
                result = poem_1D_jhp2(param,ENVR,t);
            else
                result = poem_1D_jhp2(param,ENVR,t,result);
            end
        end
    end
    save(['/Volumes/GFDL/CSV/Matlab_new_size/',simname,'/baserun_jhp2_log_' loc '.mat'],'result')
    %! Plot
    plotPoem_1D(param, result)
    print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/baserun_jhp2_log_' simname '_' loc])
    
end

toc

