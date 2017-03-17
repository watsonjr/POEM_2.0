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

%! Simname & directory
dirname
dpath = '/Volumes/GFDL/CSV/Matlab_new_size/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

%% ! Results
for L = 1:length(ids);
    ID = ids(L);
    loc = names{L};
    NX = length(ID);
    ENVR = sub_init_env(ID);
    %! Length of run (years)
    param.tEnd = 100;
    %! No Fishing
    param.F = param.dt * 0.0 *[0 0 0 0 0 0 0.1 1];
    t=0;
    for YR = 1:param.tEnd % years
        for DAY = 1:365
            t=t+1;
            ENVR = get_COBALT(COBALT,GRD,ID,DAY);
            if (YR==1 && DAY==1)
                result = poem_1D(param,ENVR,t);
            else
                result = poem_1D(param,ENVR,t,result);
            end
        end
    end
    result0 = result;
    %! Save
    save([dpath,simname,'/baserun_log_' loc '.mat'],'result')
    %! Plot
    plotPoem_1D(param, result)
    print('-dpng',[fpath,simname,'/baserun_log_' simname '_' loc])
    
    %! Fishing
    F = linspace(0.1, 1, 10);
    %Lenth of model run (years)
    param.tEnd = 50;
    Bf = NaN*ones(length(F),8);
    Y = NaN*ones(length(F),1);
    for i = 1%:length(F)
        %Fish D
        frate = F(i);
        param.F = param.dt * [0 0	0 0 0  0 0.1*F(i) F(i)];
        t=0;
        for YR = 1:param.tEnd % years
            for DAY = 1:365
                t=t+1;
                ENVR = get_COBALT(COBALT,GRD,ID,DAY);
                if (YR==1 && DAY==1)
                    result = poem_1D(param,ENVR,t,result0);
                else
                    result = poem_1D(param,ENVR,t,result);
                end
            end
        end
        lyr = (result.t(end)-364):result.t(end);
        Y(i) = sum(sum(result.Y(lyr,:)));
        %! Save
        dirname
        save([dpath,simname,'/fishing_log_' harv tfish(2:end) '_' loc '.mat'],'result')
        %! Plot
        plotPoem_1D(param, result)
        print('-dpng',[fpath,simname,'/fishing_log_' harv tfish(2:end) '_' loc])
    end
    %! Save
    dirname
    save([dpath,simname,'/fishing_log_' harv '_' loc '.mat'],'Y','F')
    
    %Plot yield vs. Fishing rate
    clf
    plot(F,Y,'bo-')
    xlabel('Fishing rate (yr^-^1)')
    ylabel([harv ' Yield'])
    title(['Repro effic = ' num2str(param.RE)])
    print('-dpng',[fpath,simname,'/fishing_log_' harv '_' simname '_' loc])
    
end


toc
