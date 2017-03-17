%
% Sweep over Cmax coefficient
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

%! Simname & directory
dirname
simname0 = simname;
dpath = '/Volumes/GFDL/CSV/Matlab_new_size/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

%% ! Results
for L = 9%1:length(ids);
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
    %! Save
    save([dpath,simname,'/baserun_log_' loc '.mat'],'result')
    %! Plot
    plotPoem_1D(param, result)
    print('-dpng',[fpath,simname,'/baserun_log_' loc])
    result0.y = result.y(end,:);
    clear result
    
    %! Cmax
    C = linspace(10, 100, 10);
    %Lenth of model run (years)
    param.tEnd = 50;
    B = NaN*ones(length(C),13);
    fl = NaN*ones(length(C),13);
    for i = 1:length(C)
        param.h = C(i);
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
        B(i,:) = mean(result.y(lyr,:));
        fl(i,:) = mean(result.f(lyr,:));
        %! Save
        dirname
        save([dpath,simname,'/baserun_log_' loc '.mat'],'result')
        %! Plot
        plotPoem_1D(param, result)
        print('-dpng',[fpath,simname,'/baserun_log_' loc])
        clear result
    end
    %! Save
    save([dpath,simname0,'/Cmax_log_' loc '.mat'],'B','fl','C')
    
    %Plot adult biomasses
    ix = cell([3,1]);
    FF = 1;
    LP = 2;
    LD = 3;
    
    ix{FF} = 6:7;
    ix{LP} = 8:10;
    ix{LD} = 11:13;
    col = cell([3,1]);
    col{FF} = 'r-';
    col{LP} = 'b-';
    col{LD} = 'k-';
    
    clf
    for i = 1:3
        plot(C, log10(B(:,ix{i}(end))), col{i},'LineWidth',2)
        hold on
    end
    legend('MF','LP','LD')
    xlabel('Cmax coeff')
    ylabel('log10 Biomass')
    print('-dpng',[fpath,'Cmax_log_' simname0 '_' loc])
    
    %Plot Feeding level
    clf
    subplot(3,1,1)
    bar(fl(:,[6,8,11])); hold on;
    colormap([1 0 0; 0 0 1; 0 0 0]);
    plot(0:11, param.fc*ones(12,1), 'k--')
    ylim([0 1])
    xlim([0 11])
    ylabel('Feeding level')
    xlabel('Cmax coeff')
    title('S')
    set(gca,'XTickLabel',C)
    
    subplot(3,1,2)
    bar(fl(:,[7,9,12])); hold on;
    colormap([1 0 0; 0 0 1; 0 0 0]);
    plot(0:11, param.fc*ones(12,1), 'k--')
    ylim([0 1])
    xlim([0 11])
    ylabel('Feeding level')
    xlabel('Cmax coeff')
    title('M')
    set(gca,'XTickLabel',C)
    
    subplot(3,1,3)
    bar(fl(:,[1,10,13])); hold on;
    %colormap([0 0 1; 0 0 0]);
    plot(0:11, param.fc*ones(12,1), 'k--')
    ylim([0 1])
    xlim([0 11])
    ylabel('Feeding level')
    xlabel('Cmax coeff')
    title('L')
    set(gca,'XTickLabel',C)
    print('-dpng',[fpath,'Cmax_flev_log_' simname0 '_' loc])
end
%%




