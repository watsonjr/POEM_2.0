%
% Sweep over Search area coefficient
%
clear all
close all

tic

%! Set parameters
baseparams_1D_lowIC_1B

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
dirname_1B
simname0 = simname;
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
    t=0;
    for YR = 1:param.tEnd % years
        for DAY = 1:365
            t=t+1;
            ENVR = get_COBALT(COBALT,GRD,ID,DAY);
            if (YR==1 && DAY==1)
                result = poem_1D_jhp2_1B(param,ENVR,t);
            else
                result = poem_1D_jhp2_1B(param,ENVR,t,result);
            end
        end
    end
    %! Save
    save([dpath,simname,'/baserun_1B_lowIC_jhp2_log_' loc '.mat'],'result')
    %! Plot
    plotPoem_1D_1B(param, result)
    print('-dpng',[fpath,simname,'/baserun_1B_lowIC_jhp2_log_' loc])
    result0.y = result.y(end,:);
    clear result
    
    %! Search
    C = 20:10:90;
    C = sort(C,'descend');
    %Lenth of model run (years)
    param.tEnd = 50;
    B = NaN*ones(length(C),11);
    fl = NaN*ones(length(C),11);
    for i = 1:length(C)
        param.gamma = C(i);
        param.V = param.dt * param.gamma * param.wc .^ param.q;
        t=0;
        for YR = 1:param.tEnd % years
            for DAY = 1:365
                t=t+1;
                ENVR = get_COBALT(COBALT,GRD,ID,DAY);
                if (YR==1 && DAY==1)
                    result = poem_1D_jhp2_1B(param,ENVR,t,result0);
                else
                    result = poem_1D_jhp2_1B(param,ENVR,t,result);
                end
            end
        end
        lyr = (result.t(end)-364):result.t(end);
        B(i,:) = mean(result.y(lyr,:));
        fl(i,:) = mean(result.f(lyr,:));
        %! Save
        dirname_1B
        save([dpath,simname,'/baserun_1B_lowIC_jhp2_log_' loc '.mat'],'result')
        %! Plot
        plotPoem_1D_1B(param, result)
        print('-dpng',[fpath,simname,'/baserun_1B_lowIC_jhp2_log_' loc])
        clear result
    end
    %! Save
    save([dpath,simname0,'/Search_1B_lowIC_jhp2_log_' loc '.mat'],'B','fl','C')
    
    %Plot adult biomasses
    ix = cell([3,1]);
    FF = 1;
    LP = 2;
    LD = 3;
    
    ix{FF} = 4:5;
    ix{LP} = 6:8;
    ix{LD} = 9:11;
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
    xlabel('Search coeff')
    ylabel('log10 Biomass')
    print('-dpng',[fpath,'Search_biom_1B_lowIC_jhp2_log_' simname0 '_' loc])
    
    %Plot Feeding level
    clf
    subplot(3,1,1)
    bar(fl(:,[4,6,9])); hold on;
    colormap([1 0 0; 0 0 1; 0 0 0]);
    plot(0:11, param.fc*ones(12,1), 'k--')
    ylim([0 1])
    xlim([0 11])
    ylabel('Feeding level')
    xlabel('Search coeff')
    title('S')
    set(gca,'XTickLabel',C)
    
    subplot(3,1,2)
    bar(fl(:,[5,7,10])); hold on;
    colormap([1 0 0; 0 0 1; 0 0 0]);
    plot(0:11, param.fc*ones(12,1), 'k--')
    ylim([0 1])
    xlim([0 11])
    ylabel('Feeding level')
    xlabel('Search coeff')
    title('M')
    set(gca,'XTickLabel',C)
    
    subplot(3,1,3)
    bar(fl(:,[1,8,11])); hold on;
    %colormap([0 0 1; 0 0 0]);
    plot(0:11, param.fc*ones(12,1), 'k--')
    ylim([0 1])
    xlim([0 11])
    ylabel('Feeding level')
    xlabel('Search coeff')
    title('L')
    set(gca,'XTickLabel',C)
    print('-dpng',[fpath,'Search_flev_1B_lowIC_jhp2_' simname0 '_' loc])
end
%%


toc

