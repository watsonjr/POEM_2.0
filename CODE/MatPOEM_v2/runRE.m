%
% Run over RE values
%
clear all
close all

%Set parameters
baseparameters

%Simname & directory
dirname
simname0 = simname;

factor = logspace(-5,0,20);
%%
B = NaN*ones(length(factor),13);
F = linspace(0, 1, 10);
YF = NaN*ones(length(F),length(factor));
YP = NaN*ones(length(F),length(factor));
YD = NaN*ones(length(F),length(factor));
for r = 1:length(factor)
    %Fishing rate
    param.F = 0.0*[0 0 0 0 0 0 0.1 1]';
    %Length of run (years)
    param.tEnd = 200;
    %Repro effic
    param.RE = factor(r);
    param.eRepro = param.RE*[1 1 1];
    result = poem(param);
    result0 = result;
    base(r) = result;
    B(r,:) = result.y(end,:);
    
    dirname
    %Plot
    clf
    plotPoem(param, result)
    print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/baserun_' simname])
    
    %Fishing
    result = result0;
    runfishingF
    YF(:,r) = Y;
    result = result0;
    runfishingP
    YP(:,r) = Y;
    result = result0;
    runfishingD
    YD(:,r) = Y;
end

%Save
save(['/Volumes/GFDL/CSV/Matlab_new_size/',simname0,'/RErun.mat'],'base','B','factor',...
    'YF','YP','YD');

%% Plot all yield curves together
cm21=[1 0.5 0;...   %orange
    0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    0 1 1;...     %c
    0 0 0.75;...  %b
    0.5 0 1;...   %purple
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0.75 0.75 0.75;... %lt grey
    0.5 0.5 0.5;...    %med grey
    49/255 79/255 79/255;... %dk grey
    0 0 0;...      %black
    1 1 0;...      %yellow
    127/255 255/255 0;... %lime green
    0 0.5 0;...    %dk green
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    188/255 143/255 143/255;... %rosy brown
    255/255 192/255 203/255;... %pink
    255/255 160/255 122/255]; %peach

set(groot,'defaultAxesColorOrder',cm21);

clf
plot(F,YF,'-')
xlabel('Fishing rate (yr^-^1)','LineWidth',1.5)
ylabel('F Yield')
legend(num2str(factor'))
legend('location','eastoutside')
print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/fishing_F_allRE_' simname0])

figure
plot(F,YP,'-')
xlabel('Fishing rate (yr^-^1)','LineWidth',1.5)
ylabel('P Yield')
legend(num2str(factor'))
legend('location','eastoutside')
print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/fishing_P_allRE_' simname0])

figure
plot(F,YD,'-')
xlabel('Fishing rate (yr^-^1)','LineWidth',1.5)
ylabel('D Yield')
legend(num2str(factor'))
legend('location','eastoutside')
print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/fishing_D_allRE_' simname0])

%% Plot biomass without fishing
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
%subplot(2,2,1)
for i = 1:3
    plot(log10(factor), log10(B(:,ix{i}(end))), col{i},'LineWidth',2)
    hold on
end
legend('MF','LP','LD')
xlabel('log10 RE')
ylabel('log10 Biomass')
print('-dpng',['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/baserun_allRE_' simname0])


