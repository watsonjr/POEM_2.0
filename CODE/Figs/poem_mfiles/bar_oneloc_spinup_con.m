%Visualize consumption output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';

dp = {npath1};

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

load('cmap_ppt_angles.mat')

cm={[1 0.5 0],...   %orange
    [0.5 0.5 0],... %tan/army
    [0 0.7 0],...   %g
    [0 1 1],...     %c
    [0 0 0.75],...  %b
    [0.5 0 1],...   %purple
    [1 0 1],...     %m
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.75 0.75 0.75],... %lt grey
    [0.5 0.5 0.5],...    %med grey
    [49/255 79/255 79/255],... %dk grey
    [0 0 0],...      %black
    [1 1 0],...      %yellow
    [127/255 255/255 0],... %lime green
    [0 0.5 0],...    %dk green
    [0/255 206/255 209/255],... %turq
    [0 0.5 0.75],...   %med blue
    [188/255 143/255 143/255],... %rosy brown
    [255/255 192/255 203/255],... %pink
    [255/255 160/255 122/255]}; %peach

L_s = 10^((log10(2)+log10(20))/2);
L_m = 10^((log10(20)+log10(200))/2);
L_l = 10^((log10(200)+log10(2000))/2);

M_s = 0.01 * (0.1*L_s)^3;
M_m = 0.01 * (0.1*L_m)^3;
M_l = 0.01 * (0.1*L_l)^3;

mass = [M_s;M_m;M_l];
mass = repmat(mass,1,length(spots));

%%
for i=1:length(dp)
    
    dpath = [datap char(dp(i))];
    fpath = [figp char(dp(i))];
    cfile = char(dp(i));
    
    load([dpath sname sname2 'lastyr_sum_mean_biom.mat'])
    
    %% Consump g/g/d --> g/d --> g/y
    Pcon = Pcon .* mass .* 365;
    Fcon = Fcon .* mass(1:2,:) .* 365;
    Dcon = Dcon .* mass .* 365;
    
    f4 = figure(4);
    subplot(3,1,1)
    plot(0.75:1:8.75,Fcon(1,:),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:1:9,Pcon(1,:),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.25:1:9.25,Dcon(1,:),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 10])
    ylabel('S')
    set(gca,'XTick',1:10,'XTickLabel',spots)
    title('Mean consumption (g y^-^1) in final year')
    
    subplot(3,1,2)
    plot(0.75:1:8.75,Fcon(2,:),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:1:9,Pcon(2,:),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.25:1:9.25,Dcon(2,:),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 10])
    ylabel('M')
    set(gca,'XTick',1:10,'XTickLabel',spots)
    
    subplot(3,1,3)
    plot(1:1:9,Pcon(3,:),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.25:1:9.25,Dcon(3,:),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 10])
    ylabel('L')
    set(gca,'XTick',1:10,'XTickLabel',spots)
    xlabel('Location')
    stamp(cfile)
    
    print(f4,'-dpng',[fpath sname sname2 'All_oneloc_consump_gyr.png'])
    
    %% Consump vs. weight
    A = 4.39;
    fc = 0.2;
    f0 = 0.6;
    epsassim = 0.6;
    n = 3/4;
    
    w = logspace(-3, 5);
    
    AvailEnergy = A*w.^n;
    Consumption = A / (epsassim*(f0-fc)) * w.^n;
    
    for s=1:length(spots)
        figure(2)
        subplot(1,3,1)
        loglog((mass(1:2,1)),(Fcon(:,s)),'.',...
            'Color',cm{s},'MarkerSize',25); hold on;
        title('F')
        xlabel('Mass (g)')
        ylabel('Mean consumption (g y^-^1) in final year')
        legend(spots)
        legend('location','northeast')
        %axis([-5 5 -1 5])
        
        subplot(1,3,2)
        loglog((mass(:,1)),(Pcon(:,s)),'.',...
            'Color',cm{s},'MarkerSize',25); hold on;
        title('P')
        xlabel('Mass (g)')
        %axis([-5 5 -1 5])
        
        subplot(1,3,3)
        loglog((mass(:,1)),(Dcon(:,s)),'.',...
            'Color',cm{s},'MarkerSize',25); hold on;
        title('D')
        xlabel('Mass (g)')
        %axis([-5 5 -1 5])
        stamp(cfile)
        print('-dpng',[fpath sname sname2 'All_oneloc_consump_gyr_vs_weight_v2.png'])
    end
    subplot(1,3,1)
    loglog(w, Consumption,'k')
    
    subplot(1,3,2)
    loglog(w, Consumption,'k')
    
    subplot(1,3,3)
    loglog(w, Consumption,'k')
    print('-dpng',[fpath sname sname2 'All_oneloc_consump_gyr_vs_weight_compare.png'])
    
end