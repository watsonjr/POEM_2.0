%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';
datap = '/Volumes/GFDL/CSV/';

npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath6 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath7 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath8 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath9 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath10 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath11 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath12 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath13 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath14 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath15 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath16 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
% npath17 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
% npath18 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
% npath19 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE025/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE075/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20/';
npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25/';
% npath12 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
% npath13 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/';
% npath14 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath15 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20/';
% npath16 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30/';
% npath18 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30/';

% dp = {npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
dp = {npath7};
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;...
% npath10;npath11;npath12;npath13;npath14;npath15;npath16;npath17;npath18;npath19;npath20};

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

load('cmap_ppt_angles.mat')

%%
for i=1:length(dp)
    close all
    
    dpath = [datap char(dp(i))];
    fpath = [figp char(dp(i))];
    cfile = char(dp(i));
    
    load([dpath sname sname2 'lastyr_sum_mean_biom']);
    
    %% Sum mean biom over stages    
    sumspec = squeeze(nansum(nansum(all_mean)));
    [ss,ix] = sort(sumspec);
    sspot = spots(ix);
    
    figure(17);
    plot(1:9,log10(ss),'k.','MarkerSize',25); hold on;
    xlim([0 10])
    %ylim([-2 1])
    set(gca,'XTick',1:9,'XTickLabel',[])
    for n=1:9
        text(n,-1.45,sspot{n},'HorizontalAlignment','center')
    end
    title('log10 mean biom of All Fishes (g m^-^2)')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_spec_horz.png'])
    
    %%
    figure(18);
    plot(1:9,log10(ss),'k.','MarkerSize',25); hold on;
    xlim([0 10])
    ylim([-2 2])
    set(gca,'XTick',1:9,'XTickLabel',[])
    for n=1:9
        text(n,-2.1,sspot{n},'HorizontalAlignment','center')
    end
    title('log10 mean biom of All Fishes (g m^-^2)')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_spec_horz_axes.png'])
    
    %%
    ss2 = flipud(ss);
    sspot2 = fliplr(sspot);
    figure(19);
    plot(log10(ss),9:-1:1,'k.','MarkerSize',25); hold on;
    %xlim([0 10])
    set(gca,'YTick',1:9,'YTickLabel',[])
    for n=1:9
        text(-1.41,n,sspot2{n},'HorizontalAlignment','right')
    end
    title('log10 mean biom of All Fishes (g m^-^2)')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_spec_vert.png'])
    
    %%
    figure(20);
    plot(log10(ss),9:-1:1,'k.','MarkerSize',25); hold on;
    xlim([-1.5 2])
    set(gca,'YTick',1:9,'YTickLabel',[])
    for n=1:9
        text(-1.55,n,sspot2{n},'HorizontalAlignment','right')
    end
    title('log10 mean biom of All Fishes (g m^-^2)')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_spec_vert_axes.png'])
    
    
end




