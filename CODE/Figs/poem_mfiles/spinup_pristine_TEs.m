%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
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
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE025/';
npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/';
npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE075/';
npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10/';
npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20/';

%cfile = 'Dc_MFeqMP_fcrit10_MZ01_NOnmort';
dp = {npath1;npath2;npath3;npath4;npath5};%;npath6;npath7;npath8;npath9};
% dp = {npath6;npath7;npath8;npath9;npath10};
% dp = {npath11;npath12;npath13;npath14;npath15};
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;...
%     npath11;npath12;npath13;npath14;npath15;npath16;npath17;npath18;npath19;...
%     npath20;npath21};

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
    
    dpath = [datap char(dp(i))];
    fpath = [figp char(dp(i))];
    cfile = char(dp(i));
    
    close all
    
    %% Zoop con
    load([dpath sname sname2 'consump.mat'],'mclev','Zcon');
    
    %%
    Con=NaN*ones(3,length(spots));
    Scon = Con;
    Prod = Con;
    Pred = Con;
    Eat = Con;
    PNmort = Con;
    PNFmort = Con;
    TEcon = NaN*ones(2,length(spots));
    TEscon = TEcon;
    TEprod = TEcon;
    TEpred = TEcon;
    TEeat = TEcon;
    TEmortPN = TEcon;
    TEmortPNF = TEcon;
    
    %%
    for s=1:length(spots)
        %%
        loc = spots{s};
        lname = [sname2 loc '_'];
        SP = csvread([dpath sname lname 'Sml_p.csv']);
        SF = csvread([dpath sname lname 'Sml_f.csv']);
        SD = csvread([dpath sname lname 'Sml_d.csv']);
        MP = csvread([dpath sname lname 'Med_p.csv']);
        MF = csvread([dpath sname lname 'Med_f.csv']);
        MD = csvread([dpath sname lname 'Med_d.csv']);
        LP = csvread([dpath sname lname 'Lrg_p.csv']);
        LD = csvread([dpath sname lname 'Lrg_d.csv']);
        C = csvread([dpath sname lname 'Cobalt.csv']);
        
        t=1:length(SP);
        lyr=t((end-365+1):end);
        
        %% Consumption efficiency
        % Biomass consumed (g/day) (I * biom)
        SP_con=(SP(lyr,14).*SP(lyr-1,1));
        SF_con=(SF(lyr,14).*SF(lyr-1,1));
        SD_con=(SD(lyr,14).*SD(lyr-1,1));
        MP_con=(MP(lyr,14).*MP(lyr-1,1));
        MF_con=(MF(lyr,14).*MF(lyr-1,1));
        MD_con=(MD(lyr,14).*MD(lyr-1,1));
        LP_con=(LP(lyr,14).*LP(lyr-1,1));
        LD_con=(LD(lyr,14).*LD(lyr-1,1));
        
        S_con=nansum(SP_con+SF_con+SD_con);
        M_con=nansum(MP_con+MF_con+MD_con);
        L_con=nansum(LP_con+LD_con);
        
        Con(:,s) = [S_con;M_con;L_con];
        
        TEcon(1,s) = nanmean(M_con./S_con);
        TEcon(2,s) = nanmean(L_con./M_con);
        
        %% Production (= nu * biom) efficiency
        SP_prod=(SP(lyr,24));
        SF_prod=(SF(lyr,24));
        SD_prod=(SD(lyr,24));
        MP_prod=(MP(lyr,24));
        MF_prod=(MF(lyr,24));
        MD_prod=(MD(lyr,24));
        LP_prod=(LP(lyr,24));
        LD_prod=(LD(lyr,24));
        
        S_prod=nansum(SP_prod+SF_prod+SD_prod);
        M_prod=nansum(MP_prod+MF_prod+MD_prod);
        L_prod=nansum(LP_prod+LD_prod);
        
        Prod(:,s) = [S_prod;M_prod;L_prod];
        
        TEprod(1,s) = M_prod./S_prod;
        TEprod(2,s) = L_prod./M_prod;
        
        %% Predation efficiency in biomass eaten (g/m2/d)
        SP_die=(SP(lyr,17));
        SF_die=(SF(lyr,17));
        SD_die=(SD(lyr,17));
        MP_die=(MP(lyr,17));
        MF_die=(MF(lyr,17));
        MD_die=(MD(lyr,17));
        LP_die=(LP(lyr,17));
        LD_die=(LD(lyr,17));
        
        S_die=nansum(SP_die+SF_die+SD_die);
        M_die=nansum(MP_die+MF_die+MD_die);
        L_die=nansum(LP_die+LD_die);
        
        Eat(:,s) = [S_die;M_die;L_die];
        
        TEeat(1,s) = M_die./S_die;
        TEeat(2,s) = L_die./M_die;
        
        %% Natural mortality biomass loss (g/m2/day)
        SP_nat=(SP(lyr,26).*SP(lyr-1,1));
        SF_nat=(SF(lyr,26).*SF(lyr-1,1));
        SD_nat=(SD(lyr,26).*SD(lyr-1,1));
        MP_nat=(MP(lyr,26).*MP(lyr-1,1));
        MF_nat=(MF(lyr,26).*MF(lyr-1,1));
        MD_nat=(MD(lyr,26).*MD(lyr-1,1));
        LP_nat=(LP(lyr,26).*LP(lyr-1,1));
        LD_nat=(LD(lyr,26).*LD(lyr-1,1));
        
        S_nat=nansum(SP_nat+SF_nat+SD_nat);
        M_nat=nansum(MP_nat+MF_nat+MD_nat);
        L_nat=nansum(LP_nat+LD_nat);
        
        %% Fishing catch (g/m2/day)
        MP_fish=(MP(lyr,28));
        MF_fish=(MF(lyr,28));
        MD_fish=(MD(lyr,28));
        LP_fish=(LP(lyr,28));
        LD_fish=(LD(lyr,28));
        
        S_fish=0;
        M_fish=nansum(MF_fish+MP_fish+MD_fish);
        L_fish=nansum(LP_fish+LD_fish);
        
        %% Total mortality w/o fishing efficiency
        Smort = S_die + S_nat;
        Mmort = M_die + M_nat;
        Lmort = L_die + L_nat;
        
        PNmort(:,s) = [Smort;Mmort;Lmort];
        
        TEmortPN(1,s) = Mmort./Smort;
        TEmortPN(2,s) = Lmort./Mmort;
        
        %% Total mortality w/ fishing efficiency
        Smortf = S_die + S_nat + S_fish;
        Mmortf = M_die + M_nat + M_fish;
        Lmortf = L_die + L_nat + L_fish;
        
        PNFmort(:,s) = [Smortf;Mmortf;Lmortf];
        
        TEmortPNF(1,s) = Mmortf./Smortf;
        TEmortPNF(2,s) = Lmortf./Mmortf;
        
    end
    
    save([dpath sname sname2 'lastyr_TEs.mat'],'TEcon','TEprod','TEeat','TEmortPN','TEmortPNF');
    
    %% Figure
    figure(1);
    subplot(2,1,1)
    plot(1:length(spots),TEcon(1,:),'.k','MarkerSize',25); hold on;
    xlim([0 length(spots)+1])
    set(gca,'XTick',1:length(spots),'XTickLabel',spots)
    ylabel('mean TE in final year')
    xlabel('Site')
    title('Medium/Small')
    stamp(cfile)
    
    subplot(2,1,2)
    plot(1:length(spots),TEcon(1,:),'.k','MarkerSize',25); hold on;
    xlim([0 length(spots)+1])
    set(gca,'XTick',1:length(spots),'XTickLabel',spots)
    ylabel('mean TE in final year')
    xlabel('Site')
    print('-dpng',[fpath sname sname2 'TEs_M.png'])
    
    figure(2);
    subplot(2,1,1)
    plot(1:length(spots),TEcon(2,:),'.k','MarkerSize',25); hold on;
    xlim([0 length(spots)+1])
    set(gca,'XTick',1:length(spots),'XTickLabel',spots)
    ylabel('mean L TE in final year')
    xlabel('Site')
    title('Large/Medium')
    stamp(cfile)
    
    subplot(2,1,2)
    plot(1:length(spots),TEcon(2,:),'.k','MarkerSize',25); hold on;
    xlim([0 length(spots)+1])
    ylim([0 1])
    set(gca,'XTick',1:length(spots),'XTickLabel',spots)
    ylabel('mean L TE in final year')
    xlabel('Site')
    print('-dpng',[fpath sname sname2 'TEs_L.png'])
    
end




