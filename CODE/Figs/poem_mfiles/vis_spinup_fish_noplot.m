% Visualize output of POEM
% Pristine historical at one location
% 145 years

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish01/';
npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish02/';
npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish03/';
npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish04/';
npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish05/';
npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish06/';
npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish01/';
npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish02/';
npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish03/';
npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish04/';
npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish05/';
npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish06/';
npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish01/';
npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish02/';
npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish03/';
npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish04/';
npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish05/';
npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish06/';
npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish01/';
npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish02/';
npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish03/';
npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish04/';
npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish05/';
npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish06/';
npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish01/';
npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish02/';
npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish03/';
npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish04/';
npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish05/';
npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish06/';
npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish01/';
npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish02/';
npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish03/';
npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish04/';
npath35 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish05/';
npath36 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish06/';
npath37 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish01/';
npath38 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish02/';
npath39 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish03/';
npath40 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish04/';
npath41 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish05/';
npath42 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish06/';


% dp = {npath3};
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% sims = {'0.0MZ01','0.025MZ01','0.05MZ01','0.10MZ01','0.20MZ01','0.30MZ01',...
%     '0.40MZ01','0.50MZ01','0.60MZ01','0.70MZ01','0.45MZ01','0.55MZ01'};
% dp = {npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9};
% sims = {'0.05','0.10','0.20','0.30','0.40','0.50','0.60','0.70'};
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;...
%     npath9};
% sims = {'0.025','0.05','0.075','0.10','0.20','0.30','0.40','0.50',...
%     '0.60','0.70'};
% dp = {npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% dp={npath10};
dp = {npath37;npath38;npath39;npath40;npath41;npath42};
% dp = {npath9;npath10;npath11;npath12;npath13;npath14;npath15;npath16};

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};

cols = {'clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

fplot=0;

%%
for i=4:length(dp)
    
    dpath = [datap char(dp(i))];
    fpath = [figp char(dp(i))];
    
    mclev=NaN*ones(length(spots),8);
    Zcon=NaN*ones(length(spots),3);
    
    %%
    for s=1:length(spots)
        %%
        loc = spots{s};
        sname = 'Spinup_';
        lname = [loc '_'];
        %lname = ['phen_' loc '_'];
        SP = csvread([dpath sname lname 'Sml_p.csv'],0,20);
        SF = csvread([dpath sname lname 'Sml_f.csv'],0,20);
        SD = csvread([dpath sname lname 'Sml_d.csv'],0,20);
        MP = csvread([dpath sname lname 'Med_p.csv'],0,20);
        MF = csvread([dpath sname lname 'Med_f.csv'],0,20);
        MD = csvread([dpath sname lname 'Med_d.csv'],0,20);
        LP = csvread([dpath sname lname 'Lrg_p.csv'],0,20);
        LD = csvread([dpath sname lname 'Lrg_d.csv'],0,20);
        C = csvread([dpath sname lname 'Cobalt.csv']);
        z(:,1) = C(:,2);
        z(:,2) = C(:,3);
        z(:,3) = C(:,4);
        z=floor(z);
        
        % Plots over time
        x=1:length(SP);
        y=x/365;
        lstd=length(SP);
        id1 = 0:365:(lstd-1);
        id2 = 365:365:(lstd);
        ID  = [id1 id2];
        
        % Mean consumption level
        c=[SF(:,1) SP(:,1) SD(:,1) MF(:,1) MP(:,1) MD(:,1) LP(:,1) LD(:,1)];
        mclev(s,:) = nanmean(c);
        
        % Zoop overconsumption 
        Zcon(s,:) = nansum(z)/lstd;
           
    end
    %%
    save([dpath sname 'consump.mat'],'mclev','Zcon');
    csvwrite([dpath sname 'clevel.csv'],mclev);
    csvwrite([dpath sname 'Zconsump.csv'],Zcon);
    
end
