% Visualize output of POEM
% Fishing spinup at one location
% 50 years

clear all
close all

datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

fsim = {'0','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
frate = {'000','01','02','03','04','05','06','07','08','09','10'};
% sims = {'1','0.5','0.1','0.05','0.01','5e-3','1e-3','5e-4','1e-4','5e-5','1e-5'};
% rfrac = {'1000','0500','0100','0050','0010','00050','00010','00005',...
%     '00001','000005','000001'};
% sims = {'4e-3','3e-3','2e-3'};
% rfrac = {'00040','00030','00020'};
% sims = {'1.75e-3','1.5e-3','1.25e-3'};
% rfrac = {'00018','00015','00012'};
sims = {'4e-3','3e-3','2e-3','1.75e-3','1.5e-3','1.25e-3'};
rfrac = {'00040','00030','00020','00018','00015','00012'};
fcrit = 40;
nmort = 'M2';
kad = 100;
pref = 'D050';

sname = 'Spinup_';
sname2 = '';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

load('cmap_ppt_angles.mat')

%%
for r = 1:length(rfrac)
    close all
    RE = rfrac{r};
    for i=1:length(frate)
        F = frate{i};
        dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort'  nmort '_BE05_RE' num2str(RE) '_MF_fish' F];
%         dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
%             '_MZ01_nmort'  nmort '_BE05_RE' num2str(RE) '_BAassim_MF_fish' F];
        dpath = [datap char(dp) '/'];
        fpath = [figp char(dp) '/'];
        cfile = char(dp);
        cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
            '_' pref '_nmort'  nmort '_BE05_RE' num2str(RE) '_MF_fishing_catch'];
%         cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
%             '_MZ01_nmort'  nmort '_BE05_RE' num2str(RE) '_BAassim_MF_fishing_catch'];
        
        %%
        Psum = NaN*ones(3,length(spots));
        Fsum = NaN*ones(2,length(spots));
        Dsum = Psum;
        Pmean = Psum;
        Fmean = Fsum;
        Dmean = Psum;
        Pmgr = Psum;
        Fmgr = Fsum;
        Dmgr = Psum;
        Pgge = Psum;
        Fgge = Fsum;
        Dgge = Psum;
        Pprod = Psum;
        Fprod = Fsum;
        Dprod = Psum;
        Pcon = Psum;
        Fcon = Fsum;
        Dcon = Psum;
        Prep = Fsum;
        Frep = Fsum;
        Drep = Fsum;
        Pmet = Psum;
        Fmet = Fsum;
        Dmet = Psum;
        Ppred = Psum;
        Fpred = Fsum;
        Dpred = Psum;
        Pnat = Psum;
        Fnat = Fsum;
        Dnat = Psum;
        Pfish = Psum;
        Ffish = Fsum;
        Dfish = Psum;
        Ptotcatch = Psum;
        Ftotcatch = Fsum;
        Dtotcatch = Psum;
        all_mean=NaN*ones(3,3,length(spots));
        z = NaN*ones(length(spots),3);
        
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
            %lyr=t((end-365+1):end);
            lyr=1:365;
            
            %% Final mean biomass in each size
            
            SP_sum=sum(SP(lyr,1));
            SF_sum=sum(SF(lyr,1));
            SD_sum=sum(SD(lyr,1));
            MP_sum=sum(MP(lyr,1));
            MF_sum=sum(MF(lyr,1));
            MD_sum=sum(MD(lyr,1));
            LP_sum=sum(LP(lyr,1));
            LD_sum=sum(LD(lyr,1));
            
            SP_mean=mean(SP(lyr,1));
            SF_mean=mean(SF(lyr,1));
            SD_mean=mean(SD(lyr,1));
            MP_mean=mean(MP(lyr,1));
            MF_mean=mean(MF(lyr,1));
            MD_mean=mean(MD(lyr,1));
            LP_mean=mean(LP(lyr,1));
            LD_mean=mean(LD(lyr,1));
            
            P_sum=[SP_sum;MP_sum;LP_sum];
            F_sum=[SF_sum;MF_sum];
            D_sum=[SD_sum;MD_sum;LD_sum];
            P_mean=[SP_mean;MP_mean;LP_mean];
            F_mean=[SF_mean;MF_mean];
            D_mean=[SD_mean;MD_mean;LD_mean];
            
            Psum(:,s) = P_sum;
            Fsum(:,s) = F_sum;
            Dsum(:,s) = D_sum;
            Pmean(:,s) = P_mean;
            Fmean(:,s) = F_mean;
            Dmean(:,s) = D_mean;
            
            
            all_mean(1:2,1,s) = F_mean;
            all_mean(:,2,s) = P_mean;
            all_mean(:,3,s) = D_mean;
            
            %% Growth rate (nu - energy for biomass production)
            SP_mgr=nanmean(SP(lyr,15));
            SF_mgr=nanmean(SF(lyr,15));
            SD_mgr=nanmean(SD(lyr,15));
            MP_mgr=nanmean(MP(lyr,15));
            MF_mgr=nanmean(MF(lyr,15));
            MD_mgr=nanmean(MD(lyr,15));
            LP_mgr=nanmean(LP(lyr,15));
            LD_mgr=nanmean(LD(lyr,15));
            
            P_mgr=[SP_mgr;MP_mgr;LP_mgr];
            F_mgr=[SF_mgr;MF_mgr];
            D_mgr=[SD_mgr;MD_mgr;LD_mgr];
            
            Pmgr(:,s) = P_mgr;
            Fmgr(:,s) = F_mgr;
            Dmgr(:,s) = D_mgr;
            
            %% Consump per biomass (I)
            SP_con=nanmean(SP(lyr,14));
            SF_con=nanmean(SF(lyr,14));
            SD_con=nanmean(SD(lyr,14));
            MP_con=nanmean(MP(lyr,14));
            MF_con=nanmean(MF(lyr,14));
            MD_con=nanmean(MD(lyr,14));
            LP_con=nanmean(LP(lyr,14));
            LD_con=nanmean(LD(lyr,14));
            
            P_con=[SP_con;MP_con;LP_con];
            F_con=[SF_con;MF_con];
            D_con=[SD_con;MD_con;LD_con];
            
            Pcon(:,s) = P_con;
            Fcon(:,s) = F_con;
            Dcon(:,s) = D_con;
            
            %% Fraction zoop losses consumed
            z(s,1) = nanmean(C(lyr,2));
            z(s,2) = nanmean(C(lyr,3));
            z(s,3) = nanmean(C(lyr,4));
            
            %% Size spectrum (sum stages)
            spec = nansum(all_mean(:,:,s),2);
            
            %% Production (= nu * biom)
            SP_prod=nanmean(SP(lyr,24));
            SF_prod=nanmean(SF(lyr,24));
            SD_prod=nanmean(SD(lyr,24));
            MP_prod=nanmean(MP(lyr,24));
            MF_prod=nanmean(MF(lyr,24));
            MD_prod=nanmean(MD(lyr,24));
            LP_prod=nanmean(LP(lyr,24));
            LD_prod=nanmean(LD(lyr,24));
            
            P_prod=[SP_prod;MP_prod;LP_prod];
            F_prod=[SF_prod;MF_prod];
            D_prod=[SD_prod;MD_prod;LD_prod];
            
            Pprod(:,s) = P_prod;
            Fprod(:,s) = F_prod;
            Dprod(:,s) = D_prod;
            
            %% Reproduction
            F_rep(1,1)=nanmean(MF(lyr,18));
            D_rep(1,1)=nanmean(LD(lyr,18));
            P_rep(1,1)=nanmean(LP(lyr,18));
            F_rep(2,1)=nanmean(MF(lyr,1).*MF(lyr,18));
            D_rep(2,1)=nanmean(LD(lyr,1).*LD(lyr,18));
            P_rep(2,1)=nanmean(LP(lyr,1).*LP(lyr,18));
            
            Prep(:,s) = P_rep;
            Frep(:,s) = F_rep;
            Drep(:,s) = D_rep;
            
            %% Metabolism
            SP_met=nanmean(SP(lyr,27));
            SF_met=nanmean(SF(lyr,27));
            SD_met=nanmean(SD(lyr,27));
            MP_met=nanmean(MP(lyr,27));
            MF_met=nanmean(MF(lyr,27));
            MD_met=nanmean(MD(lyr,27));
            LP_met=nanmean(LP(lyr,27));
            LD_met=nanmean(LD(lyr,27));
            
            P_met=[SP_met;MP_met;LP_met];
            F_met=[SF_met;MF_met];
            D_met=[SD_met;MD_met;LD_met];
            
            Pmet(:,s) = P_met;
            Fmet(:,s) = F_met;
            Dmet(:,s) = D_met;
            
            %% Predation
            SP_pred=nanmean(SP(lyr,25));
            SF_pred=nanmean(SF(lyr,25));
            SD_pred=nanmean(SD(lyr,25));
            MP_pred=nanmean(MP(lyr,25));
            MF_pred=nanmean(MF(lyr,25));
            MD_pred=nanmean(MD(lyr,25));
            LP_pred=nanmean(LP(lyr,25));
            LD_pred=nanmean(LD(lyr,25));
            
            P_pred=[SP_pred;MP_pred;LP_pred];
            F_pred=[SF_pred;MF_pred];
            D_pred=[SD_pred;MD_pred;LD_pred];
            
            Ppred(:,s) = P_pred;
            Fpred(:,s) = F_pred;
            Dpred(:,s) = D_pred;
            
            %% Natural mortality
            Pnat(1,s)=nanmean(SP(lyr,26));
            Fnat(1,s)=nanmean(SF(lyr,26));
            Dnat(1,s)=nanmean(SD(lyr,26));
            Pnat(2,s)=nanmean(MP(lyr,26));
            Fnat(2,s)=nanmean(MF(lyr,26));
            Dnat(2,s)=nanmean(MD(lyr,26));
            Pnat(3,s)=nanmean(LP(lyr,26));
            Dnat(3,s)=nanmean(LD(lyr,26));
            
            %% Fishing
            MP_fish=nanmean(MP(lyr,28));
            MF_fish=nanmean(MF(lyr,28));
            MD_fish=nanmean(MD(lyr,28));
            LP_fish=nanmean(LP(lyr,28));
            LD_fish=nanmean(LD(lyr,28));
            
            P_fish=[0;MP_fish;LP_fish];
            F_fish=[0;MF_fish];
            D_fish=[0;MD_fish;LD_fish];
            
            Pfish(:,s) = P_fish;
            Ffish(:,s) = F_fish;
            Dfish(:,s) = D_fish;
            
            MP_totcatch=nansum(MP(lyr,28));
            MF_totcatch=nansum(MF(lyr,28));
            MD_totcatch=nansum(MD(lyr,28));
            LP_totcatch=nansum(LP(lyr,28));
            LD_totcatch=nansum(LD(lyr,28));
            
            P_totcatch=[0;MP_totcatch;LP_totcatch];
            F_totcatch=[0;MF_totcatch];
            D_totcatch=[0;MD_totcatch;LD_totcatch];
            
            Ptotcatch(:,s) = P_totcatch;
            Ftotcatch(:,s) = F_totcatch;
            Dtotcatch(:,s) = D_totcatch;
            
            %% Total mortality w/o fishing
            %model lengths
            L(1) = 10^((log10(2)+log10(20))/2);
            L(2) = 10^((log10(20)+log10(200))/2);
            L(3) = 10^((log10(200)+log10(2000))/2);
            %model mass in grams
            M = 0.01 .* (0.1.*L).^3;
            %Andersen & Beyer mortality rate per year (natural + predation)
            %physiol mort * growth constant * M^-0.25
            AB = (0.35 .* 4.5 .* M.^(-0.25)) ./365;
            
            Fmort = Fpred + Fnat;
            Pmort = Ppred + Pnat;
            Dmort = Dpred + Dnat;
            
            %% Total mortality w/ fishing
            Fmortf = Fpred + Fnat + Ffish;
            Pmortf = Ppred + Pnat + Pfish;
            Dmortf = Dpred + Dnat + Dfish;
            
            %% Gross growth efficiency (= nu/consump)
            SP_gge=nanmean(SP(lyr,15)./SP(lyr,14));
            SF_gge=nanmean(SF(lyr,15)./SF(lyr,14));
            SD_gge=nanmean(SD(lyr,15)./SD(lyr,14));
            MP_gge=nanmean(MP(lyr,15)./MP(lyr,14));
            MF_gge=nanmean(MF(lyr,15)./MF(lyr,14));
            MD_gge=nanmean(MD(lyr,15)./MD(lyr,14));
            LP_gge=nanmean(LP(lyr,15)./LP(lyr,14));
            LD_gge=nanmean(LD(lyr,15)./LD(lyr,14));
            
            P_gge=[SP_gge;MP_gge;LP_gge];
            F_gge=[SF_gge;MF_gge];
            D_gge=[SD_gge;MD_gge;LD_gge];
            
            Pgge(:,s) = P_gge;
            Fgge(:,s) = F_gge;
            Dgge(:,s) = D_gge;
            
        end
        
        save([dpath sname sname2 'lastyr_sum_mean_biom'],'Psum','Fsum',...
            'Dsum','Pmean','Fmean','Dmean','all_mean',...
            'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
            'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
            'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
            'Dtotcatch','Pgge','Fgge','Dgge');
        
    end
    
    %%
    for s=1:length(spots)
        
        loc = spots{s};
        lname = [sname2 loc '_'];
        close all
        %%
        for i=1:length(frate)
            F = frate{i};
            dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE05_RE' num2str(RE) '_MF_fish' F];
%             dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
%                 '_MZ01_nmort'  nmort '_BE05_RE' num2str(RE) '_BAassim_MF_fish' F];
            dpath = [char(dp) '/'];
            load([datap dpath sname sname2 'lastyr_sum_mean_biom']);
            
            
            %% Fishing
            Ftot = sum(Ftotcatch(:,s));
            Ptot = sum(Ptotcatch(:,s));
            Dtot = sum(Dtotcatch(:,s));
            Tot = Ftot+Ptot+Dtot;
            
            f7 = figure(7);
            plot(i,Ftot,'.k','MarkerSize',20); hold on;
            xlim([0 length(frate)+1])
            if (i==length(frate))
                set(gca,'XTick',1:length(frate),'XTickLabel',fsim);
                ylabel('All F')
                title(['RE=' sims{r} ' ' loc])
                stamp(cfile2)
            end
        end
        print(f7,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'F.png'])
        
    end
    
    %%
    for s=1:length(spots)
        
        loc = spots{s};
        lname = [sname2 loc '_'];
        
        %%
        for i=1:length(frate)
            F = frate{i};
            dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE05_RE' num2str(RE) '_MF_fish' F];
%             dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
%                 '_MZ01_nmort'  nmort '_BE05_RE' num2str(RE) '_BAassim_MF_fish' F];
            dpath = [char(dp) '/'];
            load([datap dpath sname sname2 'lastyr_sum_mean_biom']);
            
            
            %% Fishing
            Ftot = sum(Ftotcatch(:,s));
            Ptot = sum(Ptotcatch(:,s));
            Dtot = sum(Dtotcatch(:,s));
            Tot = Ftot+Ptot+Dtot;
            
            f12 = figure(12);
            subplot(4,3,s)
            plot(i,Ftot,'.k','MarkerSize',25); hold on;
            xlim([0 length(frate)+1])
            if (i==length(frate))
                set(gca,'XTick',1:length(frate),'XTickLabel',fsim);
                if (s==2)
                    str = {['RE=' sims{r}], loc};
                    title(str)
                else
                    title(loc)
                end
                if (s==4)
                    ylabel('Total F catch (g) in final year')
                end
                stamp(cfile2)
            end
            
        end
        
    end
    print(f12,'-dpng',[figp sname sname2 cfile2 '_allF_locs.png'])
end
