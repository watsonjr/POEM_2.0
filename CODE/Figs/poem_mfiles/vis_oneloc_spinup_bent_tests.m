% Visualize output of POEM
% Spinup at one location
% 100 years, but only last year saved

clear all
close all

datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

RE = {'1000','0500','0100','0050','0010','00050','00010','00005','00001'};
reff = [1.0,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001];
CarCap = {'025','050','075','100','125','150','175','200','225','250','275',...
    '300','325','350','375','400','425','450','475','500'};
car = [0.25:0.25:5.0];
benteff = {'05','10','15','20','25','30'};
beff = [0.05:0.05:0.3];
fcrit = 40;
nmort = 'M0';
kad = 100;
pref = 'D100';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),3);

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);

%%
for i=5:6;%1:length(benteff)
    BE = benteff{i};
    for C = 1:length(CarCap)
        CC = CarCap{C};
        for R = 1:length(RE)
            rfrac = RE{R};
            dp = ['NoDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_CC' CC '_RE' rfrac];
            dpath = [datap char(dp) '/'];
            fpath = [figp char(dp) '/'];
            cfile = char(dp);
            cfile2 = ['NoDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_CC' CC '_RE' rfrac...
                '_BentMortTests'];
            %%
            Psum = NaN*ones(3,length(spots));
            Fsum = NaN*ones(2,length(spots));
            Dsum = Psum;
            Bmean = NaN*ones(1,length(spots));
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
            Plev = Psum;
            Flev = Fsum;
            Dlev = Psum;
            all_mean=NaN*ones(3,4,length(spots));
            z = NaN*ones(length(spots),3);
            
            %%
            for s=1:length(spots)
                %%
                loc = spots{s};
                lname = [loc '_'];
                SP = csvread([dpath sname lname 'Sml_p.csv']);
                SF = csvread([dpath sname lname 'Sml_f.csv']);
                SD = csvread([dpath sname lname 'Sml_d.csv']);
                MP = csvread([dpath sname lname 'Med_p.csv']);
                MF = csvread([dpath sname lname 'Med_f.csv']);
                MD = csvread([dpath sname lname 'Med_d.csv']);
                LP = csvread([dpath sname lname 'Lrg_p.csv']);
                LD = csvread([dpath sname lname 'Lrg_d.csv']);
                CO = csvread([dpath sname lname 'Cobalt.csv']);
                
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
                B_mean =mean(CO(lyr,1));
                
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
                Bmean(:,s) = B_mean;
                
                all_mean(1:2,1,s) = F_mean;
                all_mean(:,2,s) = P_mean;
                all_mean(:,3,s) = D_mean;
                all_mean(1,4,s) = B_mean;
                
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
                
                %% Feeding level = con/cmax
                SP_lev=nanmean(SP(lyr,21));
                SF_lev=nanmean(SF(lyr,21));
                SD_lev=nanmean(SD(lyr,21));
                MP_lev=nanmean(MP(lyr,21));
                MF_lev=nanmean(MF(lyr,21));
                MD_lev=nanmean(MD(lyr,21));
                LP_lev=nanmean(LP(lyr,21));
                LD_lev=nanmean(LD(lyr,21));
                
                P_lev=[SP_lev;MP_lev;LP_lev];
                F_lev=[SF_lev;MF_lev];
                D_lev=[SD_lev;MD_lev;LD_lev];
                
                Plev(:,s) = P_lev;
                Flev(:,s) = F_lev;
                Dlev(:,s) = D_lev;
                
                %% Fraction zoop losses consumed
                z(s,1) = nanmean(CO(lyr,2));
                z(s,2) = nanmean(CO(lyr,3));
                z(s,3) = nanmean(CO(lyr,4));
                
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
            save([dpath sname 'lastyr_sum_mean_biom'],'Psum','Fsum',...
                'Dsum','Pmean','Fmean','Dmean','all_mean',...
                'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
                'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
                'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
                'Dtotcatch','Pgge','Fgge','Dgge','Plev','Flev','Dlev','Bmean');
            
        end
    end
end

%%
ndp = length(CarCap);
for i=5:6;%1:length(benteff)
    BE = benteff{i};
    for R = 1:length(RE)
        rfrac = RE{R};
        
        close all
        MDlev   = NaN*ones(length(CarCap),3);
        LDlev   = MDlev;
        BBCC    = MDlev;
        Bbio    = NaN*ones(length(CarCap),length(spots));
        MLD     = NaN*ones(length(CarCap),length(spots));
        fishsp  = NaN*ones(4,length(spots),length(CarCap));
        for C = 1:length(CarCap)
            CC = CarCap{C};
            cc = car(C);
            dp = ['NoDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_CC' CC '_RE' rfrac];
            dpath = [datap char(dp) '/'];
            fpath = [figp char(dp) '/'];
            cfile = char(dp);
            cfile2 = ['NoDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort'  nmort '_BE' BE '_RE' rfrac...
                '_BentCCtests'];
            
            load([dpath sname 'lastyr_sum_mean_biom']);
            
            %%
            MDlev(C,:)    = Dlev(2,[1,2,6]);
            LDlev(C,:)    = Dlev(3,[1,2,6]);
            BBCC(C,:)     = all_mean(1,4,[1,2,6]) ./ cc;
            Bbio(C,:)     = squeeze(all_mean(1,4,:));
            fishsp(:,:,C) = squeeze(nansum(all_mean));
            MLD(C,:)      = nansum(squeeze(all_mean(2:3,3,:)));
            
        end %CC
        
        %%
        for s=1:length(spots)
            loc = spots{s};
            lname = [loc '_'];
            
            %% B
            f1 = figure(1);
            subplot(4,3,s)
            plot(car,log10(Bbio(:,s)),'.k','MarkerSize',25); hold on;
            xlim([0 car(end)+0.25])
            ylim([-3 1])
            if (s==4)
                ylabel('log10 Mean Biom (g m^-^2) in final year')
            end
            set(gca,'XTick',0.5:0.5:car(end),'XTickLabel',[]);
            for t=0.5:0.5:car(end)
                text(t,-3.1,num2str(t),'Rotation',45,'HorizontalAlignment','right')
            end
            if (s==11)
                text(7,0,['BE=' num2str(beff(i))]);
                text(7,-1,['RE=' num2str(reff(R))]);
            end
            stamp(cfile2)
            title([loc ' B'])
            
            %% Sum mean biom over stages
            f2=figure(2);
            subplot(4,3,s)
            plot(car-0.05,log10(squeeze(fishsp(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
                'MarkerSize',15); hold on;
            plot(car,log10(squeeze(fishsp(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
                'MarkerSize',15); hold on;
            plot(car+0.05,log10(squeeze(fishsp(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
                'MarkerSize',15); hold on;
            xlim([0 car(end)+0.25])
            ylim([-2 1])
            if (s==4)
                ylabel('log10 Mean Biom (g m^-^2) in final year')
            end
            set(gca,'XTick',0.5:0.5:car(end),'XTickLabel',[]);
            for t=0.5:0.5:car(end)
                text(t,-2.1,num2str(t),'Rotation',45,'HorizontalAlignment','right')
            end
            if (s==11)
                text(7,0,['BE=' num2str(beff(i))]);
                text(7,-1,['RE=' num2str(reff(R))]);
            end
            stamp(cfile2)
            title([loc ' All stages'])
            
        end %spots
        print(f1,'-dpng',[figp sname cfile2 '_tot_mean_bent_biomass_all_locs.png'])
        print(f2,'-dpng',[figp sname cfile2 '_tot_mean_biomass_type_all_locs.png'])
        
        %% Feeding level
        %MD
        f5=figure(5);
        bar(MDlev)
        ylim([0 1])
        xlim([0 ndp+1])
        set(gca,'XTick',2:2:ndp,'XTickLabel',car(2:2:end))
        legend('GB','EBS','NS')
        legend('location','northwest')
        ylabel('MD C / Cmax')
        xlabel('CC')
        title(['MD feeding level with BE=' num2str(beff(i)) ' and RE=' num2str(reff(R))])
        colormap(cmap3)
        print('-dpng',[figp sname cfile2 '_MD_flev.png'])
        
        %LD
        f6=figure(6);
        bar(LDlev)
        ylim([0 1])
        xlim([0 ndp+1])
        set(gca,'XTick',2:2:ndp,'XTickLabel',car(2:2:end))
        legend('GB','EBS','NS')
        legend('location','northwest')
        ylabel('LD C / Cmax')
        xlabel('CC')
        title(['LD feeding level with BE=' num2str(beff(i)) ' and RE=' num2str(reff(R))])
        colormap(cmap3)
        print('-dpng',[figp sname cfile2 '_LD_flev.png'])
        
        %% B/B_CC
        
        f7=figure(7);
        bar(BBCC)
        xlim([0 ndp+1])
        set(gca,'XTick',2:2:ndp,'XTickLabel',car(2:2:end))
        legend('GB','EBS','NS')
        legend('location','northeast')
        ylabel('B / CC')
        xlabel('CC')
        title(['B/CC with BE=' num2str(beff(i)) ' and RE=' num2str(reff(R))])
        colormap(cmap3)
        print('-dpng',[figp sname cfile2 '_BBCC.png'])
        
        %% Save values for all locs and all CC for that RE and BE combo
        save([datap 'Bent_CC_tests/' cfile2 '.mat'],'MDlev','LDlev','BBCC',...
            'Bbio','fishsp','MLD')
        
    end %RE
end %BE
