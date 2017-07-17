% Visualize output of POEM fishing tests
% Spinup at all single locations
% 150 years, monthly means saved
% All fish harvested

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/Fishing/';
%figp = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/Fishing/';

fnum = [0:0.1:1];%,1.2:0.2:2];
fsim = {'0','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};%,'1.2','1.4','1.6','1.8','2'};
frate = {'0','01','02','03','04','05','06','07','08','09','10'};%,'12','14','16','18','20'};
efn = 70;
tefn = num2str(efn);
cfn = 20;
tcfn = num2str(cfn);
fcrit = 20;
nmort = '1';
kad = 50;
Dprefs = 0.1:0.1:1;
Jprefs = 0.5:0.1:1;
Aprefs = 0.5:0.1:1;
Sprefs = 0.05:0.05:0.5;
D = 0.75;
J = 1;
Ad = 0.5;
Sm = 0.25;
td = num2str(1000+int64(100*D));
tj = num2str(1000+int64(100*J));
ta = num2str(1000+int64(100*Ad));
tsm = num2str(1000+int64(100*Sm));
BE = '05';
CC = '050';
rfrac = '00100';
rfrac2 = '00400';
sel = 'All';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),3);

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
%load('/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);

%%
for i=1:length(frate)
    F = frate{i};
    dp = ['Dc_enc',tefn,'_cmax-metab',tcfn,'_fcrit',num2str(fcrit),'_D',td(2:end),...
        '_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',nmort,...
        '_BE',BE,'_CC',CC,'_lgRE',rfrac,'_mdRE',rfrac2];
    dpath = [datap char(dp) '/'];
    %fpath = [figp char(dp) '/'];
    fpath = figp;
    cfile = char(dp);
    
    %%
    all_mean=NaN*ones(3,4,length(spots));
    z = NaN*ones(length(spots),3);
    
    if F=='0'
        load([dpath sname 'locs.mat'])
    else
        load([dpath sname 'locs','_',sel,'_fish',F,'.mat'])
    end
    
    SP = S_Sml_p;
    SF = S_Sml_f;
    SD = S_Sml_d;
    MP = S_Med_p;
    MF = S_Med_f;
    MD = S_Med_d;
    LP = S_Lrg_p;
    LD = S_Lrg_d;
    CO = S_Cobalt;
    
    t=1:size(SP,1);
    lyr=t((end-12+1):end);
    
    %% Final mean biomass in each size
    
    SP_mean=squeeze(mean(SP(lyr,1,:)))';
    SF_mean=squeeze(mean(SF(lyr,1,:)))';
    SD_mean=squeeze(mean(SD(lyr,1,:)))';
    MP_mean=squeeze(mean(MP(lyr,1,:)))';
    MF_mean=squeeze(mean(MF(lyr,1,:)))';
    MD_mean=squeeze(mean(MD(lyr,1,:)))';
    LP_mean=squeeze(mean(LP(lyr,1,:)))';
    LD_mean=squeeze(mean(LD(lyr,1,:)))';
    B_mean =squeeze(mean(CO(lyr,1,:)))';
    
    Pmean=[SP_mean;MP_mean;LP_mean];
    Fmean=[SF_mean;MF_mean];
    Dmean=[SD_mean;MD_mean;LD_mean];
    Bmean = B_mean;
    
    all_mean(1:2,1,:) = Fmean;
    all_mean(:,2,:) = Pmean;
    all_mean(:,3,:) = Dmean;
    all_mean(1,4,:) = Bmean;
    
    % Size spectrum (sum stages)
    spec = nansum(all_mean,2);
    
    %% Growth rate (nu - energy for biomass production)
    SP_mgr=squeeze(nanmean(SP(lyr,15,:)))';
    SF_mgr=squeeze(nanmean(SF(lyr,15,:)))';
    SD_mgr=squeeze(nanmean(SD(lyr,15,:)))';
    MP_mgr=squeeze(nanmean(MP(lyr,15,:)))';
    MF_mgr=squeeze(nanmean(MF(lyr,15,:)))';
    MD_mgr=squeeze(nanmean(MD(lyr,15,:)))';
    LP_mgr=squeeze(nanmean(LP(lyr,15,:)))';
    LD_mgr=squeeze(nanmean(LD(lyr,15,:)))';
    
    Pmgr=[SP_mgr;MP_mgr;LP_mgr];
    Fmgr=[SF_mgr;MF_mgr];
    Dmgr=[SD_mgr;MD_mgr;LD_mgr];
    
    %% Consump per biomass (I)
    SP_con=squeeze(nanmean(SP(lyr,14,:)))';
    SF_con=squeeze(nanmean(SF(lyr,14,:)))';
    SD_con=squeeze(nanmean(SD(lyr,14,:)))';
    MP_con=squeeze(nanmean(MP(lyr,14,:)))';
    MF_con=squeeze(nanmean(MF(lyr,14,:)))';
    MD_con=squeeze(nanmean(MD(lyr,14,:)))';
    LP_con=squeeze(nanmean(LP(lyr,14,:)))';
    LD_con=squeeze(nanmean(LD(lyr,14,:)))';
    
    Pcon=[SP_con;MP_con;LP_con];
    Fcon=[SF_con;MF_con];
    Dcon=[SD_con;MD_con;LD_con];
    
    %% Feeding level = con/cmax
    SP_lev=squeeze(nanmean(SP(lyr,20,:)))';
    SF_lev=squeeze(nanmean(SF(lyr,20,:)))';
    SD_lev=squeeze(nanmean(SD(lyr,20,:)))';
    MP_lev=squeeze(nanmean(MP(lyr,20,:)))';
    MF_lev=squeeze(nanmean(MF(lyr,20,:)))';
    MD_lev=squeeze(nanmean(MD(lyr,20,:)))';
    LP_lev=squeeze(nanmean(LP(lyr,20,:)))';
    LD_lev=squeeze(nanmean(LD(lyr,20,:)))';
    
    Plev=[SP_lev;MP_lev;LP_lev];
    Flev=[SF_lev;MF_lev];
    Dlev=[SD_lev;MD_lev;LD_lev];
    
    %% Fraction zoop losses consumed
    z(:,1) = squeeze(nanmean(CO(lyr,3,:)));
    z(:,2) = squeeze(nanmean(CO(lyr,4,:)));
    z(:,3) = squeeze(nanmean(CO(lyr,5,:)));
    
    %% Production (= nu * biom)
    SP_prod=squeeze(nanmean(SP(lyr,21,:)))';
    SF_prod=squeeze(nanmean(SF(lyr,21,:)))';
    SD_prod=squeeze(nanmean(SD(lyr,21,:)))';
    MP_prod=squeeze(nanmean(MP(lyr,21,:)))';
    MF_prod=squeeze(nanmean(MF(lyr,21,:)))';
    MD_prod=squeeze(nanmean(MD(lyr,21,:)))';
    LP_prod=squeeze(nanmean(LP(lyr,21,:)))';
    LD_prod=squeeze(nanmean(LD(lyr,21,:)))';
    
    Pprod=[SP_prod;MP_prod;LP_prod];
    Fprod=[SF_prod;MF_prod];
    Dprod=[SD_prod;MD_prod;LD_prod];
    
    %% Reproduction
    Frep(1,:)=squeeze(nanmean(MF(lyr,18,:)))';
    Drep(1,:)=squeeze(nanmean(LD(lyr,18,:)))';
    Prep(1,:)=squeeze(nanmean(LP(lyr,18,:)))';
    Frep(2,:)=squeeze(nanmean(MF(lyr,1,:).*MF(lyr,18,:)))';
    Drep(2,:)=squeeze(nanmean(LD(lyr,1,:).*LD(lyr,18,:)))';
    Prep(2,:)=squeeze(nanmean(LP(lyr,1,:).*LP(lyr,18,:)))';
    
    %% Metabolism
    SP_met=squeeze(nanmean(SP(lyr,24,:)))';
    SF_met=squeeze(nanmean(SF(lyr,24,:)))';
    SD_met=squeeze(nanmean(SD(lyr,24,:)))';
    MP_met=squeeze(nanmean(MP(lyr,24,:)))';
    MF_met=squeeze(nanmean(MF(lyr,24,:)))';
    MD_met=squeeze(nanmean(MD(lyr,24,:)))';
    LP_met=squeeze(nanmean(LP(lyr,24,:)))';
    LD_met=squeeze(nanmean(LD(lyr,24,:)))';
    
    Pmet=[SP_met;MP_met;LP_met];
    Fmet=[SF_met;MF_met];
    Dmet=[SD_met;MD_met;LD_met];
    
    %% Predation
    SP_pred=squeeze(nanmean(SP(lyr,22,:)))';
    SF_pred=squeeze(nanmean(SF(lyr,22,:)))';
    SD_pred=squeeze(nanmean(SD(lyr,22,:)))';
    MP_pred=squeeze(nanmean(MP(lyr,22,:)))';
    MF_pred=squeeze(nanmean(MF(lyr,22,:)))';
    MD_pred=squeeze(nanmean(MD(lyr,22,:)))';
    LP_pred=squeeze(nanmean(LP(lyr,22,:)))';
    LD_pred=squeeze(nanmean(LD(lyr,22,:)))';
    
    Ppred=[SP_pred;MP_pred;LP_pred];
    Fpred=[SF_pred;MF_pred];
    Dpred=[SD_pred;MD_pred;LD_pred];
    
    %% Natural mortality
    Pnat(1,:)=squeeze(nanmean(SP(lyr,23,:)))';
    Fnat(1,:)=squeeze(nanmean(SF(lyr,23,:)))';
    Dnat(1,:)=squeeze(nanmean(SD(lyr,23,:)))';
    Pnat(2,:)=squeeze(nanmean(MP(lyr,23,:)))';
    Fnat(2,:)=squeeze(nanmean(MF(lyr,23,:)))';
    Dnat(2,:)=squeeze(nanmean(MD(lyr,23,:)))';
    Pnat(3,:)=squeeze(nanmean(LP(lyr,23,:)))';
    Dnat(3,:)=squeeze(nanmean(LD(lyr,23,:)))';
    
    %% Fishing
    MP_fish=squeeze(nanmean(MP(lyr,25,:)))';
    MF_fish=squeeze(nanmean(MF(lyr,25,:)))';
    MD_fish=squeeze(nanmean(MD(lyr,25,:)))';
    LP_fish=squeeze(nanmean(LP(lyr,25,:)))';
    LD_fish=squeeze(nanmean(LD(lyr,25,:)))';
    
    Pfish=[zeros(size(MP_fish));MP_fish;LP_fish];
    Ffish=[zeros(size(MF_fish));MF_fish];
    Dfish=[zeros(size(MD_fish));MD_fish;LD_fish];
    
    MP_totcatch=squeeze(nansum(MP(lyr,25,:)))';
    MF_totcatch=squeeze(nansum(MF(lyr,25,:)))';
    MD_totcatch=squeeze(nansum(MD(lyr,25,:)))';
    LP_totcatch=squeeze(nansum(LP(lyr,25,:)))';
    LD_totcatch=squeeze(nansum(LD(lyr,25,:)))';
    
    Ptotcatch=[zeros(size(MP_totcatch));MP_totcatch;LP_totcatch];
    Ftotcatch=[zeros(size(MF_totcatch));MF_totcatch];
    Dtotcatch=[zeros(size(MD_totcatch));MD_totcatch;LD_totcatch];
    
    %% Total mortality w/o fishing
    Fmort = Fpred + Fnat;
    Pmort = Ppred + Pnat;
    Dmort = Dpred + Dnat;
    
    %% Total mortality w/ fishing
    Fmortf = Fpred + Fnat + Ffish;
    Pmortf = Ppred + Pnat + Pfish;
    Dmortf = Dpred + Dnat + Dfish;
    
    %% Gross growth efficiency (= nu/consump)
    SP_gge=squeeze(nanmean(SP(lyr,15,:)./SP(lyr,14,:)))';
    SF_gge=squeeze(nanmean(SF(lyr,15,:)./SF(lyr,14,:)))';
    SD_gge=squeeze(nanmean(SD(lyr,15,:)./SD(lyr,14,:)))';
    MP_gge=squeeze(nanmean(MP(lyr,15,:)./MP(lyr,14,:)))';
    MF_gge=squeeze(nanmean(MF(lyr,15,:)./MF(lyr,14,:)))';
    MD_gge=squeeze(nanmean(MD(lyr,15,:)./MD(lyr,14,:)))';
    LP_gge=squeeze(nanmean(LP(lyr,15,:)./LP(lyr,14,:)))';
    LD_gge=squeeze(nanmean(LD(lyr,15,:)./LD(lyr,14,:)))';
    
    Pgge=[SP_gge;MP_gge;LP_gge];
    Fgge=[SF_gge;MF_gge];
    Dgge=[SD_gge;MD_gge;LD_gge];
    
    if (F=='0')
        save([dpath sname 'lastyr_sum_mean_biom.mat'],'Pmean','Fmean','Dmean','all_mean',...
            'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
            'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
            'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
            'Dtotcatch','Pgge','Fgge','Dgge','Plev','Flev','Dlev','Bmean');
    else
        save([dpath sname 'lastyr_sum_mean_biom','_',sel,'_fish',F,'.mat'],...
            'Pmean','Fmean','Dmean','all_mean',...
            'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
            'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
            'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
            'Dtotcatch','Pgge','Fgge','Dgge','Plev','Flev','Dlev','Bmean');
    end
    
end

%%
cfile2 = [dp,'_',sel,'_fishing_catch'];
Fcaught  = NaN*ones(2,length(spots),length(frate));
Pcaught  = NaN*ones(2,length(spots),length(frate));
Dcaught  = NaN*ones(2,length(spots),length(frate));
fishsp  = NaN*ones(4,length(spots),length(frate));
for i=1:length(frate)
    F = frate{i};
    if F=='0'
        load([dpath sname 'lastyr_sum_mean_biom.mat'])
    else
        load([dpath sname 'lastyr_sum_mean_biom','_',sel,'_fish',F,'.mat'])
    end
    % Biomass of each type
    fishsp(:,:,i) = squeeze(nansum(all_mean));
    % Catch
    Fcaught(1,:,i) = Ftotcatch(2,:);
        Fcaught(2,:,i) = zeros(1,length(spots));
    Pcaught(:,:,i) = Ptotcatch(2:3,:);
    Dcaught(:,:,i) = Dtotcatch(2:3,:);
end
caught = Fcaught + Pcaught + Dcaught;
save([dpath cfile2 '.mat'],'caught','Fcaught','Pcaught','Dcaught',...
    'fnum','fishsp');

%%
np = length(frate);
for s=1:length(spots)
    
    loc = spots{s};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f1=figure(1);
    subplot(4,3,s)
    plot(1-0.25:np,log10(squeeze(fishsp(1,s,:))),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:np,log10(squeeze(fishsp(2,s,:))),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1+0.25:np+1,log10(squeeze(fishsp(3,s,:))),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 np+1])
    ylim([-2 2])
    set(gca,'XTick',1:np,'XTickLabel',[]);
    for t=1:np
        text(t,-2.1,fsim(t),'Rotation',45,'HorizontalAlignment','right')
    end
    if (s==4)
        ylabel('log10 Mean Biom (g m^-^2) in final year')
    end
    if (s==11)
        text(13,1.5,['D=' num2str(D)]);
        text(13,1,['A=' num2str(Ad)]);
        text(13,0.5,['J=' num2str(J)]);
        text(13,0,['Sm=' num2str(Sm)]);
        text(13,-0.5,['RE=' rfrac]);
        text(13,-1,['enc=' tefn]);
        text(13,-1.5,['cmax-metab=' tcfn]);
    end
    if (s==11)
        xlabel('Annual fishing rate')
    end
    title([loc ' All stages'])
    stamp(cfile2)
    
    
    %% Fishing
    %F
    f2 = figure(2);
    subplot(4,3,s)
    plot(fnum,squeeze(Fcaught(1,s,:)),'.k','MarkerSize',25); hold on;
    xlim([0 fnum(end)+0.1])
    title(loc)
    if (s==4)
        ylabel('MF catch (g) in final year')
    end
    if (s==11)
        xlabel('Annual fishing rate')
    end
    stamp(cfile2)
    
    %P
    f3 = figure(3);
    subplot(4,3,s)
    plot(fnum,squeeze(Pcaught(1,s,:)),'.k','MarkerSize',25); hold on;
    xlim([0 fnum(end)+0.1])
    title(loc)
    if (s==4)
        ylabel('MP catch (g) in final year')
    end
    if (s==11)
        xlabel('Annual fishing rate')
    end
    stamp(cfile2)
    
    f4 = figure(4);
    subplot(4,3,s)
    plot(fnum,squeeze(Pcaught(2,s,:)),'.k','MarkerSize',25); hold on;
    xlim([0 fnum(end)+0.1])
    title(loc)
    if (s==4)
        ylabel('LP catch (g) in final year')
    end
    if (s==11)
        xlabel('Annual fishing rate')
    end
    stamp(cfile2)
    
    f5 = figure(5);
    subplot(4,3,s)
    plot(fnum,squeeze(sum(Pcaught(:,s,:))),'.k','MarkerSize',25); hold on;
    xlim([0 fnum(end)+0.1])
    title(loc)
    if (s==4)
        ylabel('Total P catch (g) in final year')
    end
    if (s==11)
        xlabel('Annual fishing rate')
    end
    stamp(cfile2)
    
    %D
    f6 = figure(6);
    subplot(4,3,s)
    plot(fnum,squeeze(Dcaught(1,s,:)),'.k','MarkerSize',25); hold on;
    xlim([0 fnum(end)+0.1])
    title(loc)
    if (s==4)
        ylabel('MD catch (g) in final year')
    end
    if (s==11)
        xlabel('Annual fishing rate')
    end
    stamp(cfile2)
    
    f7 = figure(7);
    subplot(4,3,s)
    plot(fnum,squeeze(Dcaught(2,s,:)),'.k','MarkerSize',25); hold on;
    xlim([0 fnum(end)+0.1])
    title(loc)
    if (s==4)
        ylabel('LD catch (g) in final year')
    end
    if (s==11)
        xlabel('Annual fishing rate')
    end
    stamp(cfile2)
    
    f8 = figure(8);
    subplot(4,3,s)
    plot(fnum,squeeze(sum(Dcaught(:,s,:))),'.k','MarkerSize',25); hold on;
    xlim([0 fnum(end)+0.1])
    title(loc)
    if (s==4)
        ylabel('Total D catch (g) in final year')
    end
    if (s==11)
        xlabel('Annual fishing rate')
    end
    stamp(cfile2)
    
    %All
    f9 = figure(9);
    subplot(4,3,s)
    plot(fnum,squeeze(sum(caught(:,s,:))),'.k','MarkerSize',25); hold on;
    xlim([0 fnum(end)+0.1])
    title(loc)
    if (s==4)
        ylabel('Total catch (g) in final year')
    end
    if (s==11)
        xlabel('Annual fishing rate')
    end
    stamp(cfile2)
    
end %spots
print(f1,'-dpng',[figp sname cfile2 '_tot_mean_biomass_type_locs.png'])
print(f2,'-dpng',[figp sname cfile2 '_MF_locs.png'])
print(f3,'-dpng',[figp sname cfile2 '_MP_locs.png'])
print(f6,'-dpng',[figp sname cfile2 '_MD_locs.png'])
print(f4,'-dpng',[figp sname cfile2 '_LP_locs.png'])
print(f7,'-dpng',[figp sname cfile2 '_LD_locs.png'])
print(f5,'-dpng',[figp sname cfile2 '_allP_locs.png'])
print(f8,'-dpng',[figp sname cfile2 '_allD_locs.png'])
print(f9,'-dpng',[figp sname cfile2 '_all_locs.png'])



