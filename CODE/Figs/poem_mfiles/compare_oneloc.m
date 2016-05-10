% Compare continuous vs spawning season
% Pristine historical at one location
% 145 years

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

load('Oneloc_pris_phenol_all.mat')
load('Oneloc_pris_all.mat')
load('cmap_ppt_angles.mat')

%% Plots over time
x=1:length(sp);
yfrac=x/365;
y=1861:(1/365):(1861+yfrac(end));
%y=1861:(1/365):2005;
y=y(1:(end-1));
lstd=length(sp);
id1 = 0:365:(lstd-1);
id2 = 365:365:(lstd);
ID  = [id1 id2];


%% Final mean biomass size spectrum
t=1:length(sp);
%lyr=t((end-365):end);
lyr=t((end-(30*365)+1):end);
sp_sum(1,:)=sum(sp(lyr,:));
sf_sum(1,:)=sum(sf(lyr,:));
sd_sum(1,:)=sum(sd(lyr,:));
mp_sum(1,:)=sum(mp(lyr,:));
mf_sum(1,:)=sum(mf(lyr,:));
md_sum(1,:)=sum(md(lyr,:));
lp_sum(1,:)=sum(lp(lyr,:));

sp_mean(1,:)=mean(sp(lyr,:));
sf_mean(1,:)=mean(sf(lyr,:));
sd_mean(1,:)=mean(sd(lyr,:));
mp_mean(1,:)=mean(mp(lyr,:));
mf_mean(1,:)=mean(mf(lyr,:));
md_mean(1,:)=mean(md(lyr,:));
lp_mean(1,:)=mean(lp(lyr,:));

sp_sum(2,:)=sum(psp(lyr,:));
sf_sum(2,:)=sum(psf(lyr,:));
sd_sum(2,:)=sum(psd(lyr,:));
mp_sum(2,:)=sum(pmp(lyr,:));
mf_sum(2,:)=sum(pmf(lyr,:));
md_sum(2,:)=sum(pmd(lyr,:));
lp_sum(2,:)=sum(plp(lyr,:));

sp_mean(2,:)=mean(psp(lyr,:));
sf_mean(2,:)=mean(psf(lyr,:));
sd_mean(2,:)=mean(psd(lyr,:));
mp_mean(2,:)=mean(pmp(lyr,:));
mf_mean(2,:)=mean(pmf(lyr,:));
md_mean(2,:)=mean(pmd(lyr,:));
lp_mean(2,:)=mean(plp(lyr,:));

%%
for n=1:length(spots)
    loc=spots{n};
    P_mean=[sp_mean(:,n),mp_mean(:,n),lp_mean(:,n)];
    F_mean=[sf_mean(:,n),mf_mean(:,n)];
    D_mean=[sd_mean(:,n),md_mean(:,n)];
    
    figure(n)
    subplot(3,1,3)
    bar(P_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 4])
    xlabel('Stage')
    title('Pelagic Piscivores')
    ylabel('Mean Biomass (g km^-^2)')
    
    subplot(3,1,1)
    bar(F_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 3])
    title({loc; 'Forage Fishes'})
    ylabel('Mean Biomass (g km^-^2)')
    legend('Constant','Phenology')
    
    subplot(3,1,2)
    bar(D_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 3])
    title('Demersal Piscivores')
    ylabel('Mean Biomass (g km^-^2)')
    
    print('-dpng',[fpath loc '_oneloc_compare_biomass_spec.png'])
    
    
    f20=figure(20);
    subplot(2,3,n)
    bar(P_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 4])
    set(gca,'XTickLabel',{'L','J','A'})
    xlabel('Stage')
    title(loc)
    if (n==1)
    ylabel('Mean Biomass (g km^-^2)')
    end
    if (n==4)
    ylabel('Pelagic Piscivores')
    end
    
    f21=figure(21);
    subplot(2,3,n)
    bar(F_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 3])
    set(gca,'XTickLabel',{'I','A'})
    xlabel('Stage')
    title(loc)
    if (n==1)
    ylabel('Mean Biomass (g km^-^2)')
    end
    if (n==4)
    ylabel('Forage Fishes')
    end
    
    f22=figure(21);
    subplot(2,3,n)
    bar(D_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 3])
    set(gca,'XTickLabel',{'I','A'})
    xlabel('Stage')
    title(loc)
    if (n==1)
    ylabel('Mean Biomass (g km^-^2)')
    end
    if (n==4)
    ylabel('Demersal Piscivores')
    end
end
print(f20,'-dpng',[fpath 'Pisc_oneloc_compare_biomass.png'])
print(f21,'-dpng',[fpath 'Forage_oneloc_compare_biomass.png'])
print(f22,'-dpng',[fpath 'Dem_oneloc_compare_biomass.png'])

%% test color order
% figure
% plot(1:5,ones(5,1),'color',cmap_ppt(1,:),'LineWidth',2); hold on;
% plot(1:5,2*ones(5,1),'color',cmap_ppt(2,:),'LineWidth',2); hold on;
% plot(1:5,3*ones(5,1),'color',cmap_ppt(3,:),'LineWidth',2); hold on;
% plot(1:5,4*ones(5,1),'color',cmap_ppt(4,:),'LineWidth',2); hold on;
% plot(1:5,5*ones(5,1),'color',cmap_ppt(5,:),'LineWidth',2); hold on;
% plot(1:5,6*ones(5,1),'color',cmap_ppt(6,:),'LineWidth',2); hold on;
% plot(1:5,7*ones(5,1),'color',cmap_ppt(7,:),'LineWidth',2); hold on;
% ylim([0 8])

%% Recruitment
yr=t((end-(30*365)+1):end);
SPL=Mlp(yr,:);
SFL=Mmf(yr,:);
SDL=Mmd(yr,:);
PA=lp(yr,:);
FA=mf(yr,:);
DA=md(yr,:);

pSPL=pMlp(yr,:);
pSFL=pMmf(yr,:);
pSDL=pMmd(yr,:);
pPA=plp(yr,:);
pFA=pmf(yr,:);
pDA=pmd(yr,:);

st=1:365:length(yr);
en=365:365:length(yr);
PLy = NaN*ones(30,6);
FLy = PLy;
DLy = PLy;
PAy = PLy;
FAy = PLy;
DAy = PLy;
pPLy = PLy;
pFLy = PLy;
pDLy = PLy;
pPAy = PLy;
pFAy = PLy;
pDAy = PLy;
for n=1:30
    PLy(n,:) = nansum(PL(st(n):en(n),:));
    FLy(n,:) = nansum(FL(st(n):en(n),:));
    DLy(n,:) = nansum(DL(st(n):en(n),:));
    PAy(n,:) = nansum(PA(st(n):en(n),:));
    FAy(n,:) = nansum(FA(st(n):en(n),:));
    DAy(n,:) = nansum(DA(st(n):en(n),:));
    pPLy(n,:) = nansum(pPL(st(n):en(n),:));
    pFLy(n,:) = nansum(pFL(st(n):en(n),:));
    pDLy(n,:) = nansum(pDL(st(n):en(n),:));
    pPAy(n,:) = nansum(pPA(st(n):en(n),:));
    pFAy(n,:) = nansum(pFA(st(n):en(n),:));
    pDAy(n,:) = nansum(pDA(st(n):en(n),:));
end

%%
close all
for n=1:length(spots)
    loc=spots{n};
    
    figure(n)
    subplot(3,1,1)
    plot(1976:2005,PLy(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pPLy(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    title({loc; 'PP Recruitment'})
    ylabel('Biomass (g km^-^2)')
    %text(2006,max(PLy(:))+10,loc,'FontWeight','bold')
    
    subplot(3,1,2)
    plot(1976:2005,FLy(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pFLy(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    legend('Constant','Phenology')
    title('FF Recruitment')
    ylabel('Biomass (g km^-^2)')
    
    subplot(3,1,3)
    plot(1976:2005,DLy(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pDLy(:,n),'Linewidth',2); hold on;
    title('DP Recruitment')
    xlim([1976 2005])
    ylabel('Biomass (g km^-^2)')
    print('-dpng',[fpath loc '_oneloc_compare_recruitment.png'])
    
    f10=figure(10);
    subplot(2,3,n)
    plot(1976:2005,FLy(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pFLy(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    ylabel('Forage fish recruits (g km^-^2)')
    title(loc)
    
    f11=figure(11);
    subplot(2,3,n)
    plot(1976:2005,PLy(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pPLy(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    title(loc)
    ylabel('Pelagic piscivore recruits (g km^-^2)')
    
    f12=figure(12);
    subplot(2,3,n)
    plot(1976:2005,DLy(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pDLy(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    title(loc)
    ylabel('Demersal piscivore recruits (g km^-^2)')
end
print(f10,'-dpng',[fpath 'Forage_oneloc_compare_recruitment.png'])
print(f11,'-dpng',[fpath 'Pisc_oneloc_compare_recruitment.png'])
print(f12,'-dpng',[fpath 'Dem_oneloc_compare_recruitment.png'])

%% EVERYTHING SCALED TO MAX BIOMASS
Ptot = sp_mean+mp_mean+lp_mean;
Ftot = sf_mean+mf_mean;
Dtot = sd_mean+md_mean;
Tot  = Ptot+Ftot+Dtot;

% Final mean biomass size spectrum

for n=1:length(spots)
    loc=spots{n};
    P_mean=[sp_mean(:,n),mp_mean(:,n),lp_mean(:,n)];
    F_mean=[sf_mean(:,n),mf_mean(:,n)];
    P_mean = P_mean ./ repmat(Tot(:,n),1,3);
    F_mean = F_mean ./ repmat(Tot(:,n),1,2);
    
    figure(n)
    subplot(3,1,3)
    bar(P_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 4])
    xlabel('Stage')
    title('Pelagic Piscivores')
    ylabel('Fraction of Total Mean Biomass')
    
    subplot(3,1,1)
    bar(F_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 3])
    title({loc; 'Forage Fishes'})
    ylabel('Fraction of Total Mean Biomass')
    legend('Constant','Phenology')
    
    subplot(3,1,2)
    bar(D_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 3])
    xlabel('Stage')
    title('Demersal Piscivores')
    ylabel('Fraction of Total Mean Biomass')
    legend('Constant','Phenology')
    
    print('-dpng',[fpath loc '_oneloc_compare_biomass_spec_scale.png'])
    
    
    f20=figure(20);
    subplot(2,3,n)
    bar(P_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 4])
    set(gca,'XTickLabel',{'L','J','A'})
    xlabel('Stage')
    title(loc)
    if (n==1)
    ylabel('Fraction of Total Mean Biomass')
    end
    if (n==4)
    ylabel('Pelagic Piscivore')
    end
    if (n==3)
    legend('C','P')
    legend('location','northeast')
    end
    
    f21=figure(21);
    subplot(2,3,n)
    bar(F_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 3])
    set(gca,'XTickLabel',{'I','A'})
    xlabel('Stage')
    title(loc)
    if (n==1)
    ylabel('Fraction of Total Mean Biomass')
    end
    if (n==4)
    ylabel('Forage Fish')
    end
    if (n==3)
    legend('C','P')
    legend('location','northwest')
    end
    
    f22=figure(22);
    subplot(2,3,n)
    bar(D_mean')
    colormap(cmap_ppt(1:3,:))
    xlim([0 3])
    set(gca,'XTickLabel',{'I','A'})
    xlabel('Stage')
    title(loc)
    if (n==1)
    ylabel('Fraction of Total Mean Biomass')
    end
    if (n==4)
    ylabel('Demersal Piscivores')
    end
    if (n==3)
    legend('C','P')
    legend('location','northwest')
    end
end
print(f20,'-dpng',[fpath 'Pisc_oneloc_compare_biomass_scale.png'])
print(f21,'-dpng',[fpath 'Forage_oneloc_compare_biomass_scale.png'])
print(f22,'-dpng',[fpath 'Dem_oneloc_compare_biomass_scale.png'])

%% Recruitment scaled to max recruitment
RPmax(1,:) = max(PLy);
RPmax(2,:) = max(pPLy);
RFmax(1,:) = max(FLy);
RFmax(2,:) = max(pFLy);
RDmax(1,:) = max(DLy);
RDmax(2,:) = max(pDLy);

RPsc  = PLy ./ repmat(RPmax(1,:),30,1);
pRPsc = pPLy ./ repmat(RPmax(2,:),30,1);
RFsc  = FLy ./ repmat(RFmax(1,:),30,1);
pRFsc = pFLy ./ repmat(RFmax(2,:),30,1);
RDsc  = DLy ./ repmat(RDmax(1,:),30,1);
pRDsc = pDLy ./ repmat(RDmax(2,:),30,1);

close all
for n=1:length(spots)
    loc=spots{n};
    
    figure(n)
    subplot(3,1,1)
    plot(1976:2005,RFsc(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pRFsc(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    ylabel('Scaled Recruits')
    title({loc; 'Forage fishes'})
    
    subplot(3,1,3)
    plot(1976:2005,RPsc(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pRPsc(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    ylabel('Scaled Recruits')
    title('Pelagic piscivores')
    
    subplot(3,1,2)
    plot(1976:2005,RDsc(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pRDsc(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    ylabel('Scaled Recruits')
    title('Demersal piscivores')
    
    print('-dpng',[fpath loc '_oneloc_compare_recruitment_scaled.png'])
    
    
    f10=figure(10);
    subplot(2,3,n)
    plot(1976:2005,RFsc(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pRFsc(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    ylabel('Forage fish scaled recruits')
    if (n==1)
    legend('C','P')
    legend('location','southeast')
    end
    title(loc)
    
    f11=figure(11);
    subplot(2,3,n)
    plot(1976:2005,RPsc(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pRPsc(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    title(loc)
    if (n==3)
    legend('C','P')
    legend('location','southwest')
    end
    ylabel('Pelagic piscivore scaled recruits')
    
    f12=figure(12);
    subplot(2,3,n)
    plot(1976:2005,RDsc(:,n),'Linewidth',2); hold on;
    plot(1976:2005,pRDsc(:,n),'Linewidth',2); hold on;
    xlim([1976 2005])
    title(loc)
    if (n==3)
    legend('C','P')
    legend('location','southwest')
    end
    ylabel('Demersal piscivore scaled recruits')
end
print(f10,'-dpng',[fpath 'Forage_oneloc_compare_recruitment_scaled.png'])
print(f11,'-dpng',[fpath 'Pisc_oneloc_compare_recruitment_scaled.png'])
print(f12,'-dpng',[fpath 'Dem_oneloc_compare_recruitment_scaled.png'])
    

