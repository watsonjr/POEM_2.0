% Visualize output of POEM
% Pristine historical at one location
% 145 years

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

spots = {'GB','EBS','OSP','HOT','BATS','NS'};

%% Continuous
sp=NaN*ones(145*365,length(spots));
sf = sp;
sd = sp;
mf = sp;
md = sp;
mp = sp;
lp = sp;
Rmf = sp;
Rmd = sp;
Rlp = sp;

for s=1:length(spots)
    close all
    loc = spots{s};
    sname1 = 'Oneloc_pris_' ;
    sname2 = [loc '_'];
    SP = csvread([dpath sname1 sname2 'Sml_p.csv']);
    SF = csvread([dpath sname1 sname2 'Sml_f.csv']);
    SD = csvread([dpath sname1 sname2 'Sml_d.csv']);
    MP = csvread([dpath sname1 sname2 'Med_p.csv']);
    MF = csvread([dpath sname1 sname2 'Med_f.csv']);
    MD = csvread([dpath sname1 sname2 'Med_d.csv']);
    LP = csvread([dpath sname1 sname2 'Lrg_p.csv']);
    
    repMF = csvread([dpath sname1 'Rep_' sname2 'Med_f.csv']);
    repMD = csvread([dpath sname1 'Rep_' sname2 'Med_d.csv']);
    repLP = csvread([dpath sname1 'Rep_' sname2 'Lrg_p.csv']);
    
    sp(:,s) = SP;
    sf(:,s) = SF;
    sd(:,s) = SD;
    mp(:,s) = MP;
    mf(:,s) = MF;
    md(:,s) = MD;
    lp(:,s) = LP;
    
    Rmf(:,s) = repMF;
    Rmd(:,s) = repMD;
    Rlp(:,s) = repLP;
    
end

save('Oneloc_pris_all.mat','sp','sf','sd','mp','mf','md','lp',...
    'Rmf','Rmd','Rlp','spots')

%% Phenology
psp=NaN*ones(145*365,length(spots));
psf = psp;
psd = psp;
pmf = psp;
pmd = psp;
pmp = psp;
plp = psp;
pDDmf = psp;
pDDmd = psp;
pDDlp = psp;
pKmf = psp;
pKmd = psp;
pKlp = psp;
pRmf = psp;
pRmd = psp;
pRlp = psp;
for s=1:length(spots)
    close all
    loc = spots{s};
    sname1 = 'Oneloc_pris_';
    sname2 = [loc '_phenol_'];
    SP = csvread([dpath sname1 sname2 'Sml_p.csv']);
    SF = csvread([dpath sname1 sname2 'Sml_f.csv']);
    SD = csvread([dpath sname1 sname2 'Sml_d.csv']);
    MP = csvread([dpath sname1 sname2 'Med_p.csv']);
    MF = csvread([dpath sname1 sname2 'Med_f.csv']);
    MD = csvread([dpath sname1 sname2 'Med_d.csv']);
    LP = csvread([dpath sname1 sname2 'Lrg_p.csv']);
    
    ddMF = csvread([dpath sname1 'DD_' sname2 'Med_f.csv']);
    ddMD = csvread([dpath sname1 'DD_' sname2 'Med_d.csv']);
    ddLP = csvread([dpath sname1 'DD_' sname2 'Lrg_p.csv']);
    
    kMF = csvread([dpath sname1 'K_' sname2 'Med_f.csv']);
    kMD = csvread([dpath sname1 'K_' sname2 'Med_d.csv']);
    kLP = csvread([dpath sname1 'K_' sname2 'Lrg_p.csv']);
    
    repMF = csvread([dpath sname1 'Rep_' sname2 'Med_f.csv']);
    repMD = csvread([dpath sname1 'Rep_' sname2 'Med_d.csv']);
    repLP = csvread([dpath sname1 'Rep_' sname2 'Lrg_p.csv']);
    
    psp(:,s) = SP;
    psf(:,s) = SF;
    psd(:,s) = SD;
    pmp(:,s) = MP;
    pmf(:,s) = MF;
    pmd(:,s) = MD;
    plp(:,s) = LP;
    
    pDDmf(:,s) = ddMF;
    pDDmd(:,s) = ddMD;
    pDDlp(:,s) = ddLP;
    
    pKmf(:,s) = kMF;
    pKmd(:,s) = kMD;
    pKlp(:,s) = kLP;
    
    pRmf(:,s) = repMF;
    pRmd(:,s) = repMD;
    pRlp(:,s) = repLP;
end

save('Oneloc_pris_phenol_all.mat','psp','psf','psd','pmp','pmf',...
    'pmd','plp','pDDmf','pDDmd','pDDlp','pKmf','pKmd','pKlp',...
    'pRmf','pRmd','pRlp','spots')
