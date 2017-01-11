%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function make_parameters(harv,frate)
%     global FISHING 
%     global LFISHING 
%     global MFISHING 
    global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl 
    global Z_s Z_m Z_l Lambda K_l K_j K_a fcrit k h gamma q 
    global bent_eff rfrac Tu_s Tu_m Tu_l Nat_mrt MORT 
    global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE 
    global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE 

    %! Input switch for fishing (1 is yes, 0 is no)

    %! Integration parameters
    DT = 1.0; % time step

    %! Amount of fishing
    % if harv == 1
    % 	%FISHING = 80000000000 / (365/DT); % 80MT per year
    % 	LFISHING = frate/365.0;
    % 	MFISHING = LFISHING;
    % 	%MFISHING = 0.5 * LFISHING;
    % 	%MFISHING = 1.5/365.0;
    % 	%LFISHING = 0.5 * MFISHING;
    % else
    % 	MFISHING = 0;
    % 	LFISHING = 0;
    % end

    %! Benthic-pelagic coupling cutoff (depth, m)
    PI_be_cutoff = 200;
    % 0:no coupling; 1:demersal coupled only; 2:pelagic & demersal coupled
    pdc = 1;

    %! body lengths (mm)
    L_s = 10^((log10(2)+log10(20))/2); % small
    L_m = 10^((log10(20)+log10(200))/2); % medium
    L_l = 10^((log10(200)+log10(2000))/2); % large

    %%! Mass from length using Andersen & Beyer 2013
    % Convert from mm to cm and use their const coeff = 0.01g/cm3
    M_s = 0.01 * (0.1*L_s)^3;
    M_m = 0.01 * (0.1*L_m)^3;
    M_l = 0.01 * (0.1*L_l)^3;

    %! Mediant Zooplankton size in mm
    % from Charlie/COBALT
    %! Median Zooplankton body mass in g wet weight
    % eq from James Watkins, Lars Rudstam and Kristen Holeck in dry weight
    % convert to wet weight with 1g dry = 9 g wet
    L_zm = 10^((log10(0.2)+log10(2))/2); % lengths (ESD)
    L_zl = 10^((log10(2)+log10(20))/2);
    M_zm = 9.0 * exp(1.953 + (2.399*log(L_zm)))*1.0e-6; % body mass
    M_zl = 9.0 * exp(1.953 + (2.399*log(L_zl)))*1.0e-6;

    %! Ratio of initial and final body sizes per size-class
    Z_s = (0.01*(0.1*2)^3) / (0.01*(0.1*20)^3);	 %M_s./M_m
    Z_m = (0.01*(0.1*20)^3) / (0.01*(0.1*200)^3);	 %M_m./M_l
    Z_l = (0.01*(0.1*200)^3) / (0.01*(0.1*2000)^3); %NA

    %%%! Assimilation efficiency lambda (constant across everything)
    Lambda = 0.7;

    %%%! Kappa rule K as a function of body size
    % K = fraction of energy consumed diverted to somatic growth
    K_l = 1;
    K_j = 1;
    K_a = 0;

    %%%! Metabolism constants (activity and basal)
    fcrit = 0.40;	% feeding level needed to meet resting metabolic demands; 0.05-0.2
    k = 4.8;		% 10 g^(1-p)/yr at 10C; 4.8 at 10C NS mizer

    %%%! Consumption constants
    h = 60.0;  		% h=85 g^(0.25)/yr at 10C in Cmax eq; h=40 at 10C NS; h=60 at 15C?
    % tune so Cobs/Cmax ~ 0.6
    gamma = 2.9e3;     % m^3 g^(?q)/year at 10C; equiv to Andersen, Hartvig gamma = 0.8e4; mizer = 2.9e3?
    q = 0.8; 			% q=0.75-1 in beta eq in consumption

    %%%! Transfer efficiency of detritus to benthic prey
    bent_eff = 0.05;

    %%%! Reproductive efficiency
    rfrac = 0.08;

    %! Fraction of time spent swimming (from Van Leeuwen)
    Tu_s = 1.0;
    Tu_m = 1.0; %0.5
    Tu_l = 1.0; %0.1

    %%%! Background mortality
    %Currently increases from 0 to 0.01
    %Megrey et al =0.44/yr
    %Andersen & Beyer 2013 = 0.35 * 4.5 * s^(-0.25) (includes predation, excludes fishing)
    Nat_mrt = 0.0; %0.44 / 365
    %0=none, 1=constant, 2=temp-dep, 3=large only, 4=large temp-dep
    MORT = 0;

    %%%! Diet Preference Phi
    % The predator prey mass ratio is assumed 3 orders of mag, i.e. 1000, i.e. one step down
    % We don't have a pred-prey matrix anymore, we are instead explicit about who eats who:
    %-----
    %small forage fish eats medium zoo
    %small piscivores eats medium zoo
    %small detritivore eats medium zoo
    %medium forage fish eats medium & large zoo, all small fishes
    %medium piscivore eats medium & large zoo, all small fishes
    %medium detritivore eats detritus
    %large piscivore eats medium forage fish, medium piscivore, medium detritivore
    %large detritivore eats detritus, medium forage fish, medium piscivore, medium detrivore

    MF_phi_MZ = 0.1;
    MF_phi_LZ = 1.0;
    MF_phi_S = 1.0;
    MP_phi_MZ = 0.1;
    MP_phi_LZ = 1.0;
    MP_phi_S = 1.0;
    MD_phi_BE = 1.0;
    LP_phi_MF = 1.0;
    LP_phi_MP = 1.0;
    LP_phi_MD = 1.0;
    LD_phi_MF = 1.0;
    LD_phi_MP = 1.0;
    LD_phi_MD = 1.0;
    LD_phi_BE = 1.0;
    
%-----
end
