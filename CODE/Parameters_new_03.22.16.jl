#============== Parameters of the model =============#
#============ VARIABLE TYPE==========#
#! i.e. those things that change with time
type piscivore
	bio::Array{Any} # biomass g
	tmet::Array{Float64} # temperature multiplier for metabolism
	umet::Array{Any} # activity multiplier for metabolism
	tdif::Array{Float64} # fraction of time spent in the pelagic
	met::Array{Any} # metabolism g g-1 d-1
	I::Array{Any} # total ingestion g d-1
	I_z::Array{Any} # total zoo ingestion g d-1
	tau::Array{Any}  # handling time d g-1 g
	nu::Array{Any}  # total energy for growth g g-1 d-1
	gamma::Array{Any} # energy for somatic growth g g-1 d-1
	d::Array{Any} # total death g d-1
	REP::Array{Any} # total biomass reproduced g d-1
	GRW::Array{Any} # total biomass somatic growth g d-1
	MAT::Array{Any} # total biomass maturing g d-1
	MRT::Array{Any} # mortality rate (background and fishing potentially) g d-1
	enc_pi::Array{Any} # encounter rate with specific prey g g-1 d-1
	enc_pl::Array{Any} # encounter rate with specific prey g g-1 d-1
	enc_de::Array{Any} # encounter rate with specific prey g g-1 d-1
	enc_z::Array{Any} # encounter rate with specific prey g g-1 d-1
	ENC::Array{Any} # ecounter rate with total prey g g-1 d-1
end

type planktivore
	bio::Array{Any} # biomass g
	tmet::Array{Float64} # temperature multiplier for metabolism
	umet::Array{Any} # activity multiplier for metabolism
	met::Array{Any} # metabolism g g-1 d-1
	I::Array{Any} # total ingestion g d-1
	I_z::Array{Any} # total zoo ingestion g d-1
	tau::Array{Any}  # handling time d g-1 g
	nu::Array{Any}  # total energy for growth g g-1 d-1
	gamma::Array{Any} # energy for somatic growth g g-1 d-1
	d::Array{Any} # total death g d-1
	REP::Array{Any} # total biomass reproduced g d-1
	GRW::Array{Any} # total biomass somatic growth g d-1
	MAT::Array{Any} # total biomass maturing g d-1
	MRT::Array{Any} # mortality rate (background and fishing potentially) g d-1
	enc_z::Array{Any} # encounter rate with specific prey g g-1 d-1
	ENC::Array{Any} # ecounter rate with total prey g g-1 d-1
end

type detrivore
	bio::Array{Any} # biomass g
	tmet::Array{Float64} # temperature multiplier for metabolism
	umet::Array{Any} # activity multiplier for metabolism
	met::Array{Any} # metabolism g g-1 d-1
	I::Array{Any} # total ingestion g d-1
	tau::Array{Any}  # handling time d g-1 g
	nu::Array{Any}  # total energy for growth g g-1 d-1
	gamma::Array{Any} # energy for somatic growth g g-1 d-1
	d::Array{Any} # total death g d-1
	REP::Array{Any} # total biomass reproduced g d-1
	GRW::Array{Any} # total biomass somatic growth g d-1
	MAT::Array{Any} # total biomass maturing g d-1
	MRT::Array{Any} # mortality rate (background and fishing potentially) g d-1
	enc_de::Array{Any} # encounter rate with specific prey g g-1 d-1
	enc_be::Array{Any} #...
	ENC::Array{Any} # ecounter rate with total prey g g-1 d-1
end

type detritus
	bio::Array{Any} # biomass in detrital pool g m-2
	d::Array{Any} # biomass lost from predation g m-2 d-1
end

type environment
	Tp::Array{Any} # pelagic temperature
	Tb::Array{Any} # bottom temperature
	Zm::Array{Any} # medium zooplankton
	Zl::Array{Any} # large zoo
	dZm::Array{Any} # mortality rate of med zoo (COBALT)
	dZl::Array{Any} # mortality rate of large zoo (COBALT)
	det::Array{Any} # detrital flux
	K::Array{Any} # spawning flag
	T0::Array{Any} # baseline temp for spawning phenology
	Dthresh::Array{Any} # deg day threshold for spawning phenology
end

#============= PARAMETER TYPE ==========#
function make_parameters(harv)
	#! Input it switch for fishing (1 is yes, 0 is no)

	#! Grid parameters
	global GRD = load("./Data/Data_grid.jld"); # spatial information
  const global GRD_Z = GRD["Z"]
  const global GRD_A = GRD["AREA"]#[GRD["ID"]] # change this once rerun make_grid

	#! Integration parameters
	const global DT = 1.; # time step
  const global DAYS = 365; # number of days

	#! Amount of fishing
	if harv == 1
		const global FISHING = 80000000000 / (365/DT) # 80MT per year
	else
		const global FISHING = 0
	end

	#! Number of size classes (#) 3 = larvae, juveniles, adults
	const global PI_N = 3;
	const global PL_N = 2;
	const global DE_N = 3;

	#! Benthic-pelagic coupling cutoff (depth, m)
	const global PI_be_cutoff = 500

	#! Max length (L_inf)
	#! Size at maturation

	#! Min body size (mm)
	const global PI_lmin = 10^((log10(2)+log10(20))/2);
	const global PL_lmin = 10^((log10(2)+log10(20))/2);
	const global DE_lmin = 10^((log10(2)+log10(20))/2);

	#! Max body size (mm)
	const global PI_lmax = 10^((log10(200)+log10(2000))/2);
	const global PL_lmax = 10^((log10(20)+log10(200))/2);
	const global DE_lmax = 10^((log10(200)+log10(2000))/2);

	##! Length log10-linearly distributed (mm)
	#const global PI_l = linspace((PI_lmin),(PI_lmax),PI_N)
	#const global PL_l = linspace((PL_lmin),(PL_lmax),PL_N)
	#const global DE_l = linspace((DE_lmin),(DE_lmax),DE_N)
	const global PI_l = logspace(log10(PI_lmin),log10(PI_lmax),PI_N)
	const global PL_l = logspace(log10(PL_lmin),log10(PL_lmax),PL_N)
	const global DE_l = logspace(log10(DE_lmin),log10(DE_lmax),DE_N)

	##! Mass from length using Andersen & Beyer 2013
	# Convert from mm to cm and use their const coeff = 0.01g/cm3
	const global PI_s = 0.01 .* (0.1.*PI_l).^3;
	const global PL_s = 0.01 .* (0.1.*PL_l).^3;
	const global DE_s = 0.01 .* (0.1.*DE_l).^3;

	#! Median Zooplankton body mass
	# James Watkins, Lars Rudstam and Kristen Holeck
	# (from Charlie/COBALT)
	# This is in dry weight, need to convert to wet weight 0.2g dry = 1 g wet
	const Zm_ESD = 10^((log10(0.2)+log10(2))/2);
	const Zl_ESD = 10^((log10(2)+log10(20))/2);

	const global Z_s = [exp(1.953 + (2.399*log(Zm_ESD)))*1.0e-6;
						exp(1.953 + (2.399*log(Zl_ESD)))*1.0e-6];

	#! Log body mass  (log10 g)
	const global DE_slog = log10(DE_s)
	const global PI_slog = log10(PI_s)
	const global PL_slog = log10(PL_s)
	const global Z_slog  = log10(Z_s)

	#! Ratio of initial and final body sizes per size-class
	function fnc_z(s)
		ds = diff(s)/2;
		ds = [ds[1];ds]
		z = zeros(length(s))
		for i = 1:length(s)
			s_d = s[i] - ds[i] # initial body size
			s_u = s[i] + ds[i] # final body size
			z[i] = (10^s_d) / (10^s_u) # ratio of initial to final
		end
		return z
	end
	const global DE_z = fnc_z(DE_slog);
	const global PI_z = fnc_z(PI_slog);
	const global PL_z = fnc_z(PL_slog);

	###! Assimilation efficiency lambdaa
	const global DE_lambda = ones(DE_N) .* 0.7;
	const global PI_lambda = ones(PI_N) .* 0.7;
	const global PL_lambda = ones(PL_N) .* 0.7;

	###! Kappa rule K as a function of body size
	# K = fraction of energy consumed diverted to somatic growth
	# Maturity fn from Andersen & Beyer 2013, their n=0.75, nm=0.25
	# My nm from FB and RAM data nm=0.1213
	const nm = 0.1213;
	const n = 0.75;
	PI_psi = (1+(PI_s/(nm*PI_s[end]))^-10)^-1 * (PI_s/PI_s[end])^(1-n);
	PL_psi = (1+(PL_s/(nm*PL_s[end]))^-10)^-1 * (PL_s/PL_s[end])^(1-n);
	DE_psi = (1+(DE_s/(nm*DE_s[end]))^-10)^-1 * (DE_s/DE_s[end])^(1-n);
	const global PI_K = 1-PI_psi;
	const global PL_K = 1-PL_psi;
	const global DE_K = 1-DE_psi;

	###! Metabolism constants (activity and basal)
	const global PI_act = exp(0.03*(3.9*PI_s.^0.13))
	const global PL_act = exp(0.03*(3.9*PL_s.^0.13))
	const global DE_act = exp(0.03*(3.9*DE_s.^0.13))
	const global PI_bas = 0.0033*PI_s.^-0.13
	const global PL_bas = 0.0033*PL_s.^-0.13
	const global DE_bas = 0.0033*DE_s.^-0.13 ### NOTE Changed to account for slow met

	###! Swimming speed (Megrey 2007)
	# Why is swim speed a fraction from 0.8 to 0.1 as gets bigger?
	# Time spend swimming?
	const global PI_U = ((3.9*PI_s.^0.13)/100*60*60*24) .* linspace(0.8,0.1,length(PI_s))
	const global PL_U = ((3.9*PL_s.^0.13)/100*60*60*24) .* linspace(0.8,0.1,length(PL_s))
	const global DE_U = ((3.9*DE_s.^0.13)/100*60*60*24) .* linspace(0.8,0.1,length(DE_s))
  # With temp-dep
	#const global PI_U = ((3.9*PI_s.^0.13 * exp(0.149*T)) /100*60*60*24)
	#const global PL_U = ((3.9*PL_s.^0.13 * exp(0.149*T)) /100*60*60*24)
	#const global DE_U = ((3.9*DE_s.^0.13 * exp(0.149*T)) /100*60*60*24)


	###! Maximum search rate a as a function of body size
	# calculate swimming speed (m d-1)
	# factor in time spent swimming (low - Anieke's paper)
	# use body length to find area swept
	function fnc_a(s) # function to calc search rate from speed and length
		U = (3.9*s.^0.13) / 100 *60 *60 *24 #Speed: Megrey 2007
		#L = (s/2.444).^(1/3.903) / 100 #Length: Anieke's supp 2008
		L = 10*(s/0.01)^(1/3) #Length from Andersen & Beyer 2013
		p = linspace(0.8,0.1,length(s));#Fraction of time: Anieke 2008
		a = U .* (L.*2) .* p; #length x2 for visual diameter
		return a
	end
	const global DE_a = fnc_a(DE_s)./1 # Anieke says detritivores move around less
	const global PI_a = fnc_a(PI_s)./1
	const global PL_a = fnc_a(PL_s)./1


	###! Background mortality
	#Currently increases from 0 to 0.01
	#Megrey et al =0.44/yr
	#Andersen & Beyer 2013 = 0.35 * 4.5 * s^(-0.25) (includes predation, excludes fishing)
	ds_pi = abs(PI_s - (PI_s[end] * 0.1))
	pid = find(ds_pi .== minimum(ds_pi))[1]
	ds_pl = abs(PL_s - (PL_s[end] * 0.1))
	pld = find(ds_pl .== minimum(ds_pl))[1]
	ds_de = abs(DE_s - (DE_s[end] * 0.1))
	ded = find(ds_de .== minimum(ds_de))[1]
	const global PI_mrt = [zeros(pid); linspace(0,0.01,PI_N - pid)]
	const global PL_mrt = [zeros(pld); linspace(0,0.01,PL_N - pld)]
	const global DE_mrt = [zeros(ded); linspace(0,0.01,DE_N - ded)]

	###! Diet Preference Phi (j = prey, i = pred)
	#for each fish, get all prey and calculate preference
	#normalize so that most preferred prey preference == 1
	#distribute information into feeding kernels
	function fnc_pij(w1,w2,beta,sigma) # Blanchard function to calculate diet preference
		# w1: body size of pred (g)
		# w2: body size of prey (g)
		#beta = 2; # mean PPMR
		#sigma = 1; # sd PPMR
		PIJ = zeros(length(w2))

		for i = 1:length(w2)
			x1 = log10(w1)
			x2 = log10(w2[i])
			if w1 > w2[i]
				pij = exp(-((x1-x2-beta).^2) / (2.*sigma.^2));
			else
				pij = 0.
			end
			PIJ[i] = pij;
		end
		return PIJ
	end

	####! Piscivore
	const global PI_phi_PI = zeros(PI_N,PI_N);
	const global PI_phi_PL = zeros(PL_N,PI_N);
	const global PI_phi_DE = zeros(DE_N,PI_N);
	const global PI_phi_Z  = zeros(2,PI_N);
	for i = 1:PI_N
		pij_pi = fnc_pij(PI_s[i],PI_s,2,1) # prey preference for PI
		pij_pl = fnc_pij(PI_s[i],PL_s,2,1) # prey preference for PL
		pij_z  = fnc_pij(PI_s[i],Z_s,2,1) # prey preference for Zoo
		pij_de = fnc_pij(PI_s[i],DE_s,2,1) # prey preference for DE

		# normalize
		Sj = [pij_pi; pij_pl; pij_de; pij_z];
		#Sj = [pij_pi; pij_pl; pij_z];
		if sum(Sj)!=0.
			max = maximum(Sj);
			pij_pi = pij_pi / max;
			pij_pl = pij_pl / max;
			pij_de = pij_de / max;
			pij_z  = pij_z  / max;
		end

		#! put in diet preference kernels
		PI_phi_PI[:,i] = pij_pi
		PI_phi_PL[:,i] = pij_pl
		PI_phi_DE[:,i] = pij_de
		PI_phi_Z[:,i]  = pij_z
	end
	#! detritus and zoo feeding preference
	const global PI_phi_BE  = zeros(1,PI_N)
	#PI_phi_Z = zeros(2,PI_N) # Forcing piscivore to not eat zoo
	#PI_phi_Z[:,2:end] = zeros(2,PI_N-1)
	#PI_phi_Z[:,1] = 1.; ## Forcing piscivore smallest to eat zoo
	#PI_phi_DE  = zeros(DE_N,PI_N)
	#PI_phi_Z  = zeros(1,PI_N)


	####! Planktivore (eats only zooplankton)
	const global PL_phi_Z = zeros(2,PL_N);
	for i = 1:PL_N
		pij_z = fnc_pij(PL_s[i],Z_s,3,1) # prey preference for PL

		# normalize
		Sj = collect(pij_z); #use Sj = collect(pij_z) instead of Sj = [pij_z]
		if sum(Sj)!=0.
			max = maximum(Sj);
			pij_z = pij_z / max;
		end

		#! put in diet preference kernels
		PL_phi_Z[:,i]  = pij_z
	end
	#! other feeding kernels
	const global PL_phi_BE = zeros(1,PL_N)
	const global PL_phi_PI = zeros(PI_N,PL_N)
	const global PL_phi_PL = zeros(PL_N,PL_N)
	const global PL_phi_DE = zeros(DE_N,PL_N)


	####! Detritivore (eats only detritus and other detritivores)
	const global DE_phi_DE = zeros(DE_N,DE_N);
	const global DE_phi_BE = zeros(1,DE_N);
	for i = 1:DE_N
		pij_de = fnc_pij(DE_s[i],DE_s,2,1) # prey preference for DE
		pij_be = fnc_pij(DE_s[i],DE_s[1]/100,2,1) # prey preference for DE

		# normalize
		Sj = [pij_de; pij_be];
		if sum(Sj)!=0.
			max = maximum(Sj);
			pij_de = pij_de / max;
			pij_be = pij_be / max;
		end

		#! put in diet preference kernels
		DE_phi_DE[:,i]  = pij_de
		DE_phi_BE[:,i]  = pij_be
	end
	#! other feeding kernels
	const global DE_phi_PI = zeros(PI_N,DE_N)
	const global DE_phi_PL = zeros(PL_N,DE_N)
	const global DE_phi_Z  = zeros(2,DE_N)

end
