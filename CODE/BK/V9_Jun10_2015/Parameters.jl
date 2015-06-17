#============== Parameters of the model =============#
#============ VARIABLE TYPE==========#
#! i.e. those things that change with time
type piscivore
	bio::Array{Any} # biomass g 
	tmet::Array{Float64} # temperature multiplier for metabolism
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
	Tp::Array{Any}
	Tb::Array{Any}
	Zm::Array{Any}
	Zl::Array{Any}
	dZm::Array{Any}
	dZl::Array{Any}
	det::Array{Any}
end

#============= PARAMETER TYPE ==========#
function make_parameters()

	#! Grid parameters
	GRD = load("./Data/Data_grid.jld"); # spatial information
    const global GRD_Z = GRD["Z"]
    const global GRD_A = GRD["AREA"][GRD["ID"]] # change this once rerun make_grid

	#! Integration parameters
	const global DT = 1.; # time step
    const global YEARS = 1; # integration period in years
    const global DAYS = 365; # number of days

	#! Amount of fishing
	const global FISHING = 80000000000 / (365/DT) # 80MT per year

	#! Number of size classes (#)
	const global PI_N = 10;
	const global PL_N = 10;
	const global DE_N = 10;

	#! Benthic-pelagic coupling cutoff (depth, m)
	const global PI_be_cutoff = 500

	#! Min body size (g)
	const global PI_smin = 10;
	const global PL_smin = .1;
	const global DE_smin = .1;

	#! Max body size (g)
	const global PI_smax = 10000;
	const global PL_smax = 1000;
	const global DE_smax = 1000;

	##! Body mass linearly distributed (g)
	const global PI_s = linspace((PI_smin),(PI_smax),PI_N)
	const global PL_s = linspace((PL_smin),(PL_smax),PL_N)
	const global DE_s = linspace((DE_smin),(DE_smax),DE_N)

	#! Median Zooplankton body mass 
	# 2.96 + 2.73ln(ESD) McCauly 1984, James Watkins, 
	# Lars Rudstam and Kristen Holeck later
	# Zm ESD = 1.1mm, Zl ESD = 11mm (from Charlie/COBALT)
	#const global Z_s = [exp(2.96 + (2.73*log(1.1)))*1.0e-6; 
	#					exp(2.96 + (2.73*log(11)))*1.0e-6]; 
	#const global Z_s = [exp(1.953 + (2.399*log(1.1)))*1.0e-6; 
	#					exp(1.953 + (2.399*log(11)))*1.0e-6]; 
	const global Z_s = [exp(1.953 + (2.399*log(2)))*1.0e-6; 
						exp(1.953 + (2.399*log(20)))*1.0e-6]; 

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
	ds_pi = abs(PI_s - (PI_s[end] * 0.1))
	pid = find(ds_pi .== minimum(ds_pi))[1]
	ds_pl = abs(PL_s - (PL_s[end] * 0.1))
	pld = find(ds_pl .== minimum(ds_pl))[1]
	ds_de = abs(DE_s - (DE_s[end] * 0.1))
	ded = find(ds_de .== minimum(ds_de))[1]
	const global PI_K = [ones(pid); linspace(1,0,PI_N - pid)]
	const global PL_K = [ones(pld); linspace(1,0,PL_N - pld)]
	const global DE_K = [ones(ded); linspace(1,0,DE_N - ded)]

	###! Metabolism constants (activity and basal)
	const global PI_act = exp(0.03*(3.9*PI_s.^0.13))
	const global PL_act = exp(0.03*(3.9*PL_s.^0.13))
	const global DE_act = exp(0.03*(3.9*DE_s.^0.13))
	const global PI_bas = 0.0033*PI_s.^-0.13
	const global PL_bas = 0.0033*PL_s.^-0.13
	const global DE_bas = 0.00033*DE_s.^-0.13 ### NOTE Changed to account for slow met

	###! Maximum search rate a as a function of body size
	# calculate swimming speed (m d-1)
	# factor in time spent swimming (low - Anieke's paper)
	# use body length to find area swept
	function fnc_a(s) # function to calc search rate from speed and length
		U = (3.9*s.^0.13) / 100 *60 *60 *24 #Speed: Megrey 2007
		L = (s/2.444).^(1/3.903) / 100 #Length: Anieke's supp 2008
		p = linspace(0.8,0.1,length(s));#Fraction of time: Anieke 2008 
		a = U .* (L.*2) .* p; #length x2 for visual diameter 
		return a
	end
	const global DE_a = fnc_a(DE_s)./10 # Anieke says detritivores move around less
	const global PI_a = fnc_a(PI_s)./1
	const global PL_a = fnc_a(PL_s)./1


	###! Background mortality
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
		#pij_de = fnc_pij(PI_s[i],DE_s,2,1) # prey preference for DE

		# normalize
		#Sj = [pij_pi; pij_pl; pij_de; pij_z];
		Sj = [pij_pi; pij_pl; pij_z];
		if sum(Sj)!=0.
			max = maximum(Sj);
			pij_pi = pij_pi / max;
			pij_pl = pij_pl / max;
			#pij_de = pij_de / max;
			pij_z  = pij_z  / max;
		end

		#! put in diet preference kernels
		PI_phi_PI[:,i] = pij_pi
		PI_phi_PL[:,i] = pij_pl
		#PI_phi_DE[:,i] = pij_de
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
		Sj = [pij_z];
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




