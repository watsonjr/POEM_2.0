#============== Parameters of the model =============#
#============ VARIABLE TYPE==========#
#! i.e. those things that change with time
type fish
	BIOMASS::Array{Float64} # biomass density
	MET::Array{Float64} # metabolic costs
	TAU::Array{Float64} # handling times
end

#============= PARAMETER TYPE ==========#
#Note: handling times and metabolic rates are variables...
#tau::Array{Float64} # handling time (d g_j^-1 g_i)
#I::Array{Float64} # metabolic costs (g_i g_j^-1 d^-1)
type PARAM
	N::Int # Number of size classes
	s_min::Float64 # minimum size (g)
	s_max::Float64 # maximum size (g)
	s::Array{Float64} # median size (g) of earch size class
	slog::Array{Float64} # log10 median size (g) of earch size class
	z::Array{Float64} # ratio of initial and final body sizes 
	lambda::Array{Float64} # assimilation efficiency
	K::Array{Float64} # kappa, fraction of energy given to somatic growth
	a::Array{Float64} # search rate (m2 d^-1 g_i^-1)
	Phi_Z::Array{Float64} # preference for eating zooplankton
	Phi_PI::Array{Float64} # preference for eating piscivore
	Phi_PL::Array{Float64} # preference for eating planktivore
	Phi_DE::Array{Float64} # preference for eating detritivore
	Phi_W::Array{Float64} # preference for eating detritus
end


function make_parameters()

	#! Number of size classes (#)
	DE_N = 5;
	PI_N = 10;
	PL_N = 10;

	#! Min body size (g)
	DE_smin = 1;
	PI_smin = 1;
	PL_smin = 1;

	#! Max body size (g)
	DE_smax = 1000;
	PI_smax = 200;
	PL_smax = 200;

	#! Body mass loglinearly distributed (g)
	DE_s = 10.^(linspace(log10(DE_smin),log10(DE_smax),DE_N));
	PI_s = 10.^(linspace(log10(PI_smin),log10(PI_smax),PI_N));
	PL_s = 10.^(linspace(log10(PL_smin),log10(PL_smax),PL_N));
	Z_s = [128e-6; 70000e-6];

	#! Log body mass  (log10 g)
	DE_slog = log10(DE_s)
	PI_slog = log10(PI_s)
	PL_slog = log10(PL_s)
	Z_slog = log10(Z_s)

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
	DE_z = fnc_z(DE_slog);
	PI_z = fnc_z(PI_slog);
	PL_z = fnc_z(PL_slog);

	###! Assimilation efficiency lambdaa
	DE_lambda = ones(DE_N) .* 0.7;
	PI_lambda = ones(PI_N) .* 0.7;
	PL_lambda = ones(PL_N) .* 0.7;

	###! Kappa rule K as a function of body size
	DE_K = linspace(0,1,DE_N)
	PI_K = linspace(0,1,PI_N)
	PL_K = linspace(0,1,PL_N)

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
	DE_a = fnc_a(DE_s)
	PI_a = fnc_a(PI_s)
	PL_a = fnc_a(PL_s)

	###! Diet Preference Phi (j = prey, i = pred)
	#for each fish, get all prey and calculate preference
	#normalize so that most preferred prey preference == 1
	#distribute information into feeding kernels
	function fnc_pij(w1,w2) # Blanchard function to calculate diet preference
		# w1: body size of pred (g)
		# w2: body size of prey (g)
		beta = 2; # mean PPMR
		sigma = 1; # sd PPMR
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

	####! Piscivore (eats everything except zooplankton)
	PI_phi_PI = zeros(PI_N,PI_N);
	PI_phi_PL = zeros(PL_N,PI_N);
	PI_phi_DE = zeros(DE_N,PI_N);
	for i = 1:PI_N
		pij_pi = fnc_pij(PI_s[i],PI_s) # prey preference for PI
		pij_pl = fnc_pij(PI_s[i],PL_s) # prey preference for PL
		pij_de = fnc_pij(PI_s[i],DE_s) # prey preference for DE

		# normalize
		Sj = [pij_pi; pij_pl; pij_de];
		if sum(Sj)!=0.
			max = maximum(Sj);
			pij_pi = pij_pi / max;
			pij_pl = pij_pl / max;
			pij_de = pij_de / max;
		end

		#! put in diet preference kernels
		PI_phi_PI[:,i] = pij_pi
		PI_phi_PL[:,i] = pij_pl
		PI_phi_DE[:,i] = pij_de
	end
	#! detritus and zoo feeding preference
	PI_phi_Z = zeros(2,PI_N)
	PI_phi_W = zeros(PI_N)


	####! Planktivore (eats only zooplankton)
	PL_phi_Z = zeros(2,PL_N);
	for i = 1:PL_N
		pij_z  = fnc_pij(PI_s[i],Z_s) # prey preference for zoo

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
	PL_phi_W = zeros(PL_N)
	PL_phi_PI = zeros(PI_N,PL_N)
	PL_phi_PL = zeros(PL_N,PL_N)
	PL_phi_DE = zeros(DE_N,PL_N)

	####! Detritivore (eats only detritus and other detritivores)
	DE_phi_DE = zeros(DE_N,DE_N);
	for i = 1:DE_N
		pij_de = fnc_pij(DE_s[i],DE_s) # prey preference for DE

		# normalize
		Sj = [pij_de];
		if sum(Sj)!=0.
			max = maximum(Sj);
			pij_de = pij_de / max;
		end

		#! put in diet preference kernels
		DE_phi_DE[:,i]  = pij_de
	end
	#! other feeding kernels
	DE_phi_W  = ones(DE_N) # all eat detritus
	DE_phi_PI = zeros(PI_N,DE_N)
	DE_phi_PL = zeros(PL_N,DE_N)
	DE_phi_Z  = zeros(2,DE_N)


	#============ MAKE PARAMETER TYPE ==========#
	PRM_PI = PARAM(PI_N,PI_smin,PI_smax,PI_s,PI_slog,PI_z,PI_lambda,PI_K,PI_a,
				PI_phi_Z,PI_phi_PI,PI_phi_PL,PI_phi_DE,PI_phi_W);
	PRM_PL = PARAM(PL_N,PL_smin,PL_smax,PL_s,PL_slog,PL_z,PL_lambda,PL_K,PL_a,
				PL_phi_Z,PL_phi_PI,PL_phi_PL,PL_phi_DE,PL_phi_W);
	PRM_DE = PARAM(DE_N,DE_smin,DE_smax,DE_s,DE_slog,DE_z,DE_lambda,DE_K,DE_a,
				DE_phi_Z,DE_phi_PI,DE_phi_PL,DE_phi_DE,DE_phi_W);

	return PRM_PI, PRM_PL, PRM_DE
end




