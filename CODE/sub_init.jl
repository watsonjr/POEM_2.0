
#============== INITIAL CONDITIONS =============#
function sub_init_env(ID)
	###! Number of spatial cells
	NX = length(ID)
	###! environment
	ENV_Tp  = Array(Float64,NX)
	ENV_Tb  = Array(Float64,NX)
	ENV_Zm  = Array(Float64,NX)
	ENV_Zl  = Array(Float64,NX)
	ENV_dZm = Array(Float64,NX)
	ENV_dZl = Array(Float64,NX)
	ENV_det = Array(Float64,NX)
	ENV_U   = Array(Float64,NX)
	ENV_V   = Array(Float64,NX)
	ENV_T0p = Array(Float64,NX)
	ENV_T0b = Array(Float64,NX)
	ENV_Dthresh = Array(Float64,NX)
	ENV_fZm = Array(Float64,NX)
	ENV_fZl = Array(Float64,NX)
	ENV_fB = Array(Float64,NX)
	ENV_H  = Array(Float64,NX)
	ENV_A  = Array(Float64,NX)
	ENVR = environment(ENV_Tp,ENV_Tb,ENV_Zm,ENV_Zl,ENV_dZm,ENV_dZl,ENV_det,ENV_U,ENV_V,ENV_T0p,ENV_T0b,ENV_Dthresh,ENV_fZm,ENV_fZl,ENV_fB,ENV_H,ENV_A)
end

function sub_init_fish(ID,phen)

	#===== VARIABLES =====#
	###! Number of spatial cells
	NX = length(ID)

	###! fish
	#! biomass
	#X = 0.17981663628808964; # mean Zm off Spain
	X = 1.0e-5; # very small amount
	#X = 1.0;
	bio = ones(Float64,NX) * X

	# fraction of time spent in the pelagic
	tdp = ones(Float64,NX)
	#tdd = zeros(Float64,NX)

	# metabolism
	met = zeros(Float64,NX)

	#! encounter rates between fish and zoo
	enc_f = zeros(Float64,NX)
	enc_p = zeros(Float64,NX)
	enc_d = zeros(Float64,NX)
	enc_zm = zeros(Float64,NX)
	enc_zl = zeros(Float64,NX)
	enc_be = zeros(Float64,NX)

	#! consumption rates between fish and zoo
	con_f = zeros(Float64,NX)
	con_p = zeros(Float64,NX)
	con_d = zeros(Float64,NX)
	con_zm = zeros(Float64,NX)
	con_zl = zeros(Float64,NX)
	con_be = zeros(Float64,NX)

  # mass ingested (I)
  I = zeros(Float64,NX) # fish

  # mass lost to predation (g d-1)
  die = zeros(Float64,NX)

	# predation rate (d-1)
  pred = zeros(Float64,NX)

	# natural mortality rate (d-1)
  nmort = zeros(Float64,NX)

	# production
  prod = zeros(Float64,NX)

  # total energy available for growth
  nu = zeros(Float64,NX)

  # energy available for somatic growth
  gamma = zeros(Float64,NX)

	# energy available for later repro/stored biomass for repro
  egg = zeros(Float64,NX)

 	#! total biomass to reproduction
 	rep = zeros(Float64,NX)

 	#! total biomass to reproduction
 	rec = zeros(Float64,NX)

	#! degree days accumulated that year
	DD = zeros(Float64,NX)

	#! spawning flag
	# frac spawning at any given time
	if (phen == 1)
		S = zeros(Float64,NX,DAYS)
	else
		S = ones(Float64,NX,DAYS)
	end

	#! Con/Cmax
	clev = zeros(Float64,NX)

	#! Fishing harvest
	caught = zeros(Float64,NX)

	# assign to small forage fish, piscivore and detrivore
	Sml_f = fish(bio,tdp,met,enc_f,enc_p,enc_d,enc_zm,enc_zl,enc_be,con_f,con_p,con_d,con_zm,con_zl,con_be,I,nu,gamma,die,rep,rec,DD,S,egg,clev,prod,pred,nmort,caught)
	Sml_p = deepcopy(Sml_f)
	Sml_d = deepcopy(Sml_f)
	Med_f = deepcopy(Sml_f)
	Med_p = deepcopy(Sml_f)
	Lrg_p = deepcopy(Sml_f)
	Med_d = deepcopy(Sml_f)
	Lrg_d = deepcopy(Sml_f)

	###! Detritus
	mass = ones(Float64,NX) * X
  BENT = detritus(mass)

	return Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT
end
