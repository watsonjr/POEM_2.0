#### SELECT SIMULATION AND SET INPUTS
harv = 0 #0=no fishing; 1=fishing
make_parameters(harv) # make core parameters/constants

#! setup spinup (loop first year of COBALT)
COBALT = load("/Volumes/GFDL/POEM_JLD/Data_hindcast_PC_000120.jld"); # 1980
#! Add phenology params from csv file with ID as row
Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
global Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
global Sp = readdlm("./Data/Gaussian_spawn_2mo.csv",',');
YEARS = 100
global DAYS = 365

#! choose where to run the model
global GRD = load("./Data/Data_grid_hindcast_NOTflipped.jld")
XY = zeros(Int,360,200);
XY[GRD["ID"]] = collect(1:GRD["N"])
ids = [40319,42639,41782,36334,38309,42744,30051,41284,38003]
names = ["GB","EBS","OSP","HOT","BATS","NS","EEP","K2","S1"]

for L =1:9
	ID = ids[L]
	loc = names[L]

	const global NX = length(ID)
	phen=0;
	#! Initialize
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
	Med_d.td[1] = 0.0;
	Lrg_d.td[1] = 0.0;
	ENVR = sub_init_env(ID);

	#! Iterate forward in time
	for YR = 1:YEARS # years
		for DAY = 1:DT:DAYS # days
			###! Future time step
			sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
			DY+=1
		end #Days
	end #Years
end

#### THE PARAMETERS
#! Integration parameters
function make_parameters(harv)
	const global DT = 1.; # time step

	#! Amount of fishing
	if harv == 1
		const global LFISHING = 0.8/365.0
		const global MFISHING = LFISHING
		#const global MFISHING = 0.5 * LFISHING
	else
		const global MFISHING = 0
		const global LFISHING = 0
	end

	#! Benthic-pelagic coupling cutoff (depth, m)
	const global PI_be_cutoff = 200
	# 0:no coupling; 1:demersal coupled only; 2:pelagic & demersal coupled
	const global pdc = 2;

	#! body lengths (mm)
	const global L_s = 10^((log10(2)+log10(20))/2); # small
	const global L_m = 10^((log10(20)+log10(200))/2); # medium
	const global L_l = 10^((log10(200)+log10(2000))/2); # large

	##! Mass from length using Andersen & Beyer 2013
	# Convert from mm to cm and use their const coeff = 0.01g/cm3
	const global M_s = 0.01 * (0.1*L_s)^3;
	const global M_m = 0.01 * (0.1*L_m)^3;
	const global M_l = 0.01 * (0.1*L_l)^3;

	#! Ratio of initial and final body sizes per size-class
	const global Z_s = (0.01*(0.1*2)^3) / (0.01*(0.1*20)^3)
	const global Z_m = (0.01*(0.1*20)^3) / (0.01*(0.1*200)^3)
	const global Z_l = (0.01*(0.1*200)^3) / (0.01*(0.1*2000)^3)

	###! Assimilation efficiency lambda (constant across everything)
	const global Lambda = 0.7;

	###! Kappa rule K as a function of body size
	# K = fraction of energy consumed diverted to somatic growth
	const global K_l = 1
	const global K_j = 1
	const global K_a = 0

	###! Metabolism constants (activity and basal)
	const global fcrit = 0.40	# feeding level needed to meet resting metabolic demands; 0.05-0.2

	###! Consumption constants
	const global h = 60.0  		# h=40 g^(0.25)/yr at 10C in Cmax eq

	###! Transfer efficiency of detritus to benthic prey
	const global bent_eff = 0.30

	###! Reproductive efficiency
	const global rfrac = 1.0

	#! Fraction of time spent swimming (from Van Leeuwen)
	const global Tu_s = 1.0
	const global Tu_m = 1.0 #0.5
	const global Tu_l = 1.0 #0.1

	###! Background mortality
	#Currently increases from 0 to 0.01
	#Megrey et al =0.44/yr
	#Andersen & Beyer 2013 = 0.35 * 4.5 * s^(-0.25) (includes predation, excludes fishing)
	const global Nat_mrt = 0.44 / 365
	#0=none, 1=constant, 2=temp-dep, 3=large only, 4=large temp-dep
	const global MORT = 0

	###! Diet Preference Phi (j = prey, i = pred)
	#small forage fish eats medium zoo
	#small piscivores eats medium zoo
	#small detritivore eats medium zoo
	#medium forage fish eats medium & large zoo, all small fishes
	#medium piscivore eats medium & large zoo, all small fishes
	#medium detritivore eats detritus
	#large piscivore eats medium forage fish, medium piscivore, medium detritivore
	#large detritivore eats detritus, medium forage fish, medium piscivore, medium detrivore

	const global MF_phi_MZ = 0.1
	const global MF_phi_LZ = 1.0
	const global MF_phi_S = 1.0
	const global MP_phi_MZ = 0.1
	const global MP_phi_LZ = 1.0
	const global MP_phi_S = 1.0
	const global LP_phi_MF = 1.0
	const global LD_phi_MF = 1.0
	const global LP_phi_MP = 1.0
	const global LD_phi_MP = 1.0
	const global LP_phi_MD = 1.0
	const global LD_phi_MD = 1.0
end


#### THE INITIALIZATION
###! INIT FUNCTIONS
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
	X = 0.00001; # very small amount
	bio = ones(Float64,NX) * X

	# fraction of time spent in the pelagic
	tdp = ones(Float64,NX)

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
	mass = zeros(Float64,NX)
  BENT = detritus(mass)

	return Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT
end


#### THE FUNCTIONS
###! Get COBALT data
function get_COBALT!(COBALT,ID,DY,ENVR)
    ## Get data
    ENVR.Tp[:,1]  = COBALT["Tp"][ID,DY]
    ENVR.Tb[:,1]  = COBALT["Tb"][ID,DY]
    ENVR.Zm[:,1]  = COBALT["Zm"][ID,DY]
    ENVR.Zl[:,1]  = COBALT["Zl"][ID,DY]
    ENVR.det[:,1] = COBALT["det"][ID,DY]
    ENVR.dZm[:,1] = COBALT["dZm"][ID,DY]
    ENVR.dZl[:,1] = COBALT["dZl"][ID,DY]
    ENVR.U[:,1]   = COBALT["U"][ID,DY]
    ENVR.V[:,1]   = COBALT["V"][ID,DY]
    ENVR.T0p[:,1] = TrefP[ID]
    ENVR.T0b[:,1] = TrefB[ID]
    ENVR.Dthresh[:,1] = Dthresh[ID]
    ENVR.fZm[:,1] = zeros(Int64,NX)
    ENVR.fZl[:,1] = zeros(Int64,NX)
    ENVR.fB[:,1]  = zeros(Int64,NX)
    ENVR.H[:,1]   = GRD["Z"][ID]
    ENVR.A[:,1]   = GRD["AREA"][ID]
end

###! Fraction of time spent in pelagic (for piscivore)
function sub_tdif_pel(Z,bio1,bio2,biod)
  # bio1, bio2: pelagic prey
  # biod: demersal prey
	biop = bio1+bio2
	if Z < PI_be_cutoff
		tdif = biop ./ (biop+biod)
	else
		tdif = 1.0
	end
	return tdif
end

###! Fraction of time spent in pelagic (for demersal)
function sub_tdif_dem(Z,bio1,bio2,bio3,bio4)
  # bio1, bio2: pelagic prey
  # bio3, bio4: demersal prey
	biop = bio1+bio2
  biod = bio3+bio4
	if Z < PI_be_cutoff
		tdif = biop ./ (biop+biod)
	else
		tdif = 0.0
	end
	return tdif
end

###! Metabolism
function sub_met(Tp,Tb,tdif,wgt,L)
  #Tp: pelagic temp
  #Tb: bottom temp
  #tdif: frac pelagic time
  #wgt: ind weight of size class
	#fcrit: feeding level to meet resting respiration rate
  #cmax: max consumption rate
  #U: swimming speed
  temp = (Tp.*tdif) + (Tb.*(1.0-tdif))
  #Cmax
  #! Specific ingestion rate from Hartvig et al (g/g/day)
  cmax = (exp(0.063*(temp-15.0)) * h * wgt^(-0.25)) ./365.0
  #Metabolism
	bas = fcrit * cmax
  met = bas
  return met
end

###!  Encounter rates
function sub_enc(Tp,Tb,wgt,pred,prey,tpel,tprey,pref)
  # Tp: pelagic temp
  # Tb: bottom temp
  # wgt: ind weight of size class
  # pred: pred biomass density,
	# prey: prey biomass density,
	# A: predator search rate,
  # tpel: time spent in pelagic,
	# tprey: time spent in area with that prey item.
  # pref: preference for prey item
  temp = (Tp.*tpel) + (Tb.*(1.0-tpel))
  #! Specific clearance rates from Kiorboe & Hirst (m3/g/day)
  A = (exp(0.063*(temp-15.0)) * 10^(3.24) * wgt^(-0.24)) * (24e-3/9)
  #Encounter per predator, mult by biomass later
  enc = prey*A*tprey*pref
  return enc
end

###! Type I consumption
function sub_cons(Tp,Tb,tpel,wgt,enc)
  #Tp: pelagic temp
  #Tb: bottom temp
  #tpel: frac pelagic time
  #wgt: ind weight of size class
  #enc: array of all encountered food
	#! calculates consumption rate of first element of enc
  #Cmax
  temp = (Tp.*tpel) + (Tb.*(1.0-tpel))
  #! Specific ingestion rate from Hartvig et al (g/g/day)
  cmax = (exp(0.063*(temp-15.0)) * h * wgt^(-0.25)) ./365.0
  ENC = sum(enc) # total biomass encountered
	con = cmax .* enc[1] ./ (cmax + ENC) # Type II
  #con = cmax
  return con
end

###! Offline coupling
function sub_offline_zm(enc_1,enc_2,enc_3,enc_4,enc_5,bio_1,bio_2,bio_3,bio_4,bio_5,dZ)
  #Can't consume more than lost to higher predators in COBALT model
  con_1 = enc_1 * bio_1
  con_2 = enc_2 * bio_2
  con_3 = enc_3 * bio_3
  con_4 = enc_4 * bio_4
  con_5 = enc_5 * bio_5
  if (con_1 + con_2 + con_3 + con_4 + con_5) > dZ
  	frac1 = con_1 / (con_1 + con_2 + con_3 + con_4 + con_5)
    frac2 = con_2 / (con_1 + con_2 + con_3 + con_4 + con_5)
    frac3 = con_3 / (con_1 + con_2 + con_3 + con_4 + con_5)
    frac4 = con_4 / (con_1 + con_2 + con_3 + con_4 + con_5)
    frac5 = con_5 / (con_1 + con_2 + con_3 + con_4 + con_5)
  	out_1 = (frac1 * dZ) / bio_1
  	out_2 = (frac2 * dZ) / bio_2
    out_3 = (frac3 * dZ) / bio_3
    out_4 = (frac4 * dZ) / bio_4
    out_5 = (frac5 * dZ) / bio_5
  else
  	out_1 = enc_1
  	out_2 = enc_2
    out_3 = enc_3
    out_4 = enc_4
    out_5 = enc_5
  end
  zf = (out_1*bio_1 + out_2*bio_2 + out_3*bio_3 + out_4*bio_4 + out_5*bio_5) / dZ
  return out_1, out_2, out_3, out_4, out_5, zf
end

function sub_offline_zl(enc_1,enc_2,bio_1,bio_2,dZ)
  #Can't consume more than lost to higher predators in COBALT model
  con_1 = enc_1 * bio_1
  con_2 = enc_2 * bio_2
  if (con_1 + con_2) > dZ
  	frac1 = con_1 / (con_1 + con_2)
    frac2 = con_2 / (con_1 + con_2)
    out_1 = (frac1 * dZ) / bio_1
  	out_2 = (frac2 * dZ) / bio_2
  else
  	out_1 = enc_1
  	out_2 = enc_2
  end
  zf = (out_1*bio_1 + out_2*bio_2) / dZ
  return out_1, out_2, zf
end

function sub_offline_bent(enc_1,enc_2,bio_1,bio_2,B,det)
	#Can't consume more than available
  con_1 = enc_1 * bio_1
  con_2 = enc_2 * bio_2
  if (con_1 + con_2) > B
		frac1 = con_1 / (con_1 + con_2)
    frac2 = con_2 / (con_1 + con_2)
    out_1 = (frac1 * B) / bio_1
		out_2 = (frac2 * B) / bio_2
	else
		out_1 = enc_1
		out_2 = enc_2
	end
  bf = (out_1*bio_1 + out_2*bio_2) / B #/ det
	return out_1, out_2, bf
end

###! Consumption/Cmax
function sub_clev(con,Tp,Tb,tdif,wgt)
	#Cmax
  temp = (Tp.*tdif) + (Tb.*(1.0-tdif))
  #! Specific ingestion rate from Hartvig et al (g/g/day)
  cmax = (exp(0.063*(temp-15.0)) * h * wgt^(-0.25)) ./365.0
  #clev
	clev = con/cmax
  return clev
end

###! ENERGY AVAILABLE FOR GROWTH NU
function sub_nu(I,B,met)
	#Lambda: assim effic
	#I: specific ingestion rate
	#met: specific respiration rate
	#B: current biomass
	nu = (I*Lambda) - met
  prod = nu * B
  return nu, prod
end

###! ENERGY AVAILABLE FOR SOMATIC GROWTH
function sub_gamma(K,Z,nu,d,B,S)
  # convert predation mortality to biomass specific rate
	D = (d/B) + Nat_mrt
  kap=K;
	gg = ((kap*nu) - D)/(1-(Z^(1-(D/(kap*nu)))))
  if gg < 0 || isnan(gg)==true
		gamma = 0.0
	else
    gg = min(gg,nu)
		gamma = gg
	end
	return gamma
end

###! BIOMASS MADE FROM REPRODUCTION
function sub_rep(nu,K,S,egg)
  #nu: energy for growth or spawning
  #K: proportion allocated to growth
  #S: fraction of pop spawning at that time
  #egg: energy stored for later repro
  # NOTE: Still never going to accumulate biomass as muscle tissue
  # If it is spawning season, it gets spawned
  # If it is not spawning season, it gets stored as repro energy
  # Need to determine a set fraction of energy that gets converted to larvae?
  if K<1.0
      if nu > 0.0
        rho = (1.0-K) * nu  #energy available for from eating
      else
        rho = 0.0
      end
      if S>0.0
        rep = rho + S * egg         #fraction of pop reproducing now
        egg = (1.0-S) * egg         #rest gets stored for later
      else
        rep = 0.0
        egg = egg + rho
      end
  else
    rep = 0.0
    egg = 0.0
  end
	return rep, egg
end

###! Biomass recruiting to size-class (g m-2 d-1)
function sub_rec(X,bio,wgt)
	# X could be biomass of eggs (for larval class) or maturing from smaller sizes
  if (wgt==M_s)
    rec = rfrac * X * bio
  else
    rec = X * bio
  end
	return rec
end

###! Temp-dep natural mortality
function sub_nmort(Tp,Tb,tpel,wgt)
  #Tp: pelagic temp
  #Tb: bottom temp
  #tpel: frac pelagic time
  if (MORT==0) # None
    nmort = 0.0
  end
  if (MORT==1) # Constant
    nmort = Nat_mrt
  end
  if (MORT==2) # Temperature-dependent mortality
    temp = (Tp.*tpel) + (Tb.*(1.0-tpel))
    nmort = exp(0.063*(temp-15.0)) * Nat_mrt
  end
  if (MORT==3) # Large fishes only
    if (wgt == M_l)
      nmort = Nat_mrt
    else
      nmort = 0.0
    end
  end
  if (MORT==4) # Large fishes only w/ temp-dep
    if (wgt == M_l)
      temp = (Tp.*tpel) + (Tb.*(1.0-tpel))
      nmort = exp(0.063*(temp-15.0)) * Nat_mrt
    else
      nmort = 0.0
    end
  end
  return nmort
end

###! Update biomass
function sub_update_fi(bio_in,rec,nu,rep,gamma,die,egg,nmort)
	# all inputs except rec are in g g-1 d-1; rec is g d-1
	# rec = rep from smaller size class = TOTAL biomass gained from recruitment
	# grw = nu = somatic energy for growth within size class
  # store = egg = energy stored for later egg production
	# rep = rep =  energy lost to egg production
	# mat = gamma = energy lost to maturation to larger size class
	# nmort = natural mortality
	# die = predator mort = biomass lost to predation
  db = rec + ((nu - egg - rep - gamma - nmort) * bio_in) + (egg * bio_in) - die
  bio_out =  bio_in + db
end

function sub_update_be(bio_in,con,bio)
	#bio_in: benthic biomass
	#con: specific consumption by M&L demersals
	#bio: biomass of M&L demersals
	#die: benthic biomass consumed by M&L demersals
	#bio_out: remaining benthic biomass
  die = con.*bio
  bio_out = bio_in - sum(die)
end

####! Fishing
function sub_fishing_rate(bio,wgt)
	if (wgt==M_m)
    caught = bio * MFISHING
    bio -= caught
  elseif (wgt==M_l)
    caught = bio * LFISHING
    bio -= caught
  else
    bio = bio
    caught = 0.0
	end
	return bio, caught
end

###! Forward Euler checks
#Prevent biomass<0
function sub_check!(bio)
	ID = find(bio .< 0)
	bio[ID] = eps()
end


#### THE MODEL
###! DEMOGRAPHIC CALCULATIONS
function sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT)

	###! COBALT information
	get_COBALT!(COBALT,ID,DY,ENVR)

	for JD = 1:NX

		#! update benthic biomass with new detritus avail at that time step
		BENT.mass[JD] = BENT.mass[JD] + bent_eff*ENVR.det[JD]

		#! pelagic-demersal coupling
		#Lrg_p: fraction of time large piscivores spends in pelagic
		#Lrg_d: fraction of time large demersals spends in pelagic
		if (pdc == 0)
			Lrg_p.td[JD] = 1.0
			Lrg_d.td[JD] = 0.0
		elseif (pdc == 1)
			Lrg_p.td[JD] = 1.0
			Lrg_d.td[JD] = sub_tdif_dem(ENVR.H[JD],Med_f.bio[JD],Med_p.bio[JD],Med_d.bio[JD],BENT.mass[JD])
		elseif (pdc == 2)
	  	Lrg_p.td[JD] = sub_tdif_pel(ENVR.H[JD],Med_f.bio[JD],Med_p.bio[JD],Med_d.bio[JD])
			Lrg_d.td[JD] = sub_tdif_dem(ENVR.H[JD],Med_f.bio[JD],Med_p.bio[JD],Med_d.bio[JD],BENT.mass[JD])
		end

		#! metabolism
		Sml_f.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Sml_f.td[JD],M_s,L_s)
		Sml_p.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Sml_p.td[JD],M_s,L_s)
		Sml_d.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Sml_d.td[JD],M_s,L_s)
		Med_f.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m,L_m)
		Med_p.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,L_m)
		Med_d.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Med_d.td[JD],M_m,L_m)
		Lrg_p.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l,L_l)
		Lrg_d.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,L_l)

		#! encounter rates
		#sub_enc(Tp,Tb,wgt,pred,prey,td,tprey,pref)
		Sml_f.enc_zm[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_s,Sml_f.bio[JD],ENVR.Zm[JD],Sml_f.td[JD],Sml_f.td[JD],1)
		Sml_p.enc_zm[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_s,Sml_p.bio[JD],ENVR.Zm[JD],Sml_p.td[JD],Sml_f.td[JD],1)
		Sml_d.enc_zm[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_s,Sml_d.bio[JD],ENVR.Zm[JD],Sml_d.td[JD],Sml_f.td[JD],1)

		Med_f.enc_zm[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_f.bio[JD],ENVR.Zm[JD],Med_f.td[JD],Med_f.td[JD],MF_phi_MZ)
		Med_f.enc_zl[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_f.bio[JD],ENVR.Zl[JD],Med_f.td[JD],Med_f.td[JD],MF_phi_LZ)
		Med_f.enc_f[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_f.bio[JD],Sml_f.bio[JD],Med_f.td[JD],Med_f.td[JD],MF_phi_S)
		Med_f.enc_p[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_f.bio[JD],Sml_p.bio[JD],Med_f.td[JD],Med_f.td[JD],MF_phi_S)
		Med_f.enc_d[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_f.bio[JD],Sml_d.bio[JD],Med_f.td[JD],Med_f.td[JD],MF_phi_S)

		Med_p.enc_zm[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_p.bio[JD],ENVR.Zm[JD],Med_p.td[JD],Med_p.td[JD],MP_phi_MZ)
		Med_p.enc_zl[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_p.bio[JD],ENVR.Zl[JD],Med_p.td[JD],Med_p.td[JD],MP_phi_LZ)
		Med_p.enc_f[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_p.bio[JD],Sml_f.bio[JD],Med_p.td[JD],Med_p.td[JD],MP_phi_S)
		Med_p.enc_p[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_p.bio[JD],Sml_p.bio[JD],Med_p.td[JD],Med_p.td[JD],MP_phi_S)
		Med_p.enc_d[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_p.bio[JD],Sml_d.bio[JD],Med_p.td[JD],Med_p.td[JD],MP_phi_S)

		Med_d.enc_be[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,Med_d.bio[JD],BENT.mass[JD],Med_d.td[JD],1-Med_d.td[JD],1)

		Lrg_p.enc_f[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,Lrg_p.bio[JD],Med_f.bio[JD],Lrg_p.td[JD],Lrg_p.td[JD],LP_phi_MF)
		Lrg_p.enc_p[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,Lrg_p.bio[JD],Med_p.bio[JD],Lrg_p.td[JD],Lrg_p.td[JD],LP_phi_MP)
		Lrg_p.enc_d[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,Lrg_p.bio[JD],Med_d.bio[JD],Lrg_p.td[JD],1-Lrg_p.td[JD],LP_phi_MD)

		Lrg_d.enc_f[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,Lrg_d.bio[JD],Med_f.bio[JD],Lrg_d.td[JD],Lrg_d.td[JD],LD_phi_MF)
		Lrg_d.enc_p[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,Lrg_d.bio[JD],Med_p.bio[JD],Lrg_d.td[JD],Lrg_d.td[JD],LD_phi_MP)
		Lrg_d.enc_d[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,Lrg_d.bio[JD],Med_d.bio[JD],Lrg_d.td[JD],1-Lrg_d.td[JD],LD_phi_MD)
		Lrg_d.enc_be[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,Lrg_d.bio[JD],BENT.mass[JD],Lrg_d.td[JD],1-Lrg_d.td[JD],1)

		#! Consumption rates
		Sml_f.con_zm[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Sml_f.td[JD],M_s,Sml_f.enc_zm[JD])
		Sml_p.con_zm[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Sml_p.td[JD],M_s,Sml_p.enc_zm[JD])
		Sml_d.con_zm[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Sml_d.td[JD],M_s,Sml_d.enc_zm[JD])

		Med_f.con_zm[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m,[Med_f.enc_zm[JD],Med_f.enc_zl[JD],Med_f.enc_f[JD],Med_f.enc_p[JD],Med_f.enc_d[JD]])
		Med_f.con_zl[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m,[Med_f.enc_zl[JD],Med_f.enc_zm[JD],Med_f.enc_f[JD],Med_f.enc_p[JD],Med_f.enc_d[JD]])
		Med_f.con_f[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m,[Med_f.enc_f[JD],Med_f.enc_zm[JD],Med_f.enc_zl[JD],Med_f.enc_p[JD],Med_f.enc_d[JD]])
		Med_f.con_p[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m,[Med_f.enc_p[JD],Med_f.enc_zm[JD],Med_f.enc_zl[JD],Med_f.enc_f[JD],Med_f.enc_d[JD]])
		Med_f.con_d[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m,[Med_f.enc_d[JD],Med_f.enc_zm[JD],Med_f.enc_zl[JD],Med_f.enc_f[JD],Med_f.enc_p[JD]])

		Med_p.con_zm[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,[Med_p.enc_zm[JD],Med_p.enc_zl[JD],Med_p.enc_f[JD],Med_p.enc_p[JD],Med_p.enc_d[JD]])
		Med_p.con_zl[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,[Med_p.enc_zl[JD],Med_p.enc_zm[JD],Med_p.enc_f[JD],Med_p.enc_p[JD],Med_p.enc_d[JD]])
		Med_p.con_f[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,[Med_p.enc_f[JD],Med_p.enc_zm[JD],Med_p.enc_zl[JD],Med_p.enc_p[JD],Med_p.enc_d[JD]])
		Med_p.con_p[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,[Med_p.enc_p[JD],Med_p.enc_zm[JD],Med_p.enc_zl[JD],Med_p.enc_f[JD],Med_p.enc_d[JD]])
		Med_p.con_d[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,[Med_p.enc_d[JD],Med_p.enc_zm[JD],Med_p.enc_zl[JD],Med_p.enc_f[JD],Med_p.enc_p[JD]])

		Med_d.con_be[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_d.td[JD],M_m,Med_d.enc_be[JD])

		Lrg_p.con_f[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l,[Lrg_p.enc_f[JD],Lrg_p.enc_p[JD],Lrg_p.enc_d[JD]])
		Lrg_p.con_p[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l,[Lrg_p.enc_p[JD],Lrg_p.enc_f[JD],Lrg_p.enc_d[JD]])
		Lrg_p.con_d[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l,[Lrg_p.enc_d[JD],Lrg_p.enc_p[JD],Lrg_p.enc_f[JD]])

		Lrg_d.con_f[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,[Lrg_d.enc_f[JD],Lrg_d.enc_p[JD],Lrg_d.enc_d[JD],Lrg_d.enc_be[JD]])
		Lrg_d.con_p[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,[Lrg_d.enc_p[JD],Lrg_d.enc_f[JD],Lrg_d.enc_d[JD],Lrg_d.enc_be[JD]])
		Lrg_d.con_d[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,[Lrg_d.enc_d[JD],Lrg_d.enc_p[JD],Lrg_d.enc_f[JD],Lrg_d.enc_be[JD]])
		Lrg_d.con_be[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,[Lrg_d.enc_be[JD],Lrg_d.enc_f[JD],Lrg_d.enc_p[JD],Lrg_d.enc_d[JD]])


		#! Offline coupling
		#Zooplankton consumption cannot exceed amount lost to higher predation in COBALT runs
		Sml_f.con_zm[JD],Sml_p.con_zm[JD],Sml_d.con_zm[JD],Med_f.con_zm[JD],Med_p.con_zm[JD],ENVR.fZm[JD] = sub_offline_zm(Sml_f.con_zm[JD],Sml_p.con_zm[JD],Sml_d.con_zm[JD],Med_f.con_zm[JD],Med_p.con_zm[JD],Sml_f.bio[JD],Sml_p.bio[JD],Sml_d.bio[JD],Med_f.bio[JD],Med_p.bio[JD],ENVR.dZm[JD])
		Med_f.con_zl[JD],Med_p.con_zl[JD],ENVR.fZl[JD] = sub_offline_zl(Med_f.con_zl[JD],Med_p.con_zl[JD],Med_f.bio[JD],Med_p.bio[JD],ENVR.dZl[JD])
		#Benthic material consumption cannot exceed amount present
		Med_d.con_be[JD], Lrg_d.con_be[JD], ENVR.fB[JD] = sub_offline_bent(Med_d.con_be[JD],Lrg_d.con_be[JD],Med_d.bio[JD],Lrg_d.bio[JD],BENT.mass[JD],ENVR.det[JD])

		#! total consumption rates (could factor in handling times here; g m-2 d-1)
		Sml_f.I[JD] = Sml_f.con_zm[JD]
		Sml_p.I[JD] = Sml_p.con_zm[JD]
		Sml_d.I[JD] = Sml_d.con_zm[JD]
		Med_f.I[JD] = Med_f.con_zm[JD] + Med_f.con_zl[JD] + Med_f.con_f[JD] + Med_f.con_p[JD] + Med_f.con_d[JD]
		Med_p.I[JD] = Med_p.con_zm[JD] + Med_p.con_zl[JD] + Med_p.con_f[JD] + Med_p.con_p[JD] + Med_p.con_d[JD]
		Med_d.I[JD] = Med_d.con_be[JD]
		Lrg_p.I[JD] = Lrg_p.con_f[JD] + Lrg_p.con_p[JD] + Lrg_p.con_d[JD]
		Lrg_d.I[JD] = Lrg_d.con_f[JD] + Lrg_d.con_p[JD] + Lrg_d.con_d[JD] + Lrg_d.con_be[JD]

		#! consumption related to Cmax
		Sml_f.clev[JD] = sub_clev(Sml_f.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Sml_f.td[JD],M_s)
		Sml_p.clev[JD] = sub_clev(Sml_p.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Sml_p.td[JD],M_s)
		Sml_d.clev[JD] = sub_clev(Sml_d.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Sml_d.td[JD],M_s)
		Med_f.clev[JD] = sub_clev(Med_f.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m)
		Med_p.clev[JD] = sub_clev(Med_p.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m)
		Med_d.clev[JD] = sub_clev(Med_d.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Med_d.td[JD],M_m)
		Lrg_p.clev[JD] = sub_clev(Lrg_p.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l)
		Lrg_d.clev[JD] = sub_clev(Lrg_d.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l)

		#! death rates (g m-2 d-1)
		Sml_f.die[JD] = Med_p.con_f[JD]*Med_p.bio[JD] + Med_f.con_f[JD]*Med_f.bio[JD]
		Sml_p.die[JD] = Med_p.con_p[JD]*Med_p.bio[JD] + Med_f.con_p[JD]*Med_f.bio[JD]
		Sml_d.die[JD] = Med_p.con_d[JD]*Med_p.bio[JD] + Med_f.con_d[JD]*Med_f.bio[JD]
		Med_f.die[JD] = Lrg_p.con_f[JD]*Lrg_p.bio[JD] + Lrg_d.con_f[JD]*Lrg_d.bio[JD]
		Med_p.die[JD] = Lrg_p.con_p[JD]*Lrg_p.bio[JD] + Lrg_d.con_p[JD]*Lrg_d.bio[JD]
		Med_d.die[JD] = Lrg_p.con_d[JD]*Lrg_p.bio[JD] + Lrg_d.con_d[JD]*Lrg_d.bio[JD]

		#! predation rates (m-2 d-1)
		Sml_f.pred[JD] = Sml_f.die[JD] / Sml_f.bio[JD]
		Sml_p.pred[JD] = Sml_p.die[JD] / Sml_p.bio[JD]
		Sml_d.pred[JD] = Sml_d.die[JD] / Sml_d.bio[JD]
		Med_f.pred[JD] = Med_f.die[JD] / Med_f.bio[JD]
		Med_p.pred[JD] = Med_p.die[JD] / Med_p.bio[JD]
		Med_d.pred[JD] = Med_d.die[JD] / Med_d.bio[JD]

		#! natural mortality rates
		Sml_f.nmort[JD] = sub_nmort(ENVR.Tp[JD],ENVR.Tb[JD],Sml_f.td[JD],M_s)
		Sml_p.nmort[JD] = sub_nmort(ENVR.Tp[JD],ENVR.Tb[JD],Sml_p.td[JD],M_s)
		Sml_d.nmort[JD] = sub_nmort(ENVR.Tp[JD],ENVR.Tb[JD],Sml_d.td[JD],M_s)
		Med_f.nmort[JD] = sub_nmort(ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m)
		Med_p.nmort[JD] = sub_nmort(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m)
		Med_d.nmort[JD] = sub_nmort(ENVR.Tp[JD],ENVR.Tb[JD],Med_d.td[JD],M_m)
		Lrg_p.nmort[JD] = sub_nmort(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l)
		Lrg_d.nmort[JD] = sub_nmort(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l)

		#! Degree days
		Med_f.DD[JD] = sub_degday(Med_f.DD[JD],ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],ENVR.T0p[JD],Med_f.S[JD,:],DY)
		Lrg_p.DD[JD] = sub_degday(Lrg_p.DD[JD],ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],ENVR.T0p[JD],Lrg_p.S[JD,:],DY)
		#Lrg_d.DD[JD] = sub_degday(Lrg_d.DD[JD],ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],ENVR.T0b[JD],Lrg_d.S[JD,:],DY)
		#Assume demersal spawn at same time as pelagic b/c larvae also need spring bloom
		Lrg_d.DD[JD] = sub_degday(Lrg_d.DD[JD],ENVR.Tp[JD],ENVR.Tb[JD],1-Lrg_d.td[JD],ENVR.T0b[JD],Lrg_d.S[JD,:],DY)

		#! Spawning flag determined from DD, dthresh
		Med_f.S[JD,:], Med_f.DD[JD] = sub_kflag(Med_f.S[JD,:],Med_f.DD[JD],ENVR.Dthresh[JD],DY);
		Lrg_d.S[JD,:], Lrg_d.DD[JD] = sub_kflag(Lrg_d.S[JD,:],Lrg_d.DD[JD],ENVR.Dthresh[JD],DY);
		Lrg_p.S[JD,:], Lrg_p.DD[JD] = sub_kflag(Lrg_p.S[JD,:],Lrg_p.DD[JD],ENVR.Dthresh[JD],DY);

		#! energy available for somatic growth nu
		Sml_f.nu[JD], Sml_f.prod[JD] = sub_nu(Sml_f.I[JD],Sml_f.bio[JD],Sml_f.met[JD])
		Sml_p.nu[JD], Sml_p.prod[JD] = sub_nu(Sml_p.I[JD],Sml_p.bio[JD],Sml_p.met[JD])
		Sml_d.nu[JD], Sml_d.prod[JD] = sub_nu(Sml_d.I[JD],Sml_d.bio[JD],Sml_d.met[JD])
		Med_f.nu[JD], Med_f.prod[JD] = sub_nu(Med_f.I[JD],Med_f.bio[JD],Med_f.met[JD])
		Med_p.nu[JD], Med_p.prod[JD] = sub_nu(Med_p.I[JD],Med_p.bio[JD],Med_p.met[JD])
		Med_d.nu[JD], Med_d.prod[JD] = sub_nu(Med_d.I[JD],Med_d.bio[JD],Med_d.met[JD])
		Lrg_p.nu[JD], Lrg_p.prod[JD] = sub_nu(Lrg_p.I[JD],Lrg_p.bio[JD],Lrg_p.met[JD])
		Lrg_d.nu[JD], Lrg_d.prod[JD] = sub_nu(Lrg_d.I[JD],Lrg_d.bio[JD],Lrg_d.met[JD])

		#! maturation (note subscript on Kappa is larvae, juv, adult)
		Sml_f.gamma[JD] = sub_gamma(K_l,Z_s,Sml_f.nu[JD],Sml_f.die[JD],Sml_f.bio[JD],Sml_f.S[JD,DY])
		Sml_p.gamma[JD] = sub_gamma(K_l,Z_s,Sml_p.nu[JD],Sml_p.die[JD],Sml_p.bio[JD],Sml_p.S[JD,DY])
		Sml_d.gamma[JD] = sub_gamma(K_l,Z_s,Sml_d.nu[JD],Sml_d.die[JD],Sml_d.bio[JD],Sml_d.S[JD,DY])
		Med_f.gamma[JD] = sub_gamma(K_a,Z_m,Med_f.nu[JD],Med_f.die[JD],Med_f.bio[JD],Med_f.S[JD,DY])
		Med_p.gamma[JD] = sub_gamma(K_j,Z_m,Med_p.nu[JD],Med_p.die[JD],Med_p.bio[JD],Med_p.S[JD,DY])
		Med_d.gamma[JD] = sub_gamma(K_j,Z_m,Med_d.nu[JD],Med_d.die[JD],Med_d.bio[JD],Med_d.S[JD,DY])
		Lrg_p.gamma[JD] = sub_gamma(K_a,Z_l,Lrg_p.nu[JD],Lrg_p.die[JD],Lrg_p.bio[JD],Lrg_p.S[JD,DY])
		Lrg_d.gamma[JD] = sub_gamma(K_a,Z_l,Lrg_d.nu[JD],Lrg_d.die[JD],Lrg_d.bio[JD],Lrg_d.S[JD,DY])

		#! egg production (by med and large size classes only)
		Sml_f.rep[JD],Sml_f.egg[JD] = sub_rep(Sml_f.nu[JD],K_l,Sml_f.S[JD,DY],Sml_f.egg[JD])
		Sml_p.rep[JD],Sml_p.egg[JD] = sub_rep(Sml_p.nu[JD],K_l,Sml_p.S[JD,DY],Sml_p.egg[JD])
		Sml_d.rep[JD],Sml_d.egg[JD] = sub_rep(Sml_d.nu[JD],K_l,Sml_d.S[JD,DY],Sml_d.egg[JD])
		Med_f.rep[JD],Med_f.egg[JD] = sub_rep(Med_f.nu[JD],K_a,Med_f.S[JD,DY],Med_f.egg[JD])
		Med_p.rep[JD],Med_p.egg[JD] = sub_rep(Med_p.nu[JD],K_j,Med_p.S[JD,DY],Med_p.egg[JD])
		Med_d.rep[JD],Med_d.egg[JD] = sub_rep(Med_d.nu[JD],K_j,Med_d.S[JD,DY],Med_d.egg[JD])
		Lrg_p.rep[JD],Lrg_p.egg[JD] = sub_rep(Lrg_p.nu[JD],K_a,Lrg_p.S[JD,DY],Lrg_p.egg[JD])
		Lrg_d.rep[JD],Lrg_d.egg[JD] = sub_rep(Lrg_d.nu[JD],K_a,Lrg_d.S[JD,DY],Lrg_d.egg[JD])

		#! recruitment (from smaller size class)
		Sml_f.rec[JD] = sub_rec(Med_f.rep[JD],Med_f.bio[JD],M_s)
		Sml_p.rec[JD] = sub_rec(Lrg_p.rep[JD],Lrg_p.bio[JD],M_s)
		Sml_d.rec[JD] = sub_rec(Lrg_d.rep[JD],Lrg_d.bio[JD],M_s)
		Med_f.rec[JD] = sub_rec(Sml_f.gamma[JD],Sml_f.bio[JD],M_m)
		Med_p.rec[JD] = sub_rec(Sml_p.gamma[JD],Sml_p.bio[JD],M_m)
		Med_d.rec[JD] = sub_rec(Sml_d.gamma[JD],Sml_d.bio[JD],M_m)
		Lrg_p.rec[JD] = sub_rec(Med_p.gamma[JD],Med_p.bio[JD],M_l)
		Lrg_d.rec[JD] = sub_rec(Med_d.gamma[JD],Med_d.bio[JD],M_l)

		#! Mass balance
		BENT.mass[JD] = sub_update_be(BENT.mass[JD],[Med_d.con_be[JD],Lrg_d.con_be[JD]],[Med_d.bio[JD],Lrg_d.bio[JD]])

		Sml_f.bio[JD] = sub_update_fi(Sml_f.bio[JD],Sml_f.rec[JD],Sml_f.nu[JD],
								   Sml_f.rep[JD],Sml_f.gamma[JD],Sml_f.die[JD],Sml_f.egg[JD],Sml_f.nmort[JD])
		Sml_p.bio[JD] = sub_update_fi(Sml_p.bio[JD],Sml_p.rec[JD],Sml_p.nu[JD],
								   Sml_p.rep[JD],Sml_p.gamma[JD],Sml_p.die[JD],Sml_p.egg[JD],Sml_p.nmort[JD])
		Sml_d.bio[JD] = sub_update_fi(Sml_d.bio[JD],Sml_d.rec[JD],Sml_d.nu[JD],
								   Sml_d.rep[JD],Sml_d.gamma[JD],Sml_d.die[JD],Sml_d.egg[JD],Sml_d.nmort[JD])

		Med_f.bio[JD] = sub_update_fi(Med_f.bio[JD],Med_f.rec[JD],Med_f.nu[JD],
								   Med_f.rep[JD],Med_f.gamma[JD],Med_f.die[JD],Med_f.egg[JD],Med_f.nmort[JD])
		Med_p.bio[JD] = sub_update_fi(Med_p.bio[JD],Med_p.rec[JD],Med_p.nu[JD],
								   Med_p.rep[JD],Med_p.gamma[JD],Med_p.die[JD],Med_p.egg[JD],Med_p.nmort[JD])
		Med_d.bio[JD] = sub_update_fi(Med_d.bio[JD],Med_d.rec[JD],Med_d.nu[JD],
								   Med_d.rep[JD],Med_d.gamma[JD],Med_d.die[JD],Med_d.egg[JD],Med_d.nmort[JD])

		Lrg_p.bio[JD] = sub_update_fi(Lrg_p.bio[JD],Lrg_p.rec[JD],Lrg_p.nu[JD],
								   Lrg_p.rep[JD],Lrg_p.gamma[JD],Lrg_p.die[JD],Lrg_p.egg[JD],Lrg_p.nmort[JD])
	 	Lrg_d.bio[JD] = sub_update_fi(Lrg_d.bio[JD],Lrg_d.rec[JD],Lrg_d.nu[JD],
								   Lrg_d.rep[JD],Lrg_d.gamma[JD],Lrg_d.die[JD],Lrg_d.egg[JD],Lrg_d.nmort[JD])

		#! Fishing by rate
		Med_f.bio[JD], Med_f.caught[JD] = sub_fishing_rate(Med_f.bio[JD],M_m)
		Med_p.bio[JD], Med_p.caught[JD] = sub_fishing_rate(Med_p.bio[JD],M_m)
		Med_d.bio[JD], Med_d.caught[JD] = sub_fishing_rate(Med_d.bio[JD],M_m)
		Lrg_p.bio[JD], Lrg_p.caught[JD] = sub_fishing_rate(Lrg_p.bio[JD],M_l)
		Lrg_d.bio[JD], Lrg_d.caught[JD] = sub_fishing_rate(Lrg_d.bio[JD],M_l)

	end

	#! Forward Euler checks for demographics and movement
	sub_check!(Sml_f.bio);
	sub_check!(Sml_p.bio);
	sub_check!(Sml_d.bio);
	sub_check!(Med_f.bio);
	sub_check!(Med_p.bio);
	sub_check!(Med_d.bio);
	sub_check!(Lrg_p.bio);
	sub_check!(Lrg_d.bio);

end
