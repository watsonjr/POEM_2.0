#============== Parameters of the model =============#
#============ VARIABLE TYPE==========#
#! i.e. those things that change with time
type fish
	bio::Array{Float64} # biomass g
	td::Array{Float64} # fraction of time spent in the pelagic
	met::Array{Float64} # metabolism g g-1 d-1
	enc_f::Array{Float64}  # encounter rate
	enc_p::Array{Float64}  # encounter rate
	enc_d::Array{Float64}  # encounter rate
	enc_zm::Array{Float64} # encounter rate
	enc_zl::Array{Float64} # encounter rate
	enc_be::Array{Float64} # encounter rate
	con_f::Array{Float64}  # consumption rate
	con_p::Array{Float64}  # consumption rate
	con_d::Array{Float64}  # consumption rate
	con_zm::Array{Float64} # consumption rate
	con_zl::Array{Float64} # consumption rate
	con_be::Array{Float64} # consumption rate
	I::Array{Float64} # total fish ingestion g d-1
	nu::Array{Float64}  # total energy for growth g g-1 d-1
	gamma::Array{Float64} # energy for somatic growth g g-1 d-1
	die::Array{Float64} # total death g d-1
	rep::Array{Float64} # total biomass reproduced g d-1
	rec::Array{Float64} # total biomass reproduced g d-1
	DD::Array{Float64} # degree days accumulated that year, Celsius
	S::Array{Float64} # spawning flag
	egg::Array{Float64} # stored energy/biomass for reproduction
	clev::Array{Float64} # Con/Cmax
end

type detritus
	mass::Array{Float64} # biomass in detrital pool g m-2
end

type environment
	Tp::Array{Float64} # pelagic temperature
	Tb::Array{Float64} # bottom temperature
	Zm::Array{Float64} # medium zooplankton
	Zl::Array{Float64} # large zoo
	dZm::Array{Float64} # mortality rate of med zoo (COBALT)
	dZl::Array{Float64} # mortality rate of large zoo (COBALT)
	det::Array{Float64} # detrital flux
	U::Array{Float64} # U current speed
	V::Array{Float64} # U current speed
	T0p::Array{Float64}			# spawning reference temp pelagic
	T0b::Array{Float64}			# spawning reference temp benthic
	Dthresh::Array{Float64}	# spawning degree day threshold
	fZm::Array{Float64} # frac med zoop losses consumed
	fZl::Array{Float64} # frac lrg zoop losses consumed
	fB::Array{Float64} # frac detr flux consumed
end

#============= PARAMETER TYPE ==========#
function make_parameters(harv)
	#! Input it switch for fishing (1 is yes, 0 is no)

	#! Grid parameters
	#global GRD = load("./Data/JLD/Data_grid_hindcast.jld"); # spatial information
  #const global GRD_Z = GRD["Z"]
  #const global GRD_A = GRD["AREA"]

	#! Integration parameters
	const global DT = 1.; # time step
  const global DAYS = 365; # number of days

	#! Amount of fishing
	if harv == 1
		const global FISHING = 80000000000 / (365/DT) # 80MT per year
	else
		const global FISHING = 0
	end

	#! Benthic-pelagic coupling cutoff (depth, m)
	const global PI_be_cutoff = 200

	#! body lengths (mm)
	const global L_s = 10^((log10(2)+log10(20))/2); # small
	const global L_m = 10^((log10(20)+log10(200))/2); # medium
	const global L_l = 10^((log10(200)+log10(2000))/2); # large

	##! Mass from length using Andersen & Beyer 2013
	# ConverE from mm to cm and use their const coeff = 0.01g/cm3
	const global M_s = 0.01 * (0.1*L_s)^3;
	const global M_m = 0.01 * (0.1*L_m)^3;
	const global M_l = 0.01 * (0.1*L_l)^3;

	#! Mediant Zooplankton size in mm
	# from Charlie/COBALT
	#! Median Zooplankton body mass in g wet weight
	# eq from James Watkins, Lars Rudstam and Kristen Holeck in dry weight
	# convert to wet weight with 0.2g dry = 1 g wet --> 1g dry = 5g wet
	const global L_zm = 10^((log10(0.2)+log10(2))/2); # lengths (ESD)
	const global L_zl = 10^((log10(2)+log10(20))/2);
	const global M_zm = 5.0 * exp(1.953 + (2.399*log(L_zm)))*1.0e-6; # body mass
	const global M_zl = 5.0 * exp(1.953 + (2.399*log(L_zl)))*1.0e-6;

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
	const global fcrit = 0.05	# feeding level needed to meet resting metabolic demands; 0.05-0.2
	const global k = 10.0 		# 10 g^(1-p)/yr at 10C

	###! Consumption constants
	const global h = 85.0  		# h=85 g^(0.25)/yr at 10C in Cmax eq
	# tune so Cobs/Cmax ~ 0.6
	const global flev = 1.0e6	# m^3 g^(−q)/year at 10C; equiv to Andersen, Hartvig gamma
	const global q = 0.8 			# q=0.75-1 in beta eq in consumption

	###! Transfer efficiency of detritus to benthic prey
	const global bent_eff = 0.30

	#! Fraction of time spent swimming (from Van Leeuwen)
	const global Tu_s = 1.0
	const global Tu_m = 1.0 #0.5
	const global Tu_l = 1.0 #0.1

	###! Background mortality
	#Currently increases from 0 to 0.01
	#Megrey et al =0.44/yr
	#Andersen & Beyer 2013 = 0.35 * 4.5 * s^(-0.25) (includes predation, excludes fishing)
	const global Nat_mrt = 0.44 / 365

	###! Diet Preference Phi (j = prey, i = pred)
	# The predator prey mass ratio is assumed 3 orders of mag, i.e. 1000, i.e. one step down
	# We don't have a pred-prey matrix anymore, we are instead explicit about who eats who:
	#-----
	#small forage fish eats medium zoo
	#small piscivores eats medium zoo
	#small detritivore eats medium zoo
	#medium forage fish eats large zoo
	#medium piscivore eats large zoo, small forage fish, small piscvore, small detritivore
	#medium detritivore eats detritus
	#large piscivore eats medium forage fish, medium piscivore, medium detritivore
	#large detritivore eats detritus, medium forage fish, medium piscivore, medium detrivore
	#-----
end
