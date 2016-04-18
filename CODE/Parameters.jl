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
	con_f::Array{Float64}  # consumption rate
	con_p::Array{Float64}  # consumption rate
	con_d::Array{Float64}  # consumption rate
	con_zm::Array{Float64} # consumption rate
	con_zl::Array{Float64} # consumption rate
	I::Array{Float64} # total fish ingestion g d-1
	nu::Array{Float64}  # total energy for growth g g-1 d-1
	gamma::Array{Float64} # energy for somatic growth g g-1 d-1
	die::Array{Float64} # total death g d-1
	rep::Array{Float64} # total biomass reproduced g d-1
	rec::Array{Float64} # total biomass reproduced g d-1
end

type detritus
	bio::Array{Float64} # biomass in detrital pool g m-2
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
	#K::Array{Float64} #spawning flag
end

#============= PARAMETER TYPE ==========#
function make_parameters(harv)
	#! Input it switch for fishing (1 is yes, 0 is no)

	#! Grid parameters
	global GRD = load("./Data/Data_grid.jld"); # spatial information
  const global GRD_Z = GRD["Z"]
  const global GRD_A = GRD["AREA"]

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
	const global PI_be_cutoff = 500

	#! body lengths (mm)
	const global L_s = 10^((log10(2)+log10(20))/2); # small
	const global L_m = 10^((log10(20)+log10(200))/2); # medium
	const global L_l = 10^((log10(200)+log10(2000))/2); # large

	##! Mass from length using Andersen & Beyer 2013
	# ConverE from mm to cm and use their const coeff = 0.01g/cm3
	const global M_s = 0.01 * (0.1*L_s)^3;
	const global M_m = 0.01 * (0.1*L_m)^3;
	const global M_l = 0.01 * (0.1*L_l)^3;

	#! Median Zooplankton body mass
	# James Watkins, Lars Rudstam and Kristen Holeck
	# (from Charlie/COBALT)
	const global L_zm = 10^((log10(0.2)+log10(2))/2); # lengths (ESD)
	const global L_zl = 10^((log10(2)+log10(20))/2);
	const global M_zm = exp(1.953 + (2.399*log(L_zm)))*1.0e-6; # body mass
	const global M_zl = exp(1.953 + (2.399*log(L_zl)))*1.0e-6;

	#! Ratio of initial and final body sizes per size-class
	const global Z_s = (0.01*(0.1*20)^3) / (0.01*(0.1*2)^3)
	const global Z_m = (0.01*(0.1*200)^3) / (0.01*(0.1*20)^3)
	const global Z_l = (0.01*(0.1*2000)^3) / (0.01*(0.1*200)^3)

	###! Assimilation efficiency lambda (constant across everything)
	const global Lambda = 0.7;

	###! Kappa rule K as a function of body size
	# K = fraction of energy consumed diverted to somatic growth
	const global K_l = 1
	const global K_j = 0.5
	const global K_a = 0

	###! Metabolism constants (activity and basal)
	#const global Act_s = exp(0.03*(3.9*M_s.^0.13))
	#const global Act_m = exp(0.03*(3.9*M_m.^0.13))
	#const global Act_l = exp(0.03*(3.9*M_l.^0.13))
	const global Bas_s = 0.0033*M_s.^-0.13
	const global Bas_m = 0.0033*M_m.^-0.13
	const global Bas_l = 0.0033*M_l.^-0.13

	###! Swimming speed (Megrey 2007): m d-1
	# Q. Why is swim speed a fraction from 0.8 to 0.1 as gets bigger?
	# Ans: time spent swimming according to Van Leeuven 2008
	const global U_s = ((3.9*M_s.^0.13)/100*60*60*24)
	const global U_m = ((3.9*M_m.^0.13)/100*60*60*24)
	const global U_l = ((3.9*M_l.^0.13)/100*60*60*24)

  	# With temp-dep
	#const global PI_U = ((3.9*M_s.^0.13 * exp(0.149*T)) /100*60*60*24)
	#const global PL_U = ((3.9*M_m.^0.13 * exp(0.149*T)) /100*60*60*24)
	#const global DE_U = ((3.9*M_l.^0.13 * exp(0.149*T)) /100*60*60*24)

	#! Fraction of time spent swimming (from Van Leeuwen)
	const global Tu_s = 1.0
	const global Tu_m = 0.5
	const global Tu_l = 0.1

	#! Maximum search rate a as a function of body length (m2 d-1 g-1)
	# calculate swimming speed (m d-1)
	# factor in time spent swimming (low - Anieke's paper)
	# orgs can see x3 body length (just a guess, find better number)
	# divide by body mass to get weight specific search rate
	const global A_s = (U_s * ((L_s/1000)*3) * Tu_s) / M_s
	const global A_m = (U_m * ((L_m/1000)*3) * Tu_m) / M_m
	const global A_l = (U_l * ((L_l/1000)*3) * Tu_l) / M_l

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
	#small detrivore eats detritus
	#medium forage fish eats large zoo
	#medium piscivore eats large zoo, small forage fish, small piscvore
	#medium detrivore eats small detrivore
	#large piscivore eats medium forage fish, medium piscivore, medium detrivore
	#-----
end
