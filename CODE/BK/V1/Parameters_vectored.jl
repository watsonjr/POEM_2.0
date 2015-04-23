#============== Parameters of the model =============#
## Dictionaries for parameters of the major groups
#! N: number of size classes (#)
#! s_min: minimum size (g)
#! s_max: maximum size (g)
#! s: median size (g) of each size-class
#! z: ratio of initial and final body sizes for each size class
#! T: metabolic costs (g_i g_j^-1 d^01)
#! lambda: assimilation efficiency (g_i g_j^-1)
#! K: kappa, fraction of energy given to somatic growth
#! Phi: diet preference (probability)
#! a: search rate (m2 d^-1 g_i^-1)
#! tau: handling time (d g_j^-1 g_i)


##### Setup variable types
type BIOMASS{N}
	B::Array{Float64,N}
end

PISC = BIOMASS(P_PI["N"]);

##### Setup Parameter Dictionaries
P_DE = (ASCIIString=>Any)["N"=>Int[],"s_min"=>Float64[],
			    "s_max"=>Float64[],"s"=>Float64[],"z"=>Float64[],
			    "T"=>Float64[],"lambda"=>Float64[],"K"=>Float64[],
			    "Phi"=>Float64[],"a"=>Float64[],"tau"=>Float64[]];
				
P_PI = (ASCIIString=>Any)["N"=>Int[],"s_min"=>Float64[],
			    "s_max"=>Float64[],"s"=>Float64[],"z"=>Float64[],
			    "T"=>Float64[],"lambda"=>Float64[],"K"=>Float64[],
			    "Phi"=>Float64[],"a"=>Float64[],"tau"=>Float64[]];
				
P_PL = (ASCIIString=>Any)["N"=>Int[],"s_min"=>Float64[],
			    "s_max"=>Float64[],"s"=>Float64[],"z"=>Float64[],
			    "T"=>Float64[],"lambda"=>Float64[],"K"=>Float64[],
			    "Phi"=>Float64[],"a"=>Float64[],"tau"=>Float64[]];
				
##### Fill in Parameter Dictionaries
#! Number of size classes
P_DE["N"] = 5;
P_PI["N"] = 10;
P_PL["N"] = 10;

#! Min body size 
P_DE["s_min"] = 1;
P_PI["s_min"] = 1;
P_PL["s_min"] = 1;

#! Max body size 
P_DE["s_max"] = 1000;
P_PI["s_max"] = 200;
P_PL["s_max"] = 200;

#! Body sizes
P_DE["s"] = 10.^(linspace(log10(P_DE["s_min"]),log10(P_DE["s_max"]),P_DE["N"]));
P_PI["s"] = 10.^(linspace(log10(P_PI["s_min"]),log10(P_PI["s_max"]),P_PI["N"]));
P_PL["s"] = 10.^(linspace(log10(P_PL["s_min"]),log10(P_PL["s_max"]),P_PL["N"]));

#! Ratio of initial and final body sizes per size-class
### HERE
s = log10(P_DE["s"]);
ds = diff(s)/2

for i = 1:5
	si = s[i]
	zmn = 10.^(si


P_DE["z"] = ???


#! Metabolic costs T


#! Assimilation efficiency lambda


#! Kappa rule K


#! Diet Preference Phi


#! Search rate a


#! Handling Times tau









# parameters.jl
module Parameters

##### Number of size classes
P_CN  = 5; # number of cod size classes
P_FN  = 5; # number of herring size classes
P_BpN = 1; # number of detrital feeders
P_BdN = 3; # number of benthic pred
P_ZN  = 2; # number of zooplankton size-classes (from COBALT)

###### size classes (log g)
P_Cs  = linspace(1.0e2,1.0e3,P_CN); # cod size classes
P_Fs  = linspace(0.007e3,0.013e3,P_FN); # cod size classes
P_Bds = linspace(0.0002e3,0.001e3,P_BdN); # cod size classes
P_Bps = linspace(0.0002e3,0.001e3,P_BpN); # cod size classes
P_Zs  = [128e-6,70000e-6];
const sizes = [P_Zs,P_Cs,P_Fs,P_Bds,P_Bps]


#! types of organisms for feeding
# 0=zoo, 1=plankt., 2=pisc., 3=detrital fed., 4=benthic pred
orgTypes = [zeros(P_ZN); ones(P_CN)*1; ones(P_FN)*2; ones(P_BdN)*3;
			ones(P_BpN)*4];

############ GENERAL parameters
const P_noGroups = length(orgTypes) # total number of groups
const P_Tmax = 360 # total simulation time (days)
const P_current_time = 1 # start day
const P_dt = 1 # time step in days


############ SIZE SPECIFIC PARAMETERS
#! Allometric constants
const P_E = 0.69 # (eV) Energy of activations
const P_k = 8.62e-5 # (eV K-1) Boltzmann constant
const P_kau = 0.626 # 0.25-1 in Watson et al. 2014

#! Swimming speeds (figure this out better: a = 100 m d-1, b= 0.13)
P_V = 100 .* sizes .^ 0.13;
P_V[orgTypes.==0] = NaN; # zoo offline, don't need speed
P_V[orgTypes.==3] = NaN; # det. feeds don't swim
const P_V

#! Basal metabolic rates (Brown et al. 2004) and multipliers
#! could use Megrey as well
#! variable in time as a function of temp T [Kelvin] (and/or O2)
#! function = exp( (P_Meta*log(size) + P_Metb) - E/(k*t))
P_Meta = 0.71;
P_Metb = 18.47;

#! Natural Mortality (Brown et al. 2004)
#! variable in time as a function of temp T [Kelvin]
#! function = exp( (P_Morta*log(size) + P_Mortb) - E/(k*t))
P_Morta = -0.24;
P_Mortb = 26.25

#! Diet Preferences
P_PPMR  = 3. # optimum prey size - Barnes et al. 2010
P_stdev = 1.3 # standard dev.

#! Somatic Growth parameters (fraction of sizeclass allocated to somatic growth)
#! according to de Roos et al.(2008)
P_beta = [


const kP3 = 0.# proportion of production into somatic growth (no unit)-vanLeeuwen et al.(2008)
const kP2 = 0.8
const kP1 = 1.
const kF2 = 0.4 # in van Leeuwen et al.(2008) S2=0.8, S3=0
const kF1 = 1.

# for maturation calc
const matP2 = 0.125 # initial size / final size in age group (no unit)
const matP1 = 0.003 # values from van Leeuwen et al.(2008)
const matF1 = 0.05 # S1

const matAll = [0, matP2, matP1, 0, matF1, 0, 0, 0, 0, 0, 0]

# maturation is not a constant (van Leeuwen et al.(2008)); calc in code

const recP3 = 1 - kP3 # recruitment into juveniles/larvae as proportion of B produced (no unit)
const recP2 = 1 - kP2 
const recF2 = 1 - kF2 

const reprod = [recP3, recP2, 0, recF2, 0, 0, 0, 0, 0, 0, 0] # reproduction vector

const lambdaFish = 0.7 # assimilation efficiency for fish (no unit)

const lP3 = ((sP3/0.025)^(1/3))/100 # length (m); Watson et al.
const lP2 = ((sP2/0.025)^(1/3))/100
const lP1 = ((sP1/0.025)^(1/3))/100
const lF2 = ((sF2/0.025)^(1/3))/100
const lF1 = ((sF1/0.025)^(1/3))/100

const ir = 0.1 # irradiance const (no unit); 0.05-0.25 in Watson et al.
const K = 0.5 # constant (no unit); T dependant in Watson et al.
const a = 864. # (m d-1)
const b = 0.13 # (no unit)

const vP3 = K*a*sP3^b # swimming speed (m d-1 g-0.33)
const vP2 = K*a*sP2^b
const vP1 = K*a*sP1^b
const vF2 = K*a*sF2^b
const vF1 = K*a*sF1^b

const alphaP3 = (ir*pi*(lP3^2)*vP3)/sP3 # search volume (m3 d-1 g-1)
const alphaP2 = (ir*pi*(lP2^2)*vP2)/sP2
const alphaP1 = (ir*pi*(lP1^2)*vP1)/sP1
const alphaF2 = (ir*pi*(lF2^2)*vF2)/sF2
const alphaF1 = (ir*pi*(lF1^2)*vF1)/sF1

const alphaAll = [alphaP3, alphaP2, alphaP1, alphaF2, alphaF1, 0, 0, 0, 0, 0, 0]

const t0 = 93.69 # (d); for handling time

# handling time eq. t0*si^83*sj^0.5 (need to know who eats whom)

# parameters for SIZE-BASED functions

const lambdaDet = 0.7 # assimilation efficiency for detrivores (no unit) -Timmermann et al.(2012)
const lambdaPred = 0.7 # assimilation efficiency for predators (no unit) -Timmermann et al.(2012)

const imaxB2 = 0.99e-3 # maximum specific ingestion rate (g g-1 d-1)-Timmermann et al.(2012)
const imaxB13 = 2.5e-3
const imaxB12 = 2.5e-3
const imaxB11 = 2.5e-3
const imaxZ2 = (1.68*sZ2^0.703)/sZ2 # following Saiz&Calbet (2007); orig units g C ind-1 d-1 
const imaxZ1 = (1.68*sZ1^0.703)/sZ1

const imax = [0, 0, 0, 0, 0, imaxB2, imaxB13, imaxB12, imaxB11, imaxZ2, imaxZ1]

const KB2 = 75. # half-saturation constant (g m-3); Hummel 1985
const KB13 = 75.
const KB12 = 75.
const KB11 = 75.
const KZ2 = 4.8 # hansen et al.1997 - note:need to CHECK!!! (size-indep)
const KZ1 = 4.8 

const KAll = [0, 0, 0, 0, 0, KB2, KB13, KB12, KB11, KZ2, KZ1]

# the additional metabolic loss due to env. not written in yet

end
