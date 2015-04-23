# parameters.jl

module Parameters

export PP, DET # resources

export noGroups, P_Tmax, current_time, dt, Binit, sizes, orgTypes # general
export m0All, mAll
export PPMR, stdev

export reprod, matAll # size-structured
export lambdaFish, t0
export alphaAll

export lambdaDet, lambdaPred # size-based
export imax, KAll

# RESOURCES FROM BGM
const PP = 25. # primary production (g m-2 d-1)
const DET = 30. # detritus flux (g m-2 d-1)

# GENERAL parameters

const noGroups = 11 # total number of groups
const P_Tmax = 360 # total simulation time
const current_time = 1
const dt = 1

const BinitP3 = 1.# initial biomass (g)
const BinitP2 = 2.
const BinitP1 = 3.
const BinitF2 = 10.
const BinitF1 = 14.
const BinitB2 = 12.
const BinitB13 = 18.
const BinitB12 = 24.
const BinitB11 = 30.
const BinitZ2 = 30.
const BinitZ1 = 30.

const Binit = [BinitP3, BinitP2, BinitP1, BinitF2, BinitF1, BinitB2, BinitB13, 
			   BinitB12, BinitB11, BinitZ2, BinitZ1]

const orgTypes = [1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 5] # types of organisms for feeding
# 1=pisc., 2=plankt., 3=bent.pred., 4=bent.detr.feeder, 5=zpl

const sP3 = 1.0e3 # size of the largest predator group (g);cod23 ICES WGBFAS'13
const sP2 = 0.26e3
const sP1 = 1.0e2 # NO DATA IN
const sF2 = 0.013e3 # 2+ sprat, 3+ herring mean weight (weighted)
const sF1 = 0.007e3 # 1 sprat, 1-2 herring mean weight (weighted)
const sB2 = 0.0017e3 # S. entomon
const sB13 = 0.001e3 # FIMR data
const sB12 = 0.0007e3
const sB11 = 0.0002e3
const sZ2 = 1.05e-4 # P. acuspes (Renz et al.2007)
const sZ1 = 1.9e-5 # Acartia spp. (Dzierbicka-Glowacka et al.2009)

const sizes = [sP3, sP2, sP1, sF2, sF1, sB2, sB13, sB12, sB11, sZ2, sZ1]
# max swimming speed ?

const E = 0.69 # (eV)
const k = 8.62e-5 # (eV K-1)
const temp = 17 # (C)
const T = temp + 273.15 # (K)
const kau = 0.626 # 0.25-1 in Watson et al.
const ir = 0.15 # 0.05-0.25 in Watson et al.
const ab = 864 # (m d-1 g-0.33) - a in Watson et al.
const b = 0.13

const vP3 = kau*ab*sP3^b
const vP2 = kau*ab*sP2^b
const vP1 = kau*ab*sP1^b
const vF2 = kau*ab*sF2^b
const vF1 = kau*ab*sF1^b
const vB2 = kau*ab*sB2^b
const vB13 = kau*ab*sB13^b
const vB12 = kau*ab*sB12^b
const vB11 = kau*ab*sB11^b
const vZ2 = kau*ab*sZ2^b
const vZ1 = kau*ab*sZ1^b

# in case of const Temperature
const mP3 = exp((0.71*log(sP3)+18.47)-E/(k*T))
const mP2 = exp((0.71*log(sP2)+18.47)-E/(k*T))
const mP1 = exp((0.71*log(sP1)+18.47)-E/(k*T))
const mF2 = exp((0.71*log(sF2)+18.47)-E/(k*T))
const mF1 = exp((0.71*log(sF1)+18.47)-E/(k*T))
const mB2 = exp((0.71*log(sB2)+18.47)-E/(k*T))
const mB13 = exp((0.71*log(sB13)+18.47)-E/(k*T))
const mB12 = exp((0.71*log(sB12)+18.47)-E/(k*T))
const mB11 = exp((0.71*log(sB11)+18.47)-E/(k*T))
const mZ2 = exp((0.71*log(sZ2)+18.47)-E/(k*T))
const mZ1 = exp((0.71*log(sZ1)+18.47)-E/(k*T))

const mAll = [mP3, mP2, mP1, mF2, mF1, mB2, mB13, mB12, mB11, mZ2, mZ1]

const m0P3 = exp(-0.24*log(sP3)+26.25)/exp(E/(k*T))  #other nat. mortality according to Watson et al.(in press) (no unit)
const m0P2 = exp(-0.24*log(sP2)+26.25)/exp(E/(k*T))
const m0P1 = exp(-0.24*log(sP1)+26.25)/exp(E/(k*T))
const m0F2 = exp(-0.24*log(sF2)+26.25)/exp(E/(k*T))
const m0F1 = exp(-0.24*log(sF1)+26.25)/exp(E/(k*T))
const m0B2 = exp(-0.24*log(sB2)+26.25)/exp(E/(k*T))
const m0B13 = exp(-0.24*log(sB13)+26.25)/exp(E/(k*T))
const m0B12 = exp(-0.24*log(sB12)+26.25)/exp(E/(k*T))
const m0B11 = exp(-0.24*log(sB11)+26.25)/exp(E/(k*T))
const m0Z2 = exp(-0.24*log(sZ2)+26.25)/exp(E/(k*T))
const m0Z1 = exp(-0.24*log(sZ1)+26.25)/exp(E/(k*T))

const m0All = [m0P3, m0P2, m0P1, m0F2, m0F1, m0B2, m0B13, m0B12, m0B11, m0Z2, m0Z1]

const PPMR = 3. # optimum prey size - Watson et al.
const stdev = 1.3 # standard dev. - Watson et al.

# parameters for SIZE-STRUCTURED functions
# according to de Roos et al.(2008)
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
