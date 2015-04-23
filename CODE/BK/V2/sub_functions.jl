
#================= FUNCTIONS FOR SIZE-BASED MODEL ===============#

###! METABOLIC COSTS I as a function of temperature (inline function)
#! inline function - used in the main time loop
#! from Megrey 2007, units: g prey g fish−1 d−1
function fnc_met(s::Float64,U::Float64,T::Float64)
    activity = exp(0.03*U) # activity costs
    temp =  exp(0.0548*T) # temperature costs
    basal = 0.0033*s^-0.13 # basal metabolic costs
    I = basal * temp * activity * 5.258; # total resp costs
    return I
end

###! HANDLING TIMES
#! set as 1 / four times metabolic rate
function fnc_tau(I)
	tau = 1 / (4 * I)
	return tau
end



