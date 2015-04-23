
#================= FUNCTIONS FOR SIZE-BASED MODEL ===============#

###! METABOLIC COSTS I as a function of temperature (inline function)
#! inline function - used in the main time loop
#! from Megrey 2007, units: g prey g fish−1 d−1
function fnc_met(s::Float64,T::Float64)
	U = (3.9*s.^0.13)# swimming speed ind of temp
    activity = exp(0.03*U) # activity costs
    temp = exp(0.0548*T) # temperature costs ###<<<<< CORRECT THIS
    basal = 0.0033*s^-0.13 # basal metabolic costs
    #T = basal * temp * activity * 5.258; # total resp costs
    T = basal * temp * 5.258
    return T
end

###! HANDLING TIMES
#! set as 1 / four times metabolic rate
function fnc_tau(I)
	tau = 1 / (4 * I)
	return tau
end

###! CONSUMPTION
function fnc_cons(prey,enc,tot,tau)
	denom = 1 + (tau .* tot)
	numer = prey * enc
	I = numer / denom
	return I
end

###! ENERGY AVAILABLE FOR GROWTH NU
function fnc_nu(I,T,lambda)
	nu = (I*lambda) - T
	if nu < 0.
		nu = 0.
	end
	return nu
end

###! ENERGY AVAILABLE FOR SOMATIC GROWTH
function fnc_gamma(nu,d,k,z)
	numer = (k*nu) - d
	denom = 1 - (z^(1-(d/(k*nu))))
	gamma = numer./denom
	if gamma < 0
		gamma = 0.
	end
	if isnan(gamma)==true
		gamma = 0
	end
	return gamma
end

###! BIOMASS MADE FROM REPRODUCTION
function fnc_rep(nu,k,B)
	rep = (1-k) * nu * B
end











