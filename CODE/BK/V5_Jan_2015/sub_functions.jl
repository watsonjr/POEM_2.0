
#================= FUNCTIONS FOR SIZE-BASED MODEL ===============#

###! METABOLIC COSTS I as a function of temperature (inline function)
#! inline function - used in the main time loop
#! from Megrey 2007, units: g prey g fish−1 d−1
function fnc_met!(meta::Array{Float64},i::Int,s::Float64,Temp::Float64)
	U = (3.9*s.^0.13)# swimming speed ind of temp
    activity = exp(0.03*U) # activity costs
    temp = exp(0.0548*Temp) # temperature costs ###<<<<< CORRECT THIS
    basal = 0.0033*s^-0.13 # basal metabolic costs
    meta[i] = basal * temp * 5.258
    #meta[i] = basal * temp * activity * 5.258; # total resp costs
end


###! Encounter rates between pairs of species
function fnc_enc!(enc::Array{Float64},i::Int,prey::Array{Float64},phi::Array{Float64},a::Float64)
	enc[:,i] = prey.*phi.*a
end


###! HANDLING TIMES
#! set as 1 / four times metabolic rate
function fnc_tau!(tau::Array{Float64},i::Int,I::Float64)
	tau[i] = 1 / (4 * I)
end


###! CONSUMPTION
function fnc_cons!(I::Array{Float64},d::Array{Float64},i::Int,j::Int,pred::Float64,
				   prey::Float64,enco::Float64,tot::Float64,tau::Float64)
	numer = pred * enco
	denom = 1 + (tau * tot)
	I[i] += numer / denom
	d[j] += numer / denom
end


###! ENERGY AVAILABLE FOR GROWTH NU
function fnc_nu!(nu::Array{Float64},i::Int,bio::Float64,I::Float64,meta::Float64,lambda::Float64)
	nu[i] = (I/bio*lambda) - meta
end


###! ENERGY AVAILABLE FOR SOMATIC GROWTH
function fnc_gamma!(gamma::Array{Float64},i::Int,nu::Float64,
					bio::Float64,d::Float64,k::Float64,z::Float64)
	# note: divide by bio to get biomass specific units
	numer = (k*nu) - (d/bio)
	denom = 1 - (z^(1-((d/bio)/(k*nu))))
	gg = numer./denom
	if gg < 0 || isnan(gg)==true
		gamma[i] = 0.
	else
		gamma[i] = gg
	end
end


###! TOTAL BIOMASS MADE FROM REPRODUCTION
function fnc_rep!(rep::Array{Float64},i::Int,nu::Float64,k::Float64,bio::Float64)
	if nu > 0
		rep[i] = (1-k) * nu * bio
	else 
	    nu = 0.
	end
end


###! TOTAL BIOMASS SOMATIC GROWTH
function fnc_grw!(grw::Array{Float64},i::Int,nu::Float64,bio::Float64)
	grw[i] = nu * bio
end


###! TOTAL BIOMASS MATURING
function fnc_mat!(mat::Array{Float64},i::Int,gam::Float64,bio::Float64)
	mat[i] = gam * bio
end


###! TOTAL BIOMASS dying from background mortality
function fnc_mrt!(mrt::Array{Float64},i::Int,bio::Float64,m::Float64)
	mrt[i] = m * bio
end













