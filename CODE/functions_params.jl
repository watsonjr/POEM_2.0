###! METABOLISM
##! Temperature multiplier for metabolism
tmet = exp(0.0*T)
#tmet = exp(0.0548.*T)
##! Activity multiplier for metabolism
#Pisc
umet[i] = exp(0.03*(PI_U[i]*100/60/60/24))
#Plank
umet[i] = exp(0.03*(PL_U[i]*100/60/60/24))
#Detr
umet[i] = exp(0.03*(DE_U[i]*100/60/60/24))
##! Metabolism
#Pisc
met[i] = PI_bas[i] * tmet * umet[i] * 5.258
#Plank
met[i] = PL_bas[i] * tmet * umet[i] * 5.258
#Detr
met[i] = DE_bas[i] * tmet * umet[i] * 5.258

###! INGESTION
##! Fraction of time spent in pelagic (for piscivore)
if Z < PI_be_cutoff
	tdif = sum(bio1) ./ (sum(bio1).+sum(bio2))
else
	tdif = 1.0
end
##! Handling times
#Pisc
tau[i] = 1 / (4*met[i])
#Plank
tau[i] = 1 / (4*met[i])
#Detr
tau[i] = 1 / (4*met[i])
##! Piscivore encounter rates
#Pisc on Pisc
enc[j,i] = prey[j]*PI_phi_PI[j,i]*PI_a[i]*tdif
#Pisc on Plank
enc[j,i] = prey[j]*PI_phi_PL[j,i]*PI_a[i]*tdif
#Pisc on Detr
enc[j,i] = prey[j]*PI_phi_DE[j,i]*PI_a[i]*(1-tdif)
#Pisc on Zoo
enc[1,i] = Zm*PI_phi_Z[1,i]*PI_a[i]*tdif
enc[2,i] = Zl*PI_phi_Z[2,i]*PI_a[i]*tdif
##! Planktivore encounter rates
enc[1,i] = Zm*PL_phi_Z[1,i]*PL_a[i]
enc[2,i] = Zl*PL_phi_Z[2,i]*PL_a[i]
##! Detrivore encounter rates
#Detr on Detr
enc[j,i] = prey[j]*DE_phi_DE[j,i]*DE_a[i]
#Detr on Benthic material
enc[i] = prey[1]*DE_phi_BE[i]*DE_a[i]
##! Total biomass encountered
#Pisc on Pisc
ENC[i] += enc_pi[j,i]
#Pisc on Plank
ENC[i] += enc_pl[j,i]
#Pisc on Detr
ENC[i] += enc_de[j,i]
#Pisc on Zoo
ENC[i] += enc_z[j,i]
#Detr on Detr
ENC[i] += enc_de[j,i]
#Detr on Benth
ENC[i] += enc_be[i]
#Plank on Zoo
ENC[i] = 0.0 #???
ENC[i] += enc_z[1,i]
ENC[i] += enc_z[2,i]
##! type II feeding
function sub_typeII(pred,enc,tau,ENC)
	return (pred*enc) / (1 + (tau*ENC))
end
##! Piscivore consumption
#Pisc on Pisc
con = sub_typeII(bio_pi[i],enc_pipi[j,i],tau_pi[i],ENC_pi[i])
#Pisc on Plank
con = sub_typeII(bio_pi[i],enc_pipl[j,i],tau_pi[i],ENC_pi[i])
#Pisc on Detr
con = sub_typeII(bio_pi[i],enc_pide[j,i],tau_pi[i],ENC_pi[i])
##! Detrivore consumption
#Detr on Detr
con = sub_typeII(bio_de[i],enc_dede[j,i],tau_de[i],ENC_de[i])
#Detr on Benth
con = sub_typeII(bio_de[i],enc_debe[1,i],tau_de[i],ENC_de[i])
##! Consumption of zooplankton
#! by Piscivore
I_z[1,i] += sub_typeII(bio_pi[i],enc_piz[1,i],tau_pi[i],ENC_pi[i])
I_z[2,i] += sub_typeII(bio_pi[i],enc_piz[2,i],tau_pi[i],ENC_pi[i])
#! by Planktivore
I_z[1,i] += sub_typeII(bio_pl[i],enc_plz[1,i],tau_pl[i],ENC_pl[i])
I_z[2,i] += sub_typeII(bio_pl[i],enc_plz[2,i],tau_pl[i],ENC_pl[i])
##! Can't consume more zooplankton than COBALT mortality rate
ZOO = sum(pi_z,2) + sum(pl_z,2)
if ZOO[1] > dZm
	pi_z[1,:] = (pi_z[1,:]./ZOO[1]) .* dZm
	pl_z[1,:] = (pl_z[1,:]./ZOO[1]) .* dZm
end
if ZOO[2] > dZl
	pi_z[2,:] = (pi_z[2,:]./ZOO[2]) .* dZl
	pl_z[2,:] = (pl_z[2,:]./ZOO[2]) .* dZl
end

###! GROWTH & REPRODUCTION
##! ENERGY AVAILABLE FOR GROWTH (NU)
# Piscivore
nu[i] = ((I[i]/bio[i])*PI_lambda[i]) - met[i]
# Planktivore
nu[i] = ((I[i]/bio[i])*PL_lambda[i]) - met[i]
# Detr
nu[i] = ((I[i]/bio[i])*DE_lambda[i]) - met[i]
##! ENERGY AVAILABLE FOR SOMATIC GROWTH
# note: divide by bio to get biomass specific units
# Spawning flag
#Piscivore
if k==1
  kap=PI_K[i];
else
  kap=1;
end
gg = ((kap*nu[i]) - (d[i]/bio[i]))/(1-(PI_z[i]^(1-((d[i]/bio[i])/(kap*nu[i])))))
if gg < 0 || isnan(gg)==true
  gamma[i] = 0.0
else
  gamma[i] = gg
end
#Plank
if k==1
  kap=PL_K[i];
else
  kap=1;
end
gg = ((kap*nu[i]) - (d[i]/bio[i]))/(1-(PL_z[i]^(1-((d[i]/bio[i])/(kap*nu[i])))))
if gg < 0 || isnan(gg)==true
  gamma[i] = 0.0
else
  gamma[i] = gg
end
#Detr
if k==1
  kap=DE_K[i];
else
  kap=1;
end
gg = ((kap*nu[i]) - (d[i]/bio[i]))/(1-(DE_z[i]^(1-((d[i]/bio[i])/(kap*nu[i])))))
if gg < 0 || isnan(gg)==true
  gamma[i] = 0.0
else
  gamma[i] = gg
end
##! TOTAL BIOMASS MADE FROM REPRODUCTION
# Spawning flag
#Pisc
if nu[i] > 0.
  if k==1,
    kap=PI_K[i];
  else,
    kap=1;
  end
  rep[i] = (1-kap) * nu[i] * bio[i]
else
  rep[i] = 0.
end
#Plank
if nu[i] > 0.
  if k==1,
    kap=PL_K[i];
  else,
    kap=1;
  end
  rep[i] = (1-kap) * nu[i] * bio[i]
else
  rep[i] = 0.
end
#Detr
if nu[i] > 0.
  if k==1,
    kap=DE_K[i];
  else,
    kap=1;
  end
  rep[i] = (1-kap) * nu[i] * bio[i]
else
  rep[i] = 0.
end
##! TOTAL BIOMASS SOMATIC GROWTH
grw[i] = nu[i] * bio[i]
##! TOTAL BIOMASS MATURING
mat[i] = gam[i] * bio[i]
##! TOTAL BIOMASS dying from background mortality
#Pisc
mrt[i] = PI_mrt[i] * bio[i]
#Plank
mrt[i] = PL_mrt[i] * bio[i]
#Detr
mrt[i] = DE_mrt[i] * bio[i]

###! UPDATE BIOMASS
#Pisc
bio[1] += (sum(rep) + grw[1] - mat[1] - d[1] - mrt[1]) * DT
bio[i] += (mat[i-1] + grw[i] - mrt[i] - rep[i] - mat[i] - d[i]) * DT
#Plank
bio[1] += (sum(rep) + grw[1] - mat[1] - d[1] - mrt[1]) * DT
bio[i] += (mat[i-1] + grw[i] - mrt[i] - rep[i] - mat[i] - d[i]) * DT
#Detr
bio[1] += (sum(rep) + grw[1] - mat[1] - d[1] - mrt[1]) * DT
bio[i] += (mat[i-1] + grw[i] - mrt[i] - rep[i] - mat[i] - d[i]) * DT
#Benth
bio[1] += (det - d[1] - (bio[1]*0.01)) * DT # with sedimentation

###! FISHING
function sub_fishing(bio_pi,bio_pl,bio_de,AREA)
	if FISHING > 0.0
		#bio_pi = PISC.bio; bio_pl = PLAN.bio; bio_de = DETR.bio; AREA = GRD_A;
		ALL_pi  = Array(Float64,NX,PI_N)
		ALL_pl  = Array(Float64,NX,PL_N)
		ALL_de  = Array(Float64,NX,DE_N)

		for i = 1:NX
			ALL_pi[i,:] = bio_pi[i] * AREA[i]
			ALL_pl[i,:] = bio_pl[i] * AREA[i]
			ALL_de[i,:] = bio_de[i] * AREA[i]
		end

		#! Total fish biomass
		TOT = sum(ALL_pi) + sum(ALL_pl) + sum(ALL_de)
		ALL_pi -= (ALL_pi./TOT).*FISHING
		ALL_pl -= (ALL_pl./TOT).*FISHING
		ALL_de -= (ALL_de./TOT).*FISHING

		#! Calc total biomass of fish in the ocean
		for i = 1:NX
			bio_pi[i] = squeeze(ALL_pi[i,:],1) ./ AREA[i]
			bio_pl[i] = squeeze(ALL_pl[i,:],1) ./ AREA[i]
			bio_de[i] = squeeze(ALL_de[i,:],1) ./ AREA[i]
		end
	end
	return bio_pi, bio_pl, bio_de
end


###! Forward Euler checks
function sub_check_pi!(bio)
	for i = 1:PI_N
		if bio[i] <= 0.0 || isnan(bio[i])==1
			bio[i] = eps()
		end
	end
end
function sub_check_pl!(bio)
	for i = 1:PL_N
		if bio[i] <= 0.0 || isnan(bio[i])==1
			bio[i] = eps()
		end
	end
end
function sub_check_de!(bio)
	for i = 1:DE_N
		if bio[i] <= 0.0 || isnan(bio[i])==1
			bio[i] = eps()
		end
	end
end
function sub_check_be!(bio)
	if bio[1] <= 0.0 || isnan(bio[1])==1
		bio[1] = eps()
	end
end
