###! Get COBALT data
function get_COBALT!(COBALT,ID,DY,ENVR,Tref,Dthresh)
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
    ENVR.T0[:,1] = Tref[ID]
    #ENVR.T0[:,1] = minimum(COBALT["Tp"][ID,:])
    ENVR.Dthresh[:,1] = Dthresh[ID]
    ENVR.fZm[:,1] = zeros(Int64,NX)
    ENVR.fZl[:,1] = zeros(Int64,NX)
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


###! Max consumption
#function sub_cmax(Tp,Tb,tdif,wgt)
#  #Tp: pelagic temp
#  #Tb: bottom temp
#  #tdif: frac pelagic time
#  #wgt: ind weight of size class
#  temp = (Tp.*tdif) + (Tb.*(1.0-tdif))
#	 cmax = exp(0.063*(temp-15.0)) * h * wgt^(3/4) #h value for temp=15C
#  return cmax
#end

###! Swimming speed
#function sub_swim(Tp,Tb,tdif,wgt)
#  #Tp: pelagic temp
#  #Tb: bottom temp
#  #tdif: frac pelagic time
#  #wgt: ind weight of size class
#  T = (Tp.*tdif) + (Tb.*(1.0-tdif))
#  U = ((3.9*wgt.^0.13 * exp(0.149*T)) /100*60*60*24)
#  return U
#end


###! Metabolism
function sub_met(Tp,Tb,tdif,wgt,L)
  #Tp: pelagic temp
  #Tb: bottom temp
  #tdif: frac pelagic time
  #wgt: ind weight of size class
	#fcrit: feeding level to meet resting respiration rate
  #cmax: max consumption rate
  #U: swimming speed

  #Cmax
  temp = (Tp.*tdif) + (Tb.*(1.0-tdif))
  cmax = exp(0.063*(temp-15.0)) * h * wgt^(3/4) #h value for temp=15C
  #Swimming
  # NOTE: may want a new temp-dep U equation; I think too fast for larvae 8615m/d
  #U = ((3.9*wgt.^0.13 * exp(0.149*temp)) /100*60*60*24) #Megrey (m/d)
  #U = (3.9*wgt.^0.13 * exp(0.149*temp))  #Megrey (cm/s)
  U = ((L * exp(0.063*(temp-18.5))) /10) #My Q10 using Tref=median hist temp (cm/s)
  #Metabolism
	bas = fcrit * cmax
  #Megrey respiration rate (gprey/gfish/d) uses swim in (cm/s)
  met = bas * exp(0.03*U) * 5.258
  return met
end


###!  Encounter rates
function sub_enc(Tp,Tb,wgt,L,tu,pred,prey,tpel,tprey)
  # Tp: pelagic temp
  # Tb: bottom temp
  # wgt: ind weight of size class
  # L: length of size class
  # tu: frac time spent searching for food
  # pred: pred biomass density,
	# prey: prey biomass density,
	# A: predator search rate,
  # tpel: time spent in pelagic,
	# tprey: time spent in area with that prey item.
  #Swimming
  temp = (Tp.*tpel) + (Tb.*(1.0-tpel))
  #U = ((3.9*wgt.^0.13 * exp(0.149*temp)) /100*60*60*24) #Megrey
  U = ((L * exp(0.063*(temp-18.5))) /1000*60*60*24) #My Q10 using Tref=median hist temp
  #Search rate
  A = (U * ((L/1000)*3) * tu) / wgt
  #Encounter
	enc = pred*prey*A*tprey
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
  cmax = exp(0.063*(temp-15.0)) * h * wgt^(3/4) #h value for temp=15C
  #Con
	beta = flev * exp(0.063*temp-15.0) * wgt^(q)
  #beta = 1.0;
  ENC = sum(enc) # total biomass encountered
	con = cmax .* (beta .* enc[1]) ./ (cmax + beta.*ENC) # Type II
  return con
end


###! Offline coupling
function sub_offline(enc_1,enc_2,enc_3,dZ)
  # ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
	#! offline switch
	if (enc_1 + enc_2 + enc_3) > dZ
		frac1 = enc_1 / (enc_1 + enc_2 + enc_3)
    frac2 = enc_2 / (enc_1 + enc_2 + enc_3)
    frac3 = enc_3 / (enc_1 + enc_2 + enc_3)
		out_1 = (frac1 * dZ)
		out_2 = (frac2 * dZ)
    out_3 = (frac3 * dZ)
    zf = 1
	else
		out_1 = enc_1
		out_2 = enc_2
    out_3 = enc_3
    zf = 0
	end
	return out_1, out_2, out_3, zf
end

function sub_offline_bent(enc_1,enc_2,B)
  if (enc_1 + enc_2) > B
		frac1 = enc_1 / (enc_1 + enc_2)
    frac2 = enc_2 / (enc_1 + enc_2)
    out_1 = (frac1 * B)
		out_2 = (frac2 * B)
	else
		out_1 = enc_1
		out_2 = enc_2
	end
	return out_1, out_2
end


###! Consumption/Cmax
function sub_clev(con,Tp,Tb,tdif,wgt)
	#! calculates consumption rate of first element of enc
  #Cmax
  temp = (Tp.*tdif) + (Tb.*(1.0-tdif))
  cmax = exp(0.063*(temp-15.0)) * h * wgt^(3/4) #h value for temp=15C
  #clev
	clev = con/cmax
  return clev
end


###! DEGREE DAYS
function sub_degday(dd,Tp,Tb,tdif,Tref,S,dtot)
  #if (S==0.0) #Don't accumulate temp while spawning, DD represents recovery after
  #if (sum(S[1:dtot]) < dtot) #Only spawn once per year
  if (sum(S[1:dtot]) > 0.0) #Only spawn once per year
    dd = 0.0
  else
    Tavg = (Tp.*tdif) + (Tb.*(1.0-tdif))
    dd += max((Tavg-Tref),0.0)
  end
  return dd
end


###! SPAWNING FLAG
function sub_kflag(S,dd,dthresh,dtot)
  if (dd >= dthresh)
    dur=59
    #Change spawning flag
    if ((dtot+dur) <= DAYS)
      S[dtot:(dtot+dur)] = Sp
    else
      dleft = DAYS - dtot + 1
      S[dtot:DAYS] = Sp[1:dleft]
    end
    #Reset cumulative deg days
    dd = 0.0
  end
  return S, dd
end


###! ENERGY AVAILABLE FOR GROWTH NU
function sub_nu(I,B,met)
	# convert to biomass specific ingestion
	nu = ((I/B)*Lambda) - met
end


###! ENERGY AVAILABLE FOR SOMATIC GROWTH
function sub_gamma(K,Z,nu,d,B,S)
  # convert predation mortality to biomass specific rate
	D = (d/B) + Nat_mrt
  # Spawning flag
  #if S>0.0
  #  kap=min(1.0, K + (1.0-S));
  #else
  #  kap=1;
  #end
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
function sub_rec(X,bio)
	# X could be biomass of eggs (for larval class) or maturing from smaller sizes
	rec = 0.0
	for i = 1:length(X)
		rec += X[i] * bio[i]
	end
	return rec
end


###! Update biomass
function sub_update_fi(bio_in,rec,nu,rep,gamma,die,egg)
	# all inputs except rec are in g g-1 d-1; rec is g d-1
	# rec = rep from smaller size class = TOTAL biomass gained from recruitment
	# grw = nu = somatic energy for growth within size class
  # store = egg = energy stored for later egg production
	# rep = rep =  energy lost to egg production
	# mat = gamma = energy lost to maturation to larger size class
	# Nat_mrt = natural mortality
	# die = predator mort = biomass lost to predation
  db = rec + ((nu - egg - rep - gamma - Nat_mrt) * bio_in) + (egg * bio_in) - die
  bio_out =  bio_in + db
end

function sub_update_be(bio_in,die)
  bio_out = bio_in - sum(die)
end


####! Fishing
#function sub_fishing(bio_pi,bio_pl,bio_de,AREA)
#	if FISHING > 0.0
#		#bio_pi = PISC.bio; bio_pl = PLAN.bio; bio_de = DETR.bio; AREA = GRD_A;
#		ALL_pi  = Array(Float64,NX,PI_N)
#		ALL_pl  = Array(Float64,NX,PL_N)
#		ALL_de  = Array(Float64,NX,DE_N)
#
#		for i = 1:NX
#			ALL_pi[i,:] = bio_pi[i] * AREA[i]
#			ALL_pl[i,:] = bio_pl[i] * AREA[i]
#			ALL_de[i,:] = bio_de[i] * AREA[i]
#		end
#
#		#! Total fish biomass
#		TOT = sum(ALL_pi) + sum(ALL_pl) + sum(ALL_de)
#		ALL_pi -= (ALL_pi./TOT).*FISHING
#		ALL_pl -= (ALL_pl./TOT).*FISHING
#		ALL_de -= (ALL_de./TOT).*FISHING
#
#		#! Calc total biomass of fish in the ocean
#		for i = 1:NX
#			bio_pi[i] = squeeze(ALL_pi[i,:],1) ./ AREA[i]
#			bio_pl[i] = squeeze(ALL_pl[i,:],1) ./ AREA[i]
#			bio_de[i] = squeeze(ALL_de[i,:],1) ./ AREA[i]
#		end
#	end
#	return bio_pi, bio_pl, bio_de
#end


###! Forward Euler checks
function sub_check!(bio)
	ID = find(bio .< 0)
	bio[ID] = eps()
end
