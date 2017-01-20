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
    #ENVR.U[:,1]   = COBALT["U"][ID,DY]
    #ENVR.V[:,1]   = COBALT["V"][ID,DY]
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
  # use preference for benthic prey in calculation
  #biod = bio3+(LD_phi_BE*bio4)
  biod = bio3+bio4
	if Z < PI_be_cutoff
		tdif = biop ./ (biop+biod)
	else
		tdif = 0.0
	end
	return tdif
end


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
  temp = (Tp.*tdif) + (Tb.*(1.0-tdif))
  #Cmax
  #! Specific ingestion rate from Kiorboe & Hirst (g/g/day)
  #cmax = (exp(0.063*(temp-15.0)) * 10^(0.4) * wgt^(-0.51)) .* 24e-3
  #! Specific ingestion rate from Hartvig et al (g/g/day) ref to 15C
  cmax = (exp(0.063*(temp-15.0)) * 60.0 * wgt^(-0.25)) ./365.0
  #! Specific ingestion rate from Hartvig et al (g/g/day) ref to 10C
  #cmax = (exp(0.063*(temp-10.0)) * 85.0 * wgt^(-0.25)) ./365.0
  #! Specific ingestion rate from mizer (g/g/day) ref to 10C
  #cmax = (exp(0.063*(temp-10.0)) * 40.0 * wgt^(-0.25)) ./365.0
  #! Specific ingestion rate from J&C15 (g/g/day) ref to 10C
  #cmax = (exp(0.063*(temp-10.0)) * 25.0 * wgt^(-0.33)) ./365.0
  #Metabolism
	bas = fcrit * cmax
  met = bas
  #! Hartvig et al (g/g/day) k=10 ref to 10C
  #met = (exp(0.063*(temp-10.0)) * 10 * wgt^(-0.25)) ./365.0
  #! mizer et al (g/g/day) k=4.8 ref to 10C
  #met = (exp(0.063*(temp-10.0)) * 4.8 * wgt^(-0.25)) ./365.0
  #! J&C15 et al (g/g/day) k=2.5 ref to 10C
  #met = (exp(0.063*(temp-10.0)) * 2.5 * wgt^(-0.25)) ./365.0
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
  #! divide by 100m to put in m2/g/day b/c zoop is integrated over 100m depth
  A = (exp(0.063*(temp-15.0)) * 10^(3.24) * wgt^(-0.24)) * (24e-3/9) #/ 100.0
  #! Hartvig et al (m3/g/day) gamma = 0.8e4 ref to 10C
  #A = (exp(0.063*(temp-10.0)) * 0.8e4 * wgt^(q-1.0)) ./365.0
  #! mizer et al (m3/g/day) gamma = 2.9e3? ref to 10C
  #A = (exp(0.063*(temp-10.0)) * 2.9e3 * wgt^(q-1.0)) ./365.0
  #! J&C15 et al (m3/g/day) gamma = 2.9e3? ref to 10C
  #A = (exp(0.063*(temp-10.0)) * 1.37e4 * wgt^(0.9-1.0)) ./365.0
  #Encounter per predator, mult by biomass later
  enc = prey*A*tprey*pref
  return enc
end

###!  Encounter rates for demersal
function sub_enc_pel(Tp,Tb,wgt,pred,prey,tpel,tprey,pref)
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
  #! Specific clearance rates from Blanchard et al. 2008 (m3/g/day)
  A = (exp(0.063*(temp-15.0)) * 640 * wgt^(-0.18)) ./365.0
  #Encounter per predator, mult by biomass later
  enc = prey*A*tprey*pref
  return enc
end

###!  Encounter rates for demersal
function sub_enc_dem(Tp,Tb,wgt,pred,prey,tpel,tprey,pref)
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
  #! Specific clearance rates from Blanchard et al. 2008 (m3/g/day)
  A = (exp(0.063*(temp-15.0)) * 64 * wgt^(-0.25)) ./365.0
  #Encounter per predator, mult by biomass later
  enc = prey*A*tprey*pref
  return enc
end


###! Type II consumption
function sub_cons(Tp,Tb,tpel,wgt,enc)
  #Tp: pelagic temp
  #Tb: bottom temp
  #tpel: frac pelagic time
  #wgt: ind weight of size class
  #enc: array of all encountered food
	#! calculates consumption rate of first element of enc
  #Cmax
  temp = (Tp.*tpel) + (Tb.*(1.0-tpel))
  #! Specific ingestion rate from Hartvig et al (g/g/day) ref to 15C
  cmax = (exp(0.063*(temp-15.0)) * 60.0 * wgt^(-0.25)) ./365.0
  #! Hartvig et al (g/g/day) ref to 10C
  #cmax = (exp(0.063*(temp-10.0)) * 85.0 * wgt^(-0.25)) ./365.0
  #! mizer (g/g/day) ref to 10C
  #cmax = (exp(0.063*(temp-10.0)) * 40.0 * wgt^(-0.25)) ./365.0
  #! J&C15 (g/g/day) ref to 10C
  #cmax = (exp(0.063*(temp-10.0)) * 25.0 * wgt^(-0.33)) ./365.0
  #! Specific ingestion rate from Kiorboe & Hirst (g/g/day)
  #cmax = (exp(0.063*(temp-15.0)) * 10^(0.4) * wgt^(-0.51)) .* 24e-3
  ENC = sum(enc) # total biomass encountered
	con = cmax .* enc[1] ./ (cmax + ENC) # Type II
  #con = cmax
  return con
end


###! Offline coupling
function sub_offline_zm(enc_1,enc_2,enc_3,enc_4,enc_5,bio_1,bio_2,bio_3,bio_4,bio_5,dZ)
  # ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
  #! offline switch
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
  # ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
  #! offline switch
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
	#! calculates consumption rate of first element of enc
  #Cmax
  temp = (Tp.*tdif) + (Tb.*(1.0-tdif))
  #! Specific ingestion rate from Hartvig et al (g/g/day) ref to 15C
  cmax = (exp(0.063*(temp-15.0)) * 60.0 * wgt^(-0.25)) ./365.0
  #! Hartvig et al (g/g/day) ref to 10C
  #cmax = (exp(0.063*(temp-10.0)) * 85.0 * wgt^(-0.25)) ./365.0
  #! mizer (g/g/day) ref to 10C
  #cmax = (exp(0.063*(temp-10.0)) * 40.0 * wgt^(-0.25)) ./365.0
  #! J&C15 (g/g/day) ref to 10C
  #cmax = (exp(0.063*(temp-10.0)) * 25.0 * wgt^(-0.33)) ./365.0
  #! Specific ingestion rate from Kiorboe & Hirst (g/g/day)
  #cmax = (exp(0.063*(temp-15.0)) * 10^(0.4) * wgt^(-0.51)) .* 24e-3
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
	#nu = ((I/B)*Lambda) - met
  #nu = 0.5*I
  # Already in biomass specific ingestion
	nu = (I*Lambda) - met
  prod = nu * B
  return nu, prod
end

function sub_nu_MD(I,B,met)
	# convert to biomass specific ingestion
	#nu = ((I/B)*Lambda) - met
  #nu = 0.5*I
  # Already in biomass specific ingestion
	nu = (I*0.6) - met
  prod = nu * B
  return nu, prod
end

function sub_nu_LD(B,met,c1,c2,c3,c4)
	# convert to biomass specific ingestion
	#nu = ((I/B)*Lambda) - met
  #nu = 0.5*I
  # Already in biomass specific ingestion
  Ip = c1+c2+c3
  Ib = c4
	nu = (Ip*Lambda + Ib*0.6) - met
  prod = nu * B
  return nu, prod
end


###! ENERGY AVAILABLE FOR SOMATIC GROWTH
function sub_gamma(K,Z,nu,d,B,S,nmrt,dfrate,selec)
  # d = predation loss
  # nmort = natural mortality rate
  # dfrate = fishing mortality rate
  # selec = if harvested
  # convert predation mortality to biomass specific rate
  if (selec == 1)
	   D = (d/B) + nmrt + dfrate
  else
     D = (d/B) + nmrt
  end
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
            rho = (1.0-K) .* nu;  #energy available for from eating
        else
            rho = 0.0;
        end
        if S>0.0
            #rep = eggs from energy now + eggs from stored energy
            rep = S.*(rho+egg);         #fraction of pop reproducing now
            egg = (1.0-S).*(rho+egg);   #rest gets stored for later
        else
            rep = 0.0;
            egg = egg + rho;
        end
    else
        rep = 0.0;
        egg = 0.0;
    end
    if nu > 0.0
        #nu is now split into used for repro (nu) and stored (egg)
        nu = rep;
    end
	return nu, rep, egg
end


###! Biomass recruiting to larvae (g m-2 d-1)
function sub_rec_larv(X,bio,rfrac,Tp,Tb,tpel)
	#X: repro rate
  #rfrac: repro efficiency ~ sex ratio * egg survival
  #Tp: pelagic temp
  #Tb: bottom temp
  #tpel: frac pelagic time
  #! Global constant RE
  rec = rfrac * X * bio

  Tavg = (Tp.*tpel) + (Tb.*(1.0-tpel))
  #! Exponential for Demersal A = 100% or 50%
  #RE = 10^( -6.275*exp(-0.03527*Tavg) + 3.837*exp(-0.3701*Tavg) )
  #! A 3rd order polynomial for Demersal A = 100% or 50%
  #RE = 10^( -1.027e-3*Tavg^3 + 0.04538*Tavg^2 - 0.4882*Tavg - 3.016 )
  #! B 3rd order polynomial for Demersal A = 100% or 50%
  #RE = 10^( -7.864e-4*Tavg^3 + 0.03782*Tavg^2 - 0.4755*Tavg - 2.276 )
  #! A piecewise for Demersal A = 100% or 50%
  # if (Tavg < 11)
  #   RE = 10^( -3.856 + -0.05985 .* Tavg );
  # else
  #   RE = 10^( -5.607 + 0.1227 .* Tavg );
  # end
  #! B piecewise for Demersal A = 100% or 50%
  # if (Tavg < 11)
  #   RE = 10^( -2.762 + -0.152 .* Tavg );
  # else
  #   RE = 10^( -5.608 + 0.1227 .* Tavg );
  # end

  #! Exponential for Blanchard A & assim
  #RE = 10^( -4.742*exp(-0.02429*Tavg) + 2.788*exp(-0.4168*Tavg) )
  #! A 3rd order polynomial for Blanchard A & assim
  #RE = 10^( -9.865e-4*Tavg^3 + 0.04195*Tavg^2 - 0.4479*Tavg - 2.549 )
  #! B 3rd order polynomial for Blanchard A & assim
  #RE = 10^( -6.893e-4*Tavg^3 + 0.03076*Tavg^2 - 0.36*Tavg - 2.307 )
  #! A piecewise for Blanchard A & assim
  # if (Tavg < 11)
  #   RE = 10^( -3.34 + -0.04833 .* Tavg );
  # else
  #   RE = 10^( -4.591 + 0.08163 .* Tavg );
  # end
  #! B piecewise for Blanchard A & assim
  # if (Tavg < 11)
  #   RE = 10^( -2.715 + -0.09197 .* Tavg );
  # else
  #   RE = 10^( -4.591 + 0.08163 .* Tavg );
  # end
  #
  # rec = RE * X * bio

  return rec
end

###! Biomass recruiting to next size-class (g m-2 d-1)
function sub_rec(X,bio)
	#X: biomass specific maturation rate
  rec = X * bio
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
    #! Const with size
    #nmort = exp(0.063*(temp-15.0)) * Nat_mrt
    #! Hartvig
    #nmort = exp(0.063*(temp-10.0)) * 0.84 * wgt^(-0.25) /365.0;
    #! mizer
    nmort = exp(0.063*(temp-10.0)) * 3.0 * wgt^(-0.25) /365.0;
    #! J&C
    #nmort = exp(0.063*(temp-10.0)) * 0.5 * wgt^(-0.33) /365.0;
    #! Intermediate
    #nmort = exp(0.063*(temp-10.0)) * 1.5 * wgt^(-0.33) /365.0;
    #! Peterson & Wroblewski
    #nmort = exp(0.063*(temp-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25); #daily & uses dry weight
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
  #db = rec + ((nu - egg - rep - gamma - nmort) * bio_in) + (egg * bio_in) - die
  db = rec + ((nu + egg - rep - gamma - nmort) .* bio_in) - die
  bio_out =  bio_in + db
end

###! Update biomass
function sub_update_lg(bio_in,rec,nu,rep,gamma,die,egg,nmort)
	# all inputs except rec are in g g-1 d-1; rec is g d-1
	# rec = rep from smaller size class = TOTAL biomass gained from recruitment
	# grw = nu = somatic energy for growth within size class
  # store = egg = energy stored for later egg production
	# rep = rep =  energy lost to egg production
	# mat = gamma = energy lost to maturation to larger size class
	# nmort = natural mortality
	# die = predator mort = biomass lost to predation
  ### Higher predation on L by 2m predator, assume biomass = 0.1 * L bio
  # L = 2000.0;
  # M = 0.01 * (0.1*L)^3;
  # A = exp(0.063*(temp-15.0)) * 1.74e3 * M^(-0.24) * (24e-3/9)
  # enc = A * 0.1 * bio_in^2
  # cmax = exp(0.063*(temp-15.0)) * 2.5 .* M^(-0.51) * 24e-3
  # con = cmax * enc / (cmax + enc)
  # nmort = con * 0.1 * bio_in  ## Do you multiply by bio again?
  ### Simplified
  # cmax = exp(0.063*(temp-15.0)) * 2.5 .* M^(-0.51) * 24e-3
  # Khp = 6.1e-2;
  # nmort = cmax * bio_in / (Khp + bio_in)
  db = rec + ((nu - egg - rep - gamma - nmort) * bio_in) + (egg * bio_in) - die
  bio_out =  bio_in + db
end

function sub_update_be(bio_in,BE,det,CC,con,bio)
  #bio_in = benthic invert biomass
  #Tb = bottom temp
  #con = vector of biomass specific consumption by MD and LD
  #bio = vector MD and LD biomass
  #omort = other mortality

  #! Peterson & Wroblewski: daily & uses dry weight
  #nmort = exp(0.063*(temp-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25);
  # Use Wei et al mean (1.3956 gWW) or median (0.7656) wgt
  # wgt = 0.7656;
  # omort = exp(0.063*(Tb-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25)
  # eaten = con.*bio
  # bio_out = bio_in - sum(eaten) - bio_in*omort

  #! Quadratic mortality from carrying capacity
  r = BE*det
  eaten = con.*bio
  # Chemostat
  bio_out = bio_in + r * (r*CC - bio_in) - sum(eaten)
  # Logistic
  #bio_out = bio_in + r * bio_in * (1 - bio_in/CC) - sum(eaten)

end


####! Fishing
function sub_fishing_mass(MFbio,LPbio,LDbio,AREA)
	if FISHING > 0.0
    ALL_pl = MFbio .* AREA
		ALL_pi = LPbio .* AREA
		ALL_de = LDbio .* AREA

		#! Total fish biomass
		TOT = sum(ALL_pi) + sum(ALL_pl) + sum(ALL_de)
    ALL_pl -= (ALL_pl./TOT) .* FISHING
		ALL_pi -= (ALL_pi./TOT) .* FISHING
		ALL_de -= (ALL_de./TOT) .* FISHING

		#! Calc total biomass of fish in the ocean
    MFbio = ALL_pl ./ AREA
		LPbio = ALL_pi ./ AREA
		LDbio = ALL_de ./ AREA
	end
	return MFbio, LPbio, LDbio
end

function sub_fishing_rate(bio,FISHING,selec)
  #selec = if being harvested
  if (selec==1)
  	caught = bio * FISHING
    bio -= caught
  else
    bio = bio
    caught = 0.0
  end
	return bio, caught
end


###! Forward Euler checks
function sub_check!(bio)
	ID = find(bio .< 0)
	bio[ID] = eps()
  #ID2 = find(isnan(bio))
  #bio[ID2] = eps()
end
