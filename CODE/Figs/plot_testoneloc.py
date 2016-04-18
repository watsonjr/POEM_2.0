import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
import seaborn

##! Data
Sml_f = genfromtxt('../Data/CSV/Spinup_Sml_f.csv', delimiter=',')
Sml_p = genfromtxt('../Data/CSV/Spinup_Sml_p.csv', delimiter=',')
Sml_d = genfromtxt('../Data/CSV/Spinup_Sml_d.csv', delimiter=',')
Med_f = genfromtxt('../Data/CSV/Spinup_Med_f.csv', delimiter=',')
Med_p = genfromtxt('../Data/CSV/Spinup_Med_p.csv', delimiter=',')
Med_d = genfromtxt('../Data/CSV/Spinup_Med_d.csv', delimiter=',')
Lrg_p = genfromtxt('../Data/CSV/Spinup_Lrg_p.csv', delimiter=',')

##! Plot
plt.close('all')

#x = np.arange(1,PISC.shape[0]+1,1);
#lookback = 10*365 # number of years to look back at
lookback = Sml_f.shape[0]-1
x = (np.arange(0,lookback))/365. # plot last two years

plt.figure(figsize=(10,5))
ax = plt.axes()
plt.plot(x,(Sml_f[-lookback-1:-1]))
plt.plot(x,(Sml_p[-lookback-1:-1]))
plt.plot(x,(Sml_d[-lookback-1:-1]))
plt.xlabel('Time (years)',fontsize=14)
plt.ylabel('Biomass (g m-2)',fontsize=14)
plt.savefig('./PDF/Fig_POEM_plank.pdf',dpi=200)

plt.figure(figsize=(10,5))
ax = plt.axes()
plt.plot(x,(Med_f[-lookback-1:-1]))
plt.plot(x,(Med_p[-lookback-1:-1]))
plt.plot(x,(Med_d[-lookback-1:-1]))
plt.xlabel('Time (years)',fontsize=14)
plt.ylabel('Biomass (g m-2)',fontsize=14)
plt.savefig('./PDF/Fig_POEM_pisc.pdf',dpi=200)

plt.figure(figsize=(10,5))
ax = plt.axes()
plt.plot(x,(Lrg_p[-lookback-1:-1]))
plt.xlabel('Time (years)',fontsize=14)
plt.ylabel('Biomass (g m-2)',fontsize=14)
plt.savefig('./PDF/Fig_POEM_detr.pdf',dpi=200)

plt.show()

