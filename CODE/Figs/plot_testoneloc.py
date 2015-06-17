import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt


##! Data
PISC = genfromtxt('../Data/CSV/Spinup_PISC.csv', delimiter=',')
PLAN = genfromtxt('../Data/CSV/Spinup_PLAN.csv', delimiter=',')
DETR = genfromtxt('../Data/CSV/Spinup_DETR.csv', delimiter=',')
BENT = genfromtxt('../Data/CSV/Spinup_BENT.csv', delimiter=',')

##! Plot
#x = np.arange(1,PISC.shape[0]+1,1);
lookback = 10*365 # number of years to look back at
#lookback = PISC.shape[0]-1
x = (np.arange(0,lookback))/365. # plot last two years

plt.figure(figsize=(10,5))
ax = plt.axes()
ax.set_color_cycle([plt.cm.copper_r(i) for i in np.linspace(0, 1, PISC.shape[1])])
for i in np.arange(0,PISC.shape[1]-1):
    plt.plot(x,(PISC[-lookback-1:-1,i]))
#plt.axis([10200,12000,0.,0.00015])
plt.legend(['small','','','','','','','','big'],bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.)
pos1 = ax.get_position()
pos2 = [pos1.x0, pos1.y0,  pos1.width / 1.2, pos1.height]
ax.set_position(pos2) 
plt.xlabel('Time (years)',fontsize=14)
plt.ylabel('Piscivore biomass (g m-2)',fontsize=14)
plt.savefig('./PDF/Fig_POEM_pisc.pdf',dpi=200)


plt.figure(figsize=(10,5))
ax = plt.axes()
ax.set_color_cycle([plt.cm.cool(i) for i in np.linspace(0, 1, PLAN.shape[1])])
for i in np.arange(0,PLAN.shape[1]-1):
    plt.plot(x,(PLAN[-lookback-1:-1,i]))
#plt.axis([10200,12000,0.,0.00015])
plt.legend(['small','','','','','','','','','','big'],bbox_to_anchor=(1.05, 1),loc=2,                borderaxespad=0.)
pos1 = ax.get_position()
pos2 = [pos1.x0, pos1.y0,  pos1.width / 1.2, pos1.height]
ax.set_position(pos2)
plt.xlabel('Time (years)',fontsize=14)
plt.ylabel('Planktivore biomass (g m-2)',fontsize=14)
plt.savefig('./PDF/Fig_POEM_plan.pdf',dpi=200)


plt.figure(figsize=(10,5))
ax = plt.axes()
ax.set_color_cycle([plt.cm.winter_r(i) for i in np.linspace(0, 1, DETR.shape[1])])
for i in np.arange(0,DETR.shape[1]-1):
    plt.plot(x,(DETR[-lookback-1:-1,i]))
plt.legend(['small','','','','big'],bbox_to_anchor=(1.05, 1),loc=2,                            borderaxespad=0.)
pos1 = ax.get_position()
pos2 = [pos1.x0, pos1.y0,  pos1.width / 1.2, pos1.height]
ax.set_position(pos2)
plt.xlabel('Time (years)')
plt.ylabel('Detritivore biomass (g m-2)')
plt.savefig('./PDF/Fig_POEM_detr.pdf',dpi=200)


plt.figure(); ax = plt.axes()
plt.plot(x,BENT[-lookback-1:-1])
plt.xlabel('Time (years)')
plt.ylabel('Detrital pool biomass (g m-2)')


plt.show()

