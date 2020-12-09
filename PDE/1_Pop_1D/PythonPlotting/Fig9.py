import numpy as np
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
del matplotlib.font_manager.weight_dict['roman']
matplotlib.font_manager._rebuild()
from math import floor, log10
from matplotlib import rc
plt.rcParams["font.family"] = "Times New Roman"
rc('font',**{'family':'serif','serif':['Time New Roman']})
rc('text', usetex=True)

#function for rounding to significant figure
def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)

#load data
lowk = np.array(np.load("./data/1D_sync_lowk.npy")) #.87
midk = np.array(np.load("./data/1D_sync_midk.npy"))#.88
highk = np.array(np.load("./data/1D_sync_highk.npy")) #.92
highhighk = np.array(np.load("./data/1D_sync_highhighk.npy")) #1.05

#plotting parameters
sizer1 = np.size(lowk[:, 1])
sizer2 = np.size(lowk[1,:])
data = np.zeros((sizer1,sizer2,4))
data[:,:,0] = lowk
data[:,:,1] = midk
data[:,:,2] = highk
data[:,:,3] = highk
T = 250
Y2 = (np.linspace(-5*np.pi, 5*np.pi, sizer1))
Y1 = (np.linspace(0, T , sizer2))
Y1, Y2 = np.meshgrid(Y1, Y2)
levels = np.linspace(np.amin(data), np.amax(data), 200)
fig, ax = plt.subplots(1,3, sharex='col', sharey='row',figsize=(20,6))

#plotting
for j in range(3):
    im = ax[j].contourf(Y1,Y2,data[:,:,j], levels = levels, cmap = "inferno")
    ax[j].tick_params(which='both', labelsize=15)
    for c in im.collections:
        c.set_edgecolor("face")
    if j == 0:
        ax[j].set_ylabel("x", fontsize=20, rotation=0, labelpad=1)
        ax[j].yaxis.set_ticks(np.linspace(-5 * np.pi, 5 * np.pi, 2))
        ax[j].set_yticklabels(['$-5 \pi$', '$5 \pi$'])

    ax[j].set_xlabel("time (ms)", fontsize=20)

ax[0].text(125, 16.5, 'I',va='center',ha='center',fontsize = 18)
ax[1].text(125, 16.5, 'II', va='center',ha='center',fontsize = 18)
ax[2].text(125, 16.5,'III', va='center',ha='center',fontsize = 18)


cbar = fig.colorbar(im, ax=ax.ravel().tolist(),ticks = np.linspace(np.amin(data),np.amax(data),2))
cbar.ax.set_ylabel('$|$Z$|$',fontsize = 20,rotation = 0)
cbar.ax.tick_params(labelsize=15)
DirSave = '/home/james/PhD_Work/Figures/BrainTopPaper/'
plt.savefig(DirSave + 'Fig9.pdf')
plt.show()
