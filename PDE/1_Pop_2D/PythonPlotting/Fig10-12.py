import numpy as np
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from math import log10, floor
del matplotlib.font_manager.weight_dict['roman']
matplotlib.font_manager._rebuild()
from matplotlib import rc
plt.rcParams["font.family"] = "Times New Roman"
rc('font',**{'family':'serif','serif':['Time New Roman']})
rc('text', usetex=True)

def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)

Z = np.array(np.load("./data/2D_TuringHopf_kappaV=0.695_NEW.npy")) #beta = 7

print(np.shape(Z))

t = np.linspace(0,500,np.size(Z[1,1,:]))
sync = np.abs(Z)
W = (1-np.conj(Z))/(1+np.conj(Z))
R = np.real(W)/np.pi
V = np.imag(W)
phase = np.angle(Z)

s_pos1 = 100
s_pos2 = 120
cpos1 = s_pos1*((np.pi)/(200))-np.pi
cpos2 = s_pos2*((np.pi)/(200))-np.pi




sizer1 = np.size(sync[:, 1, 1])
sizer2 = np.size(sync[1, :, 1])
data = np.zeros((sizer1,sizer2,4))
types = np.array([[80,100], [120,140]])
time  = 14
data[:,:,0] = R[:,:,time]
data[:,:,1] = V[:,:,time]
data[:,:,2] = sync[:,:,time]
data[:,:,3] = phase[:,:,time]
T = sizer2*0.1
Y2 = (np.linspace(-np.pi, np.pi, sizer1))
Y1 = (np.linspace(-np.pi, np.pi, sizer2))
Y1, Y2 = np.meshgrid(Y1, Y2)
nlevels = 100
levels = np.zeros((nlevels,4))
cbarticks = np.zeros((2,4))


def r1(x):
    return round(x, -int(floor(log10(abs(x)))))

fig, ax = plt.subplots(2,4,sharey="col",figsize=(18,10))


for i in range(4):
        print(np.amin(data[:,:,i]))
        levels[:,i] = np.linspace(np.amin(data[:,:,i]), np.amax(data[:,:,i]), nlevels)
        cbarticks[:,i] = np.linspace(np.amin(data[:,:,i]), np.amax(data[:,:,i]), 2)

labels = ["R","V","$|$Z$|$",r"$\theta$"]

for i in range(4):
        im = ax[0,i].contourf(Y1,Y2,data[:,:,i], levels = levels[:,i], cmap = "inferno")
        #ax[i].text(6.0,16,] ,fontsize = 20)
        cbar = plt.colorbar(im, ax=ax[0, i],orientation="horizontal",ticks = cbarticks[:,i])
        cbar.set_label(labels[i], fontsize=20)
        cbar.ax.tick_params(labelsize=20)
        ax[0,i].tick_params(axis='both', which='both', labelsize=20)
        ax[0,i].yaxis.set_ticks(np.linspace(round_sig(-np.pi, 3), round_sig( np.pi, 3), 2))
        ax[0,i].set_yticklabels(['$-\pi$', '$\pi$'])
        ax[0, i].xaxis.set_ticks(np.linspace(round_sig(-np.pi, 3), round_sig(np.pi, 3), 2))
        ax[0, i].set_xticklabels(['$-\pi$', '$\pi$'])
        for c in im.collections:
            c.set_edgecolor("face")
        if i == 0:
            ax[0,i].set_ylabel("y", fontsize=25,rotation = 0 )
        else:
            ax[0,i].set_yticks([], [])

        ax[0,i].set_xlabel("x", fontsize=25,labelpad=-5)
        #ax[0,i].xaxis.set_ticks_position("top")


gs = ax[1, 0].get_gridspec()
for j in range(4):
    fig.delaxes(ax[1,j])
axbig = fig.add_subplot(gs[1:, 0:2])
axbig2 = fig.add_subplot(gs[1:, 2:4])

lpad = -15
stime = -0
axbig.plot(t,R[s_pos1,s_pos2,:],label="R",color='teal')
axbig.tick_params(axis='both', which='both', labelsize=20)
axbig.set_yticks(np.linspace(round_sig(np.amin(R),2),round_sig(np.amax(R),2),2))
axbig.set_ylabel('R', color='teal',rotation = 0, fontsize = 25,labelpad = lpad)
axbigshare = axbig.twinx()
axbigshare.plot(t,V[s_pos1,s_pos2,:],'--',label="V",color='orange')
axbigshare.set_ylabel('V', color='orange',rotation = 0, fontsize = 25,labelpad = lpad)
axbigshare.set_yticks(np.linspace(round_sig(np.amin(V),3),round_sig(np.amax(V),3),2))
axbigshare.tick_params(axis='both', which='both', labelsize=20)
axbig.legend(fontsize=16,loc="lower left",bbox_to_anchor=(0.6, -0.2))
axbigshare.legend(fontsize=16,loc="lower left",bbox_to_anchor=(0.8, -0.2))


axbig2.plot(t,phase[s_pos1,s_pos2,:],'--',label=r"$\theta$",color='red')
axbig2.tick_params(axis='both', which='both', labelsize=20)
axbig2.set_ylabel(r'$\theta$', color='red',rotation = 0, fontsize = 25,labelpad = lpad)
axbig2.set_yticks(np.linspace(round_sig(-np.pi,3),round_sig(np.pi,3),2))
axbig2.set_yticklabels(['$-\pi$','$\pi$'],fontsize = 25)
axbig2share = axbig2.twinx()
axbig2share.plot(t,sync[s_pos1,s_pos2,:],label="$|$Z$|$",color='green')
axbig2share.set_ylabel('$|$Z$|$', color='green',rotation = 0, fontsize = 25,labelpad = lpad+15)
axbig2share.set_yticks(np.linspace(np.amin(sync),np.amax(sync),2))
axbig2share.tick_params(axis='both', which='both', labelsize=20)
axbig2.legend(fontsize=16,loc="lower left",bbox_to_anchor=(0.6, -0.2))
axbig2share.legend(fontsize=16,loc="lower left",bbox_to_anchor=(0.8, -0.2))


#axbig.legend(fontsize = 20,loc = 'upper right')

#axbig2.legend(fontsize = 20,loc = 'upper right')
#fig.subplots_adjust(right=0.8)

axbig.set_xlabel("time (ms)",fontsize = 20)
axbig2.set_xlabel("time (ms)",fontsize = 20)

axbig.margins(x=0)
axbig2.margins(x=0)


circle0 = plt.Circle((cpos1,cpos2), radius= 0.1,color='lime',fill=True,linewidth=1)
circle1 = plt.Circle((cpos1,cpos2), radius= 0.1,color='lime',fill=True,linewidth=1)
circle2 = plt.Circle((cpos1,cpos2), radius= 0.1,color='lime',fill=True,linewidth=1)
circle3 = plt.Circle((cpos1,cpos2), radius= 0.1,color='lime',fill=True,linewidth=1)
ax[0,0].add_artist(circle0)
ax[0,1].add_artist(circle1)
ax[0,2].add_artist(circle2)
ax[0,3].add_artist(circle3)

#cbar_ax = fig.add_axes([0.82, 0.11, 0.02, 0.77])
#fig.colorbar(im, cax=cbar_ax,ticks=cbarticks[:,i])
#cbar = plt.colorbar(im,pad=2, ticks=cbarticks[:,i])
#cbar.ax.set_ylabel('|Z|',labelpad=20, fontsize=20, rotation=0)
plt.tight_layout()
DirSave = '/home/james/PhD_Work/Figures/BrainTopPaper/'
plt.savefig(DirSave + 'Fig10.pdf',dots = 2000)
plt.show()
