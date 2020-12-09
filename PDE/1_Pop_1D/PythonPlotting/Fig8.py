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
#rc('font',**{'family':'sans-serif','sans-serif':['Times New Roman']})
rc('font',**{'family':'serif','serif':['Time New Roman']})
rc('text', usetex=True)

def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)
turing=open('../bifurcation_Diagrams/hopf_kappaV_v_eta0_0.3_kappaS_5_delta_0.5_alpha_3_tau_1.dat',"r")
lines1=turing.readlines()

turinghopf1=open('../bifurcation_Diagrams/turing_hopf_kappaV_v_eta0_0.3_kappaS_5_delta_0.5_alpha_3_tau_1.dat',"r")
lines2=turinghopf1.readlines()

hopf1 = np.array(np.load("../bifurcation_Diagrams/1D_hopf_curve.npy"))[::-1]
turinghopf2 = np.array(np.load("../bifurcation_Diagrams/1D_turinghopf_curve.npy"))[::-1]




def extractdata(lines,col,st,en):
    data = []
    for x in lines:
        data.append(x.split(' ')[col])

    for i in range(np.size(data)):
        data[i] = float(data[i])

    return data[st:en]

st1 = 0
en1 = 50

st2= 0
en2 = 0

st22= 165
en22 = -1850


v1 = extractdata(lines1,0,st1,en1)

v2 = extractdata(lines2,0,st2,en2)
v22 = extractdata(lines2,0,st22,en22)


curve1 = extractdata(lines1,1,st1,en1)
curve2 = extractdata(lines2,1,st2,en2)
curve22 = extractdata(lines2,1,st22,en22)






#print(np.arange(0,501,2))

#g = np.array(np.load("u_1D.npy"))
en = 150
hopf = np.array(np.load("data/1D_sync_hopf.npy"))[:,0:-1] #v = 0.1 kv = 0.85
thopf1 = np.array(np.load("data/1D_sync_turinghopf1.npy"))[:,0:-1] #v = 0.1  kv = 0.87
thopf2 = np.array(np.load("data/1D_sync_turinghopf2.npy"))[:,0:-1] #v = 0.8  kv = 0.875




print(np.size(hopf,1))
print(np.size(hopf))
sizer1 = np.size(hopf[:, 1])
sizer2 = np.size(hopf[1,:])
data = np.zeros((sizer1,sizer2,3))
data[:,:,0] = hopf
data[:,:,1] = thopf1
data[:,:,2] = thopf2



T = sizer2*0.1
Y2 = (np.linspace(-5*np.pi, 5*np.pi, sizer1))
Y1 = (np.linspace(0, 250 , sizer2))
print(T/2)


Y1, Y2 = np.meshgrid(Y1, Y2)
levels = np.linspace(np.amin(data), np.amax(data), 200)

fig, ax = plt.subplots(1,4,figsize=(20,6))


for j in range(1,4):
    im = ax[j].contourf(Y1,Y2,data[:,:,j-1], levels = levels, cmap = "inferno")
    ax[j].tick_params(which = 'both',labelsize = 15)
    for c in im.collections:
        c.set_edgecolor("face")
    if j == 1:
        ax[j].set_ylabel("x",fontsize = 20,rotation = 0,labelpad = -10)
        ax[j].yaxis.set_ticks(np.linspace(round_sig(-5*np.pi,3), round_sig(5*np.pi,3), 2))
        ax[j].set_yticklabels(['$-5\pi$','$5\pi$'])
    if j > 1:
        ax[j].axes.yaxis.set_visible(False)

    ax[j].set_xlabel("time (ms)", fontsize=20)

ax[1].text(125, 16.5, 'I', va='center',ha='center',fontsize = 20)
ax[2].text(125, 16.5, 'II', va='center',ha='center',fontsize = 20)
ax[3].text(125, 16.5, 'III', va='center',ha='center',fontsize = 20)

v_space = np.linspace(0.02,1.2,1000)

ax[0].plot(v_space,hopf1,'--',color='blue',label = 'Hopf')
ax[0].plot(v_space,turinghopf2,color='red',label = 'Turing-Hopf')
ax[0].tick_params(which = 'both',labelsize = 15)
ax[0].set_ylabel('$\kappa_v$',fontsize = 20,rotation=0,labelpad = 10)
ax[0].set_xlabel('v',fontsize = 20,rotation=0)
ax[0].text(0.17,0.87,'II',fontsize=20)
ax[0].text(.16,.845,'I',fontsize=20)
ax[0].text(0.9,0.87,'III',fontsize=20)
ax[0].legend(loc = 'upper left',fontsize=15)
ax[0].margins(x=0)



im = ax[2].contourf(Y1,Y2,data[:,:,1], levels = levels, cmap = "inferno")
cbar = fig.colorbar(im, ax=ax.ravel().tolist(),ticks = np.linspace(np.amin(data[:,:,:]),round_sig(np.amax(data[:,:,:]-0.001),3),2))
cbar.ax.set_ylabel('$|$Z$|$',fontsize = 20,rotation = 0,labelpad=-20)
cbar.ax.tick_params(labelsize=15)
DirSave = '/home/james/PhD_Work/Figures/BrainTopPaper/'
plt.savefig(DirSave + 'Fig8.pdf')
plt.show()
