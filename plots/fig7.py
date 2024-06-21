

import matplotlib.pyplot as plt
import numpy as np
import re
 
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14


# plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=16)  # fontsize of the figure title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('legend', fontsize=17)    # legend fontsize


species = ['PMGARD-HB', 'PSZ3', 'PSZ3-delta']

qois = ["VTOT", "T", "C", "Mach", "PT", "mu"]
prefix = "./run_tests/ge_logs/"

n_eb = 20
n_qoi = len(qois)

pmgard_br = np.zeros([n_qoi,n_eb])
sz3_br = np.zeros([n_qoi,n_eb])
sz3delta_br = np.zeros([n_qoi,n_eb])
vr = np.array([626.443, 236.349, 127.804, 1.93198, 196480, 1.04299e-05])


for i in range(len(qois)):
    qoi = qois[i]
    pmgard_path = prefix + "pmgard_" + qoi + ".log"
    with open(pmgard_path, 'r') as file:
        content = file.read()
    pmgard_cr = [float(match) for match in re.findall(r'aggregated cr = (-?[\d.]+)', content)]
    pmgard_br[i] = 64 / np.array(pmgard_cr)
    
    sz3_path = prefix + "psz3_" + qoi + ".log"
    with open(sz3_path, 'r') as file:
        content = file.read()
    sz3_cr = [float(match) for match in re.findall(r'aggregated cr = (-?[\d.]+)', content)]
    sz3_br[i] = 64 / np.array(sz3_cr)
    
    sz3delta_path = prefix + "psz3delta_" + qoi + ".log"
    with open(sz3delta_path, 'r') as file:
        content = file.read()
    sz3delta_cr = [float(match) for match in re.findall(r'aggregated cr = (-?[\d.]+)', content)]
    sz3delta_br[i] = 64 / np.array(sz3delta_cr)

 
pmgard_br = np.transpose(pmgard_br)
sz3_br = np.transpose(sz3_br)
sz3delta_br = np.transpose(sz3delta_br)

request_e = np.geomspace(0.1, 0.1*(0.5**19), num=20)

styles=['r-o', 'g:^', 'b--*', 'c-o', 'm:^', 'y--*']
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(16,15))
ax=[axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1], axs[2, 0], axs[2, 1]]

i = 0
while i < 6:
                names = ["VTOT", "T", "C", "Mach", "PT", "mu"]
                p, = ax[i].plot(pmgard_br[:, i], request_e, styles[0], label='{}'.format(species[0]))
                p, = ax[i].plot(sz3_br[:, i], request_e, styles[1], label='{}'.format(species[1]))
                p, = ax[i].plot(sz3delta_br[:, i], request_e, styles[2], label='{}'.format(species[2]))

                ax[i].set_ylabel('Relative QoI Error',fontsize=30)
                ax[i].set_yscale('log')
                ax[i].set_yticks([1e-1, 1e-2, 1e-3, 1e-4, 1e-5,1e-6])
                # ax_t.set_ylim(0, 2000)
                ax[i].set_xlim(4, 25)
                # ax_t.set_xlim(0, 30)
                # ax[i].grid(which='major', axis='y')
                ax[i].set_title(names[i],fontsize=30)
                ax[i].set_xlabel('Bitrate',fontsize=30)
                ax[i].legend(loc='lower left')
                i += 1

plt.tight_layout()
plt.savefig('./plots/figures/fig7.png')
plt.show()
