

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
plt.rc('legend', fontsize=15)    # legend fontsize


qois = ['x1x3', 'x4x5', 'x0x4', 'x3x5']
prefix = "./run_tests/s3d_logs/"

n_qoi = 4
n_eb = 20

pmgard_br = np.zeros([n_qoi,n_eb])
sz3_br = np.zeros([n_qoi,n_eb])
sz3delta_br = np.zeros([n_qoi,n_eb])

a = 0.1
r = 0.5
n = 20
end = a * (r ** (n - 1))
request_e = np.geomspace(a, end, num=n)

vr = np.array([2.32979e-04, 2.17227e-04, 5.85201e-05, 3.10047e-05])
 
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

styles=['r-o', 'g:^', 'b--*', 'c-o', 'm:^', 'y--*']
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,8))
ax=[axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]

labels = ['PMGARD-HB', 'PSZ3', 'PSZ3-delta']
 
i = 0
while i < 4:
                names = [r'$x_1x_3$', r'$x_4x_5$', r'$x_0x_4$', r'$x_3x_5$']
                p, = ax[i].plot(pmgard_br[:, i][4:], request_e[4:], styles[0], label='{}'.format(labels[0]))
                p, = ax[i].plot(sz3_br[:, i][4:], request_e[4:], styles[1], label='{}'.format(labels[1]))
                p, = ax[i].plot(sz3delta_br[:, i][4:], request_e[4:], styles[2], label='{}'.format(labels[2]))
                ax[i].set_ylabel('Relative QoI Error',fontsize=22)
                ax[i].set_yscale('log')
                # ax_t.set_ylim(0, 2000)
                ax[i].set_xlim(0, 15)
                if i == 1 or i == 3:
                    ax[i].set_xlim(0, 10)
                # ax_t.set_xlim(0, 30)
                # ax[i].grid(which='major', axis='y')
                ax[i].set_yticks([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7])
                ax[i].set_title(names[i],fontsize=22)
                ax[i].set_xlabel('Bitrate',fontsize=22)
                ax[i].legend(loc='upper right')
                i += 1

plt.tight_layout()
plt.savefig('./plots/figures/fig8.png')
plt.show()
