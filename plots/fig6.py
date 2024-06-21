

import matplotlib.pyplot as plt
import numpy as np
import re
from os import system

 
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
plt.rc('legend', fontsize=16)    # legend fontsize

qois = ['x1x3', 'x4x5', 'x0x4', 'x3x5']
prefix = "./run_tests/s3d_logs/"

a = 0.1
r = 0.5
n = 20
end = a * (r ** (n - 1))
request_e = np.geomspace(a, end, num=n)

n_qoi = 4
n_eb = 20

br = np.zeros([n_qoi,n_eb])
estimate_e = np.zeros([n_qoi,n_eb])
max_e = np.zeros([n_qoi,n_eb])

for i in range(len(qois)):
    qoi = qois[i]
    pmgard_path = prefix + "pmgard_" + qoi + ".log"
    with open(pmgard_path, 'r') as file:
        content = file.read()
    pmgard_error_est = [float(match) for match in re.findall(r'max_est_error = (-?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content)]
    pmgard_error_act = [float(match) for match in re.findall(r'max_act_error = (-?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content)]
    pmgard_cr = [float(match) for match in re.findall(r'aggregated cr = (-?[\d.]+)', content)]
    br[i] = 64 / np.array(pmgard_cr)
    estimate_e[i] = np.array(pmgard_error_est)
    max_e[i] = np.array(pmgard_error_act)
    
br = np.transpose(br)
estimate_e = np.transpose(estimate_e)
max_e = np.transpose(max_e)
vr = np.array([2.32979e-04, 2.17227e-04, 5.85201e-05, 3.10047e-05])
    

styles=['r-o', 'g:^', 'b--*', 'c-o', 'm:^', 'y--*']
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,8))
ax=[axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]

labels = ['Requested tolerance', 'Max estimated error', 'Max real error']
 
i = 0
index = [0, 1, 2, 3]
while i < 4:
                names = [r'$x_1x_3$', r'$x_4x_5$', r'$x_0x_4$', r'$x_3x_5$']
                p, = ax[i].plot(br[:, index[i]][4:], request_e[4:], styles[0], label='{}'.format(labels[0]))
                p, = ax[i].plot(br[:, index[i]][4:], (estimate_e[:, index[i]]/vr[index[i]])[4:], styles[1], label='{}'.format(labels[1]))
                p, = ax[i].plot(br[:, index[i]][4:], (max_e[:, index[i]]/vr[index[i]])[4:], styles[2], label='{}'.format(labels[2]))
                ax[i].set_ylabel('Relative QoI Error', fontsize=24)
                ax[i].set_yscale('log')
                # ax_t.set_ylim(0, 2000)
                ax[i].set_xlim(0, 10)
                if i==1 or i ==3:
                    ax[i].set_xlim(0, 8)
                # ax_t.set_xlim(0, 30)
                # ax[i].grid(which='major', axis='y')
                ax[i].set_yticks([1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7])
                ax[i].set_title(names[index[i]],fontsize=24)
                ax[i].set_xlabel('Bitrate',fontsize=20)
                ax[i].legend(loc='upper right')
                if i==0 or i==2:
                    ax[i].legend(loc='lower left')
                i += 1

plt.tight_layout()
plt.savefig('./plots/figures/fig6.png')
plt.show()
