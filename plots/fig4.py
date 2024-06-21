import matplotlib.pyplot as plt
import numpy as np
import re
 
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14


# plt.rc('font', size=BIGGER_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=18)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=22)    # fontsize of the tick labels
plt.rc('ytick', labelsize=22)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=16)  # fontsize of the figure title
plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('legend', fontsize=18)    # legend fontsize


qois = ["VTOT", "T", "C", "Mach", "PT", "mu"]
prefix = "./run_tests/ge_logs/"

n_eb = 20
n_qoi = len(qois)

br = np.zeros([n_qoi,n_eb])
estimate_e = np.zeros([n_qoi,n_eb])
max_e = np.zeros([n_qoi,n_eb])

for i in range(len(qois)):
    qoi = qois[i]
    pmgard_path = prefix + "pmgard_" + qoi + ".log"
    with open(pmgard_path, 'r') as file:
        content = file.read()
    pmgard_error_act = [float(match) for match in re.findall(r'max_act_error = (-?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content)]
    pmgard_error_est = [float(match) for match in re.findall(r'max_est_error = (-?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', content)]
    pmgard_cr = [float(match) for match in re.findall(r'aggregated cr = (-?[\d.]+)', content)]
    br[i] = 64 / np.array(pmgard_cr)
    estimate_e[i] = np.array(pmgard_error_est)
    max_e[i] = np.array(pmgard_error_act)
    
br = np.transpose(br)
estimate_e = np.transpose(estimate_e)
max_e = np.transpose(max_e)
vr = np.array([626.443, 236.349, 127.804, 1.93198, 196480, 1.04299e-05])

a = 0.1
r = 0.5
n = 20
end = a * (r ** (n - 1))
request_e = np.geomspace(a, end, num=n)

styles=['r-o', 'g:^', 'b--*', 'c-o', 'm:^', 'y--*']
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(16,15))
ax=[axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1], axs[2, 0], axs[2, 1]]
labels = ['Requested tolerance', 'Max estimated error', 'Max real error']


i = 0
while i < 6:
                names = ["VTOT", "T", "C", "Mach", "PT", "mu"]
                p, = ax[i].plot(br[:, i], request_e, styles[0], label='{}'.format(labels[0]))
                p, = ax[i].plot(br[:, i], estimate_e[:, i]/vr[i], styles[1], label='{}'.format(labels[1]))
                p, = ax[i].plot(br[:, i], max_e[:, i]/vr[i], styles[2], label='{}'.format(labels[2]))

                ax[i].set_ylabel('Relative QoI Error',fontsize=30)
                ax[i].set_yscale('log')
                ax[i].set_yticks([1e-1, 1e-3, 1e-5, 1e-7, 1e-9])
                # ax_t.set_ylim(0, 2000)
                ax[i].set_xlim(4, 25)
                # ax[i].set_xticks(np.arange(5,26,5))
                # ax_t.set_xlim(0, 30)
                # ax[i].grid(which='major', axis='y')
                ax[i].set_title(names[i], fontsize=30)
                ax[i].set_xlabel('Bitrate',fontsize=30)
                ax[i].legend(loc='upper right')
                if i >= 4:
                    ax[i].legend(loc='lower left')                
                i += 1

plt.tight_layout()
plt.savefig('./plots/figures/fig4.png')
plt.show()