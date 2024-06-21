import matplotlib.pyplot as plt
import numpy as np
import re
# from os import system

 
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
plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
plt.rc('legend', fontsize=14)    # legend fontsize


labels = ['Requested tolerance', 'Max estimated error', 'Max real error']
data = ["NYX", "Hurricane"]
prefix = ["./run_tests/NYX_logs/", "./run_tests/Hurricane_logs/"]

a = 0.1
r = 0.1
n_eb = 7
end = a * (r ** (n_eb - 1))
request_e = np.geomspace(a, end, num=n_eb)
n_qoi = 2

br = np.zeros([n_qoi,n_eb])
estimate_e = np.zeros([n_qoi,n_eb])
max_e = np.zeros([n_qoi,n_eb])

for i in range(len(data)):
    pmgard_path = prefix[i] + "pmgard_VTOT_" + data[i] + ".log"
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
vr = np.array([6.05364e+07, 56.286])

print(request_e)
print(max_e[:,1])

styles=['r-o', 'g:^', 'b--*', 'c-o', 'm:^', 'y--*']
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12,4))
ax=[axs[0], axs[1]]

i = 0
while i < 2:
                names = ["VTOT - NYX", "VTOT - Hurricane"]
                p, = ax[i].plot(br[:, i], request_e, styles[0], label='{}'.format(labels[0]))
                p, = ax[i].plot(br[:, i], estimate_e[:, i]/vr[i], styles[1], label='{}'.format(labels[1]))
                p, = ax[i].plot(br[:, i], max_e[:, i]/vr[i], styles[2], label='{}'.format(labels[2]))

                ax[i].set_ylabel('Relative VTOT Error')
                ax[i].set_yscale('log')
                # ax_t.set_ylim(0, 2000)
                if i == 0:
                    ax[i].set_xlim(4, 20)
                else:
                    ax[i].set_xlim(4, 20)
                # ax_t.set_xlim(0, 30)
                # ax[i].grid(which='major', axis='y')
                ax[i].set_yticks([1e-1, 1e-3, 1e-5, 1e-7])
                ax[i].set_title(names[i],fontsize=22)
                ax[i].set_xlabel('Bitrate',fontsize=22)
                ax[i].legend(loc='upper right')
                i += 1

plt.tight_layout()
plt.savefig('./plots/figures/fig5.png')
plt.show()
