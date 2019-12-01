import numpy as np
import matplotlib.pyplot as plt

hbar = 1
kB = 1.38064852e-23
hartree = 4.359744650e-18

m = 1728
omega = 0.017709
Ts = [1000, 300, 100, 0]

crd = np.linspace(-1, 1, 200)

fig, ax = plt.subplots(2, 2)
fig.suptitle(f"Density distribution")
for i in range(2):
    for j in range(2):
        if Ts[i * 2 + j] == 0:
            tanh = 1
        else:
            beta = 1 / (kB * Ts[i * 2 + j] / hartree)
            tanh = np.tanh(beta * hbar * omega / 2)
        val = np.sqrt(m * omega * tanh / (np.pi * hbar)) * np.exp(-m * omega * np.power(crd, 2) * tanh / hbar)
        ax[i, j].plot(crd, val, c='dodgerblue', lw=0.5, ls='-', alpha=0.8, label=f'$T = {Ts[i * 2 + j]}$')
        ax[i, j].set_xlabel("$q$ (a.u.)")
        ax[i, j].set_ylabel("$\\rho$ (a.u.)")
        ax[i, j].legend(loc='upper right', prop={'size':6})
fig.subplots_adjust(hspace=0.5, wspace=0.4)
fig.savefig(f'sho_rho_x.png', dpi=500, bbox_inches='tight')