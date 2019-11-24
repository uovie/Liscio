import numpy as np
import matplotlib.pyplot as plt

crd = np.fromfile("crd.dat", dtype=np.float64)
val = np.fromfile("val.dat", dtype=np.float64)

oval = open('val.txt', 'wt')
for v in val:
    oval.write(f'{v}\n')
oval.close()

ocrd = open('crd.txt', 'wt')
for q in crd:
    ocrd.write(f'{q}\n')
ocrd.close()

fig, ax = plt.subplots()
fig.suptitle(f"Density distribution")
ax.plot(crd, val, c='dodgerblue', lw=0.5, ls='-', alpha=0.8)
ax.set_xlabel("q (a.u.)")
ax.set_ylabel("$\\rho$ (a.u.)")
fig.savefig(f'sho_rho.png', dpi=500, bbox_inches='tight')