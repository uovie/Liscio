import numpy as np
import matplotlib.pyplot as plt

val = np.fromfile("val.dat", dtype=np.float64)
oval = open('val.txt', 'wt')
for v in val:
    oval.write(f'{v}\n')
oval.close()

crd = np.fromfile("crd.dat", dtype=np.float64)
for i in range(val.size):
    vec = np.fromfile(f"vec_{i}.dat", dtype=np.float64)
    fig, ax = plt.subplots()
    fig.suptitle(f"State {i}")
    ax.plot(crd, vec, c='dodgerblue', lw=0.5, ls='-', alpha=0.8)
    ax.set_xlabel("crd (a.u.)")
    ax.set_ylabel("vec (a.u.)")
    fig.savefig(f'vec_{i}.png', dpi=500, bbox_inches='tight')