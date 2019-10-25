import numpy as np
import matplotlib.pyplot as plt

time = np.fromfile("sho.time", dtype=np.float64)
etot = np.fromfile("sho.etot", dtype=np.float64)
pos = np.fromfile("sho.pos", dtype=np.float64)
mom = np.fromfile("sho.mom", dtype=np.float64)

fig, ax = plt.subplots()
fig.suptitle('Energy Evolution')
ax.plot(time, etot, c='dodgerblue', lw=0.5, ls='-', alpha=0.8)
ax.set_xlabel("t (a.u.)")
ax.set_ylabel("Etot (a.u.)")
fig.savefig('eto_ev.png', dpi=500, bbox_inches='tight')

fig, ax = plt.subplots()
fig.suptitle('Phase Space Trajectory')
ax.plot(pos, mom, c='dodgerblue', lw=0.5, ls='-', alpha=0.8)
ax.set_xlabel("q (a.u.)")
ax.set_ylabel("p (a.u.)")
fig.savefig('phs_tr.png', dpi=500, bbox_inches='tight')