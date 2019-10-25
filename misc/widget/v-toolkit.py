### Versatile Toolkit ###
import numpy as np

''' Contants '''
h_bar = 1

''' SHO: Determine Suitable Well Width '''
mass = 1
omega = np.sqrt(2)
nbasis = 2000
nstate = 20
width_max = np.sqrt((nbasis * np.pi * h_bar)**2 / (2 * mass) / ((nstate - 0.5) * h_bar * omega))
print(f'width < {width_max}')

''' SHO: Exact Energy Levels '''
for i in range(nstate):
    print((i + 0.5) * h_bar * omega)


''' ISW: Determine Suitable Angular Frequency '''
mass = 1
width = 2
nbasis = 5000
nstate = 20
omega_min = (nstate * np.pi * h_bar)**2 / (2 * mass * width**2) / ((nbasis - 0.5) * h_bar)
print(f'omega > {omega_min}')

''' ISW: Exact Energy Levels '''
for i in range(nstate):
    print(((i + 1) * np.pi * h_bar)**2 / (2 * mass * width**2))


''' MRS: Change units to a.u. '''
a0 = 5.2917721067e-11
Avo = 6.022140857e23
Cal = 4184
print(f'De = {116.09 * Cal / Avo}')
print(f'alpha = {2.287e10 * a0}')
print(f'xeq = {9.419e-11 / a0}')

''' MRS: Determine Suitable Mass '''
De = 8.0656e-19
xeq = 1.7799
nbasis = 150
nstate = 5
mass = (nbasis * np.pi * h_bar)**2 / (2 * 20 **2 * De)
print(f'mass > {mass}')


''' EKT: Change units to a.u. '''
eV = 1.6021766208e-19
Ht = 4.359744650e-18 # Hartree energy
print(f'V0 = {0.425 * eV / Ht}')

''' MRS: Determine Suitable Mass '''
V0 = 0.0156
c = 1.3624
nbasis = 150
nstate = 5
mass = (nbasis * np.pi * h_bar)**2 / (2 * 10 **2 * V0)
print(f'mass > {mass}')