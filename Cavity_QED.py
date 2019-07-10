# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 16:21:30 2018

@author: Administrator
"""

#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
N = 7

wr = 2.0*10**1    # resonator frequency
wq = 2.0*10**1     # qubit frequency
chi = 2   # parameter in the dispersive hamiltonian

delta = abs(wr - wq)        # detuning
g = 20  # coupling strength that is consistent with chi

# cavity operators
a = tensor(destroy(N), qeye(2))
nc = a.dag() * a
xc = a + a.dag()

# atomic operators
sm = tensor(qeye(N), destroy(2))
sz = tensor(qeye(N), sigmaz())
sx = tensor(qeye(N), sigmax())
nq = sm.dag() * sm
xq = sm + sm.dag()

I = tensor(qeye(N), qeye(2))
# dispersive hamiltonian
H = wr * (a.dag() * a + I/2.0) - (wq / 2.0) * sz + g* (a.dag()+ a) * sx

#psi0 = tensor(coherent(N, sqrt(6)), (basis(2,0)+basis(2,1)).unit())
#psi0 = tensor(thermal_dm(N, 3), ket2dm(basis(2,0)+basis(2,1))).unit()
psi0 = tensor(coherent(N, np.sqrt(4)), (basis(2,0)+basis(2,1)).unit())

tlist = np.linspace(0, 10, 100001)
corr_vec = correlation(H, psi0, None, tlist, [], a.dag(), a)


w, S = spectrum_correlation_fft(tlist, corr_vec)
fig, ax = plt.subplots(figsize=(9,3))
ax.plot(w / (2), abs(S))
ax.set_xlabel(r'$\omega$', fontsize=18)
#ax.set_xlim(wr/(2)-.5, wr/(2)+.5);