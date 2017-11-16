import numpy as np

#Test for 4b) analytical values
T=1.
J=1.
beta = 1.
N=2.
L=N*N

Z=2*(np.exp(-8*J*beta)+np.exp(8*J*beta)+6)

E=16*J*(np.exp(-8*J*beta)-np.exp(8*J*beta))/Z/L
E_sq = 128*J*(np.exp(-8*J*beta)+np.exp(8*J*beta))/Z/L
M = 8*(np.exp(8*J*beta)+2)/Z/L
M_sq = 32*(np.exp(8*J*beta) +1)/Z/L

Evariance = (E_sq-E**2*L)
Mvariance = (M_sq - M**2*L)

print 'E = ', E, ', E_sq = ', E_sq, ', E_variance', Evariance
print 'M = ', M, ', M_sq = ', M_sq, ', M_variance', Mvariance
Z=4*(np.cosh(8*J*beta) + 3)/L
E=-8*J*np.sinh(8*J*beta)/(np.cosh(8*J*beta) + 3)/L
E_sq = 64*J**2*np.cosh(8*J*beta)/(np.cosh(8*J*beta)+3)/L
M_abs = 2*(np.exp(8*J*beta) +2)/(np.cosh(8*J*beta) +3)/L
M_sq = 8*(np.exp(8*J*beta) +1)/(np.cosh(8*J*beta) + 3)/L

Evariance = (E_sq-E**2*L)*L#64*J**2*(3*np.cosh(8*J*beta)+1)/(np.cosh(8*J*beta) +3)**2/T/L
Mvariance = (M_sq - M_abs**2*L)*L#(8*(np.exp(8*J*beta)+1)/(np.cosh(8*J*beta)+3)-(2*(np.exp(8*J*beta)+2)/(np.cosh(8*J*beta)+3))**2)/T/L
#Mvariance = 4*(2*(np.exp(8*J*beta) +1)*(np.cosh(8*J*beta) + 1)-np.exp(16*J*beta))/(np.cosh(8*J*beta) + 3)**2/T/L

print 'E = ', E, ', E_sq = ', E_sq, ', E_variance', Evariance
print 'M = ', M_abs, ', M_sq = ', M_sq, ', M_variance', Mvariance
