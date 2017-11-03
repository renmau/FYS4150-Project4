import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from scipy.interpolate import interp1d

T40 = np.genfromtxt('temperature_e_40.txt')
Cv40 = np.genfromtxt('Cv_values_e_40.txt')
Cv60 = np.genfromtxt('Cv_values_e_60.txt')
Cv80 = np.genfromtxt('Cv_values_e_80.txt')
Cv100 = np.genfromtxt('Cv_values_e_100.txt')

M40 = np.genfromtxt('Mabs_values_e_40.txt')
M60 = np.genfromtxt('Mabs_values_e_60.txt')
M80 = np.genfromtxt('Mabs_values_e_80.txt')
M100 = np.genfromtxt('Mabs_values_e_100.txt')

f40 = interp1d(T40,Cv40, kind='quadratic')
f60 = interp1d(T40,Cv60, kind='quadratic')
f80 = interp1d(T40,Cv80, kind='quadratic')
f100 = interp1d(T40,Cv100, kind='quadratic')
T_new = np.linspace(2.2,2.3,100000)

plt.plot(T_new,f40(T_new))
plt.plot(T_new,f60(T_new))
plt.plot(T_new,f80(T_new))
plt.plot(T_new,f100(T_new))
plt.show()

maxT40 = T_new[np.argmax(f40(T_new))]
maxT60 = T_new[np.argmax(f60(T_new))]
maxT80 = T_new[np.argmax(f80(T_new))]
maxT100 = T_new[np.argmax(f100(T_new))]

a1 = (maxT40-maxT60)/(1./40 - 1./60)
a2 = (maxT40-maxT80)/(1./40 - 1./80)
a3 = (maxT40-maxT100)/(1./40 - 1./100)
a4 = (maxT60-maxT80)/(1./60 - 1./80)
a5 = (maxT80-maxT100)/(1./80 - 1./100)
a6 = (maxT80-maxT40)/(1./80 - 1./40)
a7 = (maxT80-maxT60)/(1./80 - 1./60)
a8 = (maxT60-maxT40)/(1./60 - 1./40)
a9 = (maxT60-maxT100)/(1./60 - 1./100)
a10 = (maxT100-maxT40)/(1./100 - 1./40)
a11 = (maxT100-maxT60)/(1./100 - 1./60)
a12 = (maxT100-maxT80)/(1./100 - 1./80)


a = (a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12)/12.
#a = (a1+a2+a4+a6+a7+a8)/6.
print a

'''
plt.plot(T40, Cv40)
plt.plot(T40, Cv60)
plt.plot(T40, Cv80)
plt.plot(T40, Cv100)
plt.show()
'''
plt.plot(T40, M40)
plt.plot(T40, M60)
plt.plot(T40, M80)
plt.plot(T40, M100)
plt.show()

nu = 1.
L = np.array([40.,60.,80.,100.])
T_inf = np.zeros(len(L))
Cv = [Cv40, Cv60, Cv80, Cv100]
Tc = [maxT40,maxT60,maxT80,maxT100]

for i in range(len(L)):
	T_inf[i] = Tc[i] -a*L[i]**(-1./nu)

print np.mean(T_inf)

'''
#N3 = np.genfromtxt('accepted_changes_c_T24_random.txt')
N4 = np.genfromtxt('probability_d_T24_random.txt')
N4 = N4[0::4]/sum(N4)
N6 = np.genfromtxt('Evariance_d_T24_random.txt')
#N2 = range(0,len(N3)) 
indx = np.argmax(N4)

N5 = np.linspace(-800,800,len(N4))

mean = N5[indx]
#should be shifted slightly to the right, as not a perfect gaussian
Y=scipy.stats.norm.pdf(N5,loc=mean,scale=np.sqrt(N6*400))
plt.plot(N5,N4)
plt.plot(N5,Y/sum(Y))
#plt.axis([0,100, -2,-.4])
#plt.xlabel('Number of MC cycles')
#plt.ylabel('Energy')
#plt.hist(N1[20:],bins=np.arange(min(N1[20:]),max(N1[20:]) + 0.0001, 0.0001))
plt.show()
'''