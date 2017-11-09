import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

def plot_formatting(fam='serif',fam_font='Computer Modern Roman',font_size=14,tick_size=14):
	""" you get to define what font and size of xlabels and axis ticks you"""
	"""like, if you want bold text or not.								  """
	
	plt.rc('text',usetex=True)
	axis_font={'family': fam,'serif':[fam_font],'size':font_size}
	plt.rc('font',**axis_font)
	plt.rc('font',weight ='bold')
	plt.xticks(fontsize=tick_size)
	plt.yticks(fontsize=tick_size)

plot_formatting()


N1 = np.genfromtxt('accepted_changes_c_T24_fixed1.txt')
print 'read1'
N2 = np.genfromtxt('accepted_changes_c_T24_random1.txt')
print 'read2'
MC = np.linspace(0,10**7,len(N1))

count1 = 0
count2 = 0
accepted1 = np.zeros(len(N1))
accepted2 = np.zeros(len(N1))

for i in range(len(N1)):
	accepted1[i] = count1
	accepted2[i] = count2
	count1 += N1[i]
	count2 += N2[i]


plt.plot(MC[:2000],accepted1[:2000])
plt.plot(MC[:2000],accepted2[:2000])
plt.xlabel('Number of MC cycles')
plt.ylabel('Accepted configurations')
plt.legend(['Ordered initial spin state, T=2.4','Random initial spin state, T=2.4'],loc='best')
plt.show()

'''
#N3 = np.genfromtxt('accepted_changes_c_T24_random.txt')
N4 = np.genfromtxt('probability_d_T24_random.txt')
N4 = N4[0::4]/sum(N4)
N6 = np.genfromtxt('Evariance_d_T24_random.txt')
#N2 = range(0,len(N3)) 
indx = np.argmax(N4)
print N6

N5 = np.linspace(-800,800,len(N4))

mean = N5[indx]
#should be shifted slightly to the right, as not a perfect gaussian, is only perfect if no correlation
Y=scipy.stats.norm.pdf(N5,loc=mean,scale=np.sqrt(N6*400))
plt.plot(N5/20**2,N4)
plt.plot(N5/20**2,Y/sum(Y))
plt.xlabel(r'$\langle E\rangle$ per spin')
plt.ylabel(r'$P(E)$')
plt.legend(['Simulation, $T=2.4$','Normal distribution, $T=2.4$'])
#plt.axis([0,100, -2,-.4])
#plt.xlabel('Number of MC cycles')
#plt.ylabel('Energy')
#plt.hist(N1[20:],bins=np.arange(min(N1[20:]),max(N1[20:]) + 0.0001, 0.0001))
plt.show()
'''