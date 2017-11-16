import matplotlib as mpl
mpl.use('TkAgg') # Ensure that the Tkinter backend is used for generating figures
import matplotlib.pyplot as plt
import numpy as np

def plot_formatting(fam='serif',fam_font='Computer Modern Roman',font_size=12,tick_size=14):
	""" you get to define what font and size of xlabels and axis ticks you"""
	"""like, if you want bold text or not.								  """
	
	plt.rc('text',usetex=True)
	axis_font={'family': fam,'serif':[fam_font],'size':font_size}
	plt.rc('font',**axis_font)
	plt.rc('font',weight ='bold')
	plt.xticks(fontsize=tick_size)
	plt.yticks(fontsize=tick_size)

plot_formatting()

E = np.array([-8.0,-8.0,-8.0,-7.994,-7.9898,-7.98372,-7.98406,-7.984])
Mabs = np.array([4.0,4.0,4.0,3.9984,3.99489,3.99459,3.99467,3.99466])
Cv = np.array([0,0,0,0.0447686,0.121369,0.129975,0.127298,0.127737])
chi = np.array([0,0,0,0.00399744,0.0157334,0.0161388,0.0160076,0.0159878])
Esq = Cv+E*2
Msq = chi +Mabs**2

E_real = -7.9839283437467605
Esq_real = 63.8714411255
Msq_real = 15.9732151042
Mabs_real = 3.9946429309943987
Cv_real = 0.1283006340818
chi_real = 0.01604295806491


rel_E = abs((E-E_real)/E_real)
rel_Esq = abs(Esq-Esq_real)/Esq_real
rel_Mabs = abs(Mabs-Mabs_real)/Mabs_real
rel_Msq = abs(Msq-Msq_real)/Msq_real
rel_Cv = abs(Cv-Cv_real)/Cv_real
rel_chi = abs(chi-chi_real)/chi_real

MC = np.array([10,100,1000,10000,100000,1000000,10000000,100000000])

ax1=plt.subplot(2,2,1)
plt.semilogx(MC,rel_E,'x-')
#plt.xlabel('Monte Carlo cycles')
plt.ylabel(r'Relative error for $\langle E\rangle$')
plt.setp(ax1.get_xticklabels(), visible=False)
#plt.annotate('$\lambda = $'+np.str(int(wav[k]*1e8))+' \AA', xy=(0.75, 0.12), xycoords='axes fraction')

#plt.show()
ax2=plt.subplot(2,2,2)
plt.semilogx(MC,rel_Mabs,'x-')
#plt.xlabel('Monte Carlo cycles')
plt.ylabel(r'Relative error for $\langle |\mathcal{M}|\rangle$')
plt.setp(ax2.get_xticklabels(), visible=False)
#plt.show()

ax3=plt.subplot(2,2,3)
plt.semilogx(MC,rel_Cv,'x-')
plt.xlabel('Monte Carlo cycles')
plt.ylabel(r'Relative error for $C_V$')
#plt.setp(ax.get_xticklabels(), visible=False)
#plt.show()

ax4=plt.subplot(2,2,4)
plt.semilogx(MC,rel_chi,'x-')
plt.xlabel('Monte Carlo cycles')
plt.ylabel(r'Relative error for $\mathcal{X}$')
plt.subplots_adjust(hspace=.1)
plt.tight_layout()
plt.show()