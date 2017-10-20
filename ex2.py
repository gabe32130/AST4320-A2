import pylab as pl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.stats
from scipy.special import erf

def sigmasq(Sc):
    return np.pi / Sc**4

def Pnc(x, sigma = 1.0, deltacrit = 1.0):
    f = scipy.stats.norm.pdf
    return f(x, scale = sigma) - f(2 * deltacrit - x, scale = sigma)

    
N = 100000 #Number of random walks
Sc0 = ( np.pi / 1e-4 )**0.25 #sigma(Sc)^2 < 1e-4
epsilon = 0.99 #Epsilon factor to reduce steps

#Compute number of steps needed for Sc ~ 1
#solve Sc0 epsilon^n = 1 for n
nsteps = np.ceil( -np.log(Sc0)/np.log(epsilon) )
D=nsteps + 1
zeromatrix=(259, 100000)


delta = np.zeros(zeromatrix)
i = 0
Sc = Sc0 * epsilon**np.arange(D)
delta[0,:] = np.random.normal(size=N, scale=np.sqrt(sigmasq(Sc0)))
#(S1,del)


while(i < nsteps):
    sigma12 = np.sqrt( sigmasq(Sc[i + 1]) - sigmasq(Sc[i]) )
    beta = np.random.normal(size = N, scale = sigma12 )
    delta[i + 1, :] = delta[i,:] + beta #new del2 (S2,del2)
    i += 1


plt.clf()
plt.hist(delta[i,:], 100, histtype='stepfilled', normed=True, label=r"$\delta( Sc = 1.0)$")
xlst = np.linspace(plt.xlim()[0],plt.xlim()[1], 1000)
ylst = scipy.stats.norm.pdf(xlst, scale = np.sqrt(sigmasq(Sc[i])))
plt.plot(xlst,ylst, label=r"$N(\sigma=1)$")
#plt.xlabel(r"$\delta$")
#plt.legend(loc='best')
#plt.show()


#only the ones that never cross the threshold
deltacrit = 1.0
notcrossed = np.max(delta, axis=0) < deltacrit


plt.hist(delta[i, notcrossed], 100, histtype='stepfilled', normed=True, label=r"$\delta_{nc} (Sc = 1.0)$")
xlst = np.linspace(plt.xlim()[0], deltacrit, 1000)
ylst = Pnc(xlst, sigma = np.sqrt(sigmasq(Sc[i])), deltacrit = deltacrit)
normf = 1.0/erf(deltacrit / np.sqrt(2 * np.pi))
plt.plot(xlst,ylst * normf, label=r"$P_{nc}(\delta)$")
plt.xlabel(r"$\delta$")
plt.legend(loc='best')
plt.savefig("fig.pdf")
plt.show()
