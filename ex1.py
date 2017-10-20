import pylab as pl
import numpy as np
import cmath as m
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.interpolate import UnivariateSpline
import pylab as pl

#plot the step function
step=1000
x=np.linspace(-10, 10, step)
xn=np.zeros(len(x))
xp=np.zeros(len(x))
Wx=np.zeros(len(x))
Wxn=np.zeros(len(x))
Wxp=np.zeros(len(x))
R=6.5

for i in range (len(x)):
    if x[i] <0:
        xn[i]=x[i]
        if abs(xn[i]) < R:
            Wxn[i]=1
        else:
            Wxn[i]=0
    else:
        xn[i]=0

for i in range (len(x)):
    if x[i] >0:
        xp[i]=x[i]
        if abs(xp[i]) < R:
            Wxp[i]=1
        else:
            Wxp[i]=0
    else:
        xp[i]=0
x= xn+xp
Wx=Wxn+Wxp


plt.plot(x,Wx, label=r'linewidth')
plt.xlabel(r'x', size=14)
plt.ylabel(r'W(x)', size=14)
plt.ylim([0,2])
plt.legend(fontsize=14)
plt.savefig("fig1.pdf",bbox_inches='tight')
plt.show()


################################################################################
#Fourier Transform
W_f=np.zeros(len(x))
k=x
W_f = np.sin(2.0*R*k)/(2.0*np.pi*k)


plt.plot(x,W_f, label=r'linewidth')
plt.xlabel(r'x', size=14)
plt.ylabel(r'W(f)', size=14)
plt.ylim([-0.5,2.5])
plt.legend(fontsize=14)
plt.savefig("fig2.pdf",bbox_inches='tight')
plt.show()


################################################################################
#FWHM

half_max=np.max(W_f)/2
print (half_max)
#max_x = x[W_f.index(half_max)]
#print (max_x)
#indx=x.index(-0.14695)
#print (indx)


x_curve = UnivariateSpline(x, W_f, s=0)
r=x_curve.roots()
L=len(r)
#print (L)
max= (L/2)-2
min= (L/2)-1
r1=r[40]
r2=r[41]
FWHM=abs(r1-r2)
print(FWHM)


pl.plot(x, W_f)
pl.axvspan(r1, r2, facecolor='g', alpha=0.5)
plt.savefig("fig3.pdf",bbox_inches='tight')
pl.show()

#-0.14695
