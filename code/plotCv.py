
import matplotlib.pyplot as plt
import numpy as np

savePlots = 1
showPlots = 0

T_Ar,CV_Ar= np.genfromtxt("../data/Ar_cv.data",unpack=True)
T_Ne,CV_Ne= np.genfromtxt("../data/Ne_cv.data",unpack=True)
T_Xe,CV_Xe= np.genfromtxt("../data/Xe_cv.data",unpack=True)
T_Kr,CV_Kr= np.genfromtxt("../data/Kr_cv.data",unpack=True)

plt.figure("CV_comparison")
plt.plot(T_Ar,CV_Ar,label="Ar")
plt.plot(T_Ne,CV_Ne,label="Ne")
plt.plot(T_Xe,CV_Xe,label="Xe")
plt.plot(T_Kr,CV_Kr,label="Kr")
plt.legend()
if savePlots:
    plt.savefig("../figs/CV_comparison.pdf")

T_Ar,CV_Ar= np.genfromtxt("../data/Ar_cv_long.data",unpack=True)
T_Ne,CV_Ne= np.genfromtxt("../data/Ne_cv_long.data",unpack=True)
T_Xe,CV_Xe= np.genfromtxt("../data/Xe_cv_long.data",unpack=True)
T_Kr,CV_Kr= np.genfromtxt("../data/Kr_cv_long.data",unpack=True)

plt.figure("CV_comparison_HighT")
plt.plot(T_Ar,CV_Ar,label="Ar")
plt.plot(T_Ne,CV_Ne,label="Ne")
plt.plot(T_Xe,CV_Xe,label="Xe")
plt.plot(T_Kr,CV_Kr,label="Kr")
plt.legend()
if savePlots:
    plt.savefig("../figs/CV_HighT.pdf")
if showPlots:
    plt.show()


A = 10**-10
kb = 1.38064903*10**-23
rNe = 3.1562*A
rAr =  3.7477*A
rKr =   3.9922*A
rXe =    4.3346*A

VNe = (0.2*2**0.5*np.pi/rNe)**3/2
VAr = (0.2*2**0.5*np.pi/rAr)**3/2
VKr = (0.2*2**0.5*np.pi/rKr)**3/2
VXe = (0.2*2**0.5*np.pi/rXe)**3/2

CvNe = 3 * kb * VNe /(2*np.pi)**3
CvAr = 3 * kb * VAr /(2*np.pi)**3
CvKr = 3 * kb * VKr /(2*np.pi)**3
CvXe = 3 * kb * VXe /(2*np.pi)**3

print(CvNe)
print(CV_Ne[-1])
print(CvAr)
print(CV_Ar[-1])
print(CvKr)
print(CV_Kr[-1])
print(CvXe)
print(CV_Xe[-1])
