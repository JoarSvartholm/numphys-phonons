
import matplotlib.pyplot as plt
import numpy as np

savePlots = 0
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

VNe = (2*rNe/2**0.5)**3
VAr = (2**0.5*rAr)**3
VKr = (2**0.5*rKr)**3
VXe = (2**0.5*rXe)**3

CvNe = 3 * kb * 4 / VNe
CvAr = 3 * kb * 4 / VAr
CvKr = 3 * kb * 4 / VKr
CvXe = 3 * kb * 4 / VXe

print(CvNe)
print(CV_Ne[-1])
print((CvNe-CV_Ne[-1])/CvNe)
print(CvAr)
print(CV_Ar[-1])
print((CvAr-CV_Ar[-1])/CvAr)
print(CvKr)
print(CV_Kr[-1])
print((CvKr-CV_Kr[-1])/CvKr)
print(CvXe)
print(CV_Xe[-1])
print((CvXe-CV_Xe[-1])/CvXe)
