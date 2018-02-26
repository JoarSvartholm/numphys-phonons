
import matplotlib.pyplot as plt
import numpy as np

savePlots = 1
showPlots = 0

q1,q2,q3,g1,g2,g3= np.genfromtxt("../data/Ne_gamma_100.data",unpack=True)

plt.figure("100")
plt.plot(q1,g1,label="$\gamma_1$")
plt.plot(q1,g2,'--',label="$\gamma_2$")
plt.plot(q1,g3,label="$\gamma_3$")
plt.ylabel("Gruneisen parameter [dimensionless]")
plt.xlabel("Normalized q-vector")
plt.legend()
if savePlots:
    plt.savefig("../figs/Ne_gamma_100.pdf")

q1,q2,q3,o1,o2,o3= np.genfromtxt("../data/Ne_gamma_110.data",unpack=True)

plt.figure("110")
plt.plot(q1,o1,label="$\gamma_1$")
plt.plot(q1,o2,'--',label="$\gamma_2$")
plt.plot(q1,o3,label="$\gamma_3$")
plt.ylabel("Gruneisen parameter [dimensionless]")
plt.xlabel("Normalized q-vector")
plt.legend()
if savePlots:
    plt.savefig("../figs/Ne_gamma_110.pdf")

q1,q2,q3,o1,o2,o3= np.genfromtxt("../data/Ne_gamma_111.data",unpack=True)

plt.figure("111")
plt.plot(q1,o1,label="$\gamma_1$")
plt.plot(q1,o2,'--',label="$\gamma_2$")
plt.plot(q1,o3,label="$\gamma_3$")
plt.ylabel("Gruneisen parameter [dimensionless]")
plt.xlabel("Normalized q-vector")
plt.legend()
if savePlots:
    plt.savefig("../figs/Ne_gamma_111.pdf")
if showPlots:
    plt.show()
