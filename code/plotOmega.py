
import matplotlib.pyplot as plt
import numpy as np

savePlots = 0
showPlots = 1

q1,q2,q3,o1,o2,o3= np.genfromtxt("../data/phonons_Kr_omega_100_n=200.data",unpack=True)

plt.figure("100")
plt.plot(q1,o1,label="$\omega_1$")
plt.plot(q1,o2,label="$\omega_2$")
plt.plot(q1,o3,label="$\omega_3$")
plt.legend()
if savePlots:
    plt.savefig("../figs/Kr_omega_100.pdf")

q1,q2,q3,o1,o2,o3= np.genfromtxt("../data/phonons_Kr_omega_110_n=200.data",unpack=True)

plt.figure("110")
plt.plot(q1,o1,label="$\omega_1$")
plt.plot(q1,o2,label="$\omega_2$")
plt.plot(q1,o3,label="$\omega_3$")
plt.legend()
if savePlots:
    plt.savefig("../figs/Kr_omega_110.pdf")

q1,q2,q3,o1,o2,o3= np.genfromtxt("../data/phonons_Kr_omega_111_n=200.data",unpack=True)

plt.figure("111")
plt.plot(q1,o1,label="$\omega_1$")
plt.plot(q1,o2,label="$\omega_2$")
plt.plot(q1,o3,label="$\omega_3$")
plt.legend()
if savePlots:
    plt.savefig("../figs/Kr_omega_111.pdf")
if showPlots:
    plt.show()
