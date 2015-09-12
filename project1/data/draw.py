import numpy as np
import matplotlib.pyplot as plt

nbs = np.arange(4, 12)          # nbs
nbs3 = nbs**3
nbs4 = nbs**4
E = np.zeros((8, 5))            # average
V = np.zeros((8, 5))            # standard deviation
for i in range(8):
    data = np.loadtxt('nbs_' + str(i+4) + '.txt')
    E[i] = np.mean(data, axis=0)
    V[i] = np.std(data, axis=0)

case = 4

if case == 1:
    """
    nbs vs time plot
    """
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax.errorbar(nbs, E[:, 0], yerr=V[:, 0], fmt='-o', label=r'$t_{01}$')
    ax.errorbar(nbs, E[:, 1], yerr=V[:, 1], fmt='-s', label=r'$t_{12}$')
    ax.errorbar(nbs, E[:, 2], yerr=V[:, 2], fmt='-^', label=r'$t_{23}$')
    ax.errorbar(nbs, E[:, 3], yerr=V[:, 3], fmt='-+', label=r'$t_{34}$')
    ax.set_xlabel('nbs')
    ax.set_ylabel('time')
    fig.tight_layout(pad=0)
    plt.legend(loc='best')
    plt.show()


if case == 2:
    """
    nbs^3 vs time plot
    """
    
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax.errorbar(nbs3, E[:, 0], yerr=V[:, 0], fmt='-o', label=r'$t_{01}$')
    ax.errorbar(nbs3, E[:, 1], yerr=V[:, 1], fmt='-s', label=r'$t_{12}$')
    ax.errorbar(nbs3, E[:, 2], yerr=V[:, 2], fmt='-^', label=r'$t_{23}$')
    ax.errorbar(nbs3, E[:, 3], yerr=V[:, 3], fmt='-+', label=r'$t_{34}$')
    ax.set_xlabel('number of balls (nbs*nbs*nbs)')
    ax.set_ylabel('time')
    fig.tight_layout(pad=0)
    plt.legend(loc='best')
    plt.show()

if case == 3:
    """
    collision time per frame
    """
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax.errorbar(nbs, E[:, 4], yerr=V[:, 4], fmt='-o', label=r'collisions per frame')
    ax.set_xlabel('nbs')
    ax.set_ylabel('collisions/frame')
    ax.grid('on')
    fig.tight_layout(pad=0)
    plt.legend(loc='best')
    plt.show()

if case == 4:
    """
    forth order fit curve
    """
    p = np.polyfit(nbs4, E[:, 4], 1)
    fig = plt.figure(figsize=(6, 4))
    x = np.arange(4, 12)**4
    ax = fig.add_subplot(111)
    ax.plot(nbs4, E[:, 4], '-o', lw=1.5, label=r'experimental result')
    ax.plot(x, p[1]+p[0]*x, 'r--', lw=2, label='analytical relation')
    ax.set_xlabel(r'$nbs^4$')
    ax.set_ylabel('collisions/frame')
    ax.grid('on')
    fig.tight_layout(pad=0)
    plt.legend(loc='best')
    plt.show()
