import numpy as np
from math import ceil
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


#params TO SET BEFORE RUNNING
n = 5
L = 2.3
timesteps = 1000
numP = n*n*n



#sphere coordinates
u = np.linspace(0, np.pi, 30)
v = np.linspace(0, 2 * np.pi, 30)

xs = L*np.outer(np.sin(u), np.sin(v))
ys = L*np.outer(np.sin(u), np.cos(v))
zs = L*np.outer(np.cos(u), np.ones_like(v))


# ax.plot_wireframe(x, y, z, rstride=3, cstride=3, linewidth=1)
# ax.plot(atoms[:, 0], atoms[:, 1], atoms[:, 2], 'mo', markersize=5)

#argon crystal
with open('./avs_T500.dat', 'r') as f:
    particleData = []
    for line in f:
        line = line.split()
        particleData.append(line)

x = [float(item[0]) for item in particleData]
y = [float(item[1]) for item in particleData]
z = [float(item[2]) for item in particleData]


#time, temperature, pressure data
data = np.loadtxt('./out_T500.txt', skiprows=1)
T = data[:, 3]
P = data[:, 4]
t = data[:, 0]



fig = plt.figure()
ax = p3.Axes3D(fig)


border = ceil(L)
ax.set_xlim3d([-border, border])
ax.set_ylim3d([-border, border])
ax.set_zlim3d([-border, border])


def animate(i):
    ax.clear()
    ax.set_xlim3d([-border, border])
    ax.set_ylim3d([-border, border])
    ax.set_zlim3d([-border, border])
    idx0 = i*numP
    idx1 = numP*(i+1)
    ax.plot_wireframe(xs, ys, zs, rstride=3, cstride=3, linewidth=1)
    ax.scatter(x[idx0:idx1],y[idx0:idx1],z[idx0:idx1], c = 'm', s = 48, label = "T = " + str(T[i]) + "\nP = " + str(P[i]) + "\nt = " + str(t[i]))
    ax.legend()

ani = animation.FuncAnimation(fig, animate, frames=timesteps, interval=50, blit=False, repeat=False)
plt.show()
