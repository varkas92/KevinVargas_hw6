from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import re
import matplotlib.pyplot as plt
import sys

E = sys.argv[1]

alpha = sys.argv[2]

def coor(K, a):
    
    s = str(K)
    
    E = s
    
    if "e-" in s:
        
        E2 = s.split("e-")
    
        E = format(K, "." + E2[1] + "f")
    
    alpha = str(a)
    
    data = np.loadtxt("trayectoria_" + E + "_" + alpha + ".dat")

    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]

t, x, y, z = coor(E, alpha)

plt.figure(figsize = (10, 10))
plt.plot(x, y)
plt.title("$xy$")
plt.xlabel("$x (km)$", size = 20)
plt.ylabel("$y (km)$", size = 20)
            
plt.savefig("particle_xy.pdf")
            
plt.close()

fig = plt.figure(figsize = (11, 7), dpi = 100)
ax = fig.gca(projection = "3d")
ax.set_title("$yzx$", size = 30)
ax.set_xlabel("$y (km)$", size = 20)
ax.set_ylabel("$z (km)$", size = 20)
ax.set_zlabel("$x (km)$", size = 20)
wire1 = ax.plot_wireframe(y, z, x, cmap=cm.coolwarm)

plt.savefig("particle_xyz.pdf")
            
plt.close()
