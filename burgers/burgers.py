from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot as plt
import numpy as np
import os

data = np.loadtxt("salida.dat")

nx = 41
ny = 41
x = np.linspace(0,2,nx)
y = np.linspace(0,2,ny)

nt = 500

for i in range(0, (nt + 1)):
    
    k = i * ny
    
    u = data[k:(k + ny), 0:nx]
    v = data[k:(k + ny), nx:(2 * nx)]

    fig = plt.figure(figsize = (11,7), dpi = 100)
    ax = fig.gca(projection = '3d')
    X,Y = np.meshgrid(x,y)

    ax.set_zlim(1, 2.0)

    wire1 = ax.plot_wireframe(X,Y,u[:], cmap=cm.coolwarm)
    wire2 = ax.plot_wireframe(X,Y,v[:], cmap=cm.coolwarm)
    
    print i

    if i < 10:

        plt.savefig("Burgers00" + str(i) + ".png")
    
        plt.close(fig)

    elif (i >= 10) and (i < 100):

        plt.savefig("Burgers0" + str(i) + ".png")
    
        plt.close(fig)

    else:
    
        plt.savefig("Burgers" + str(i) + ".png")
    
        plt.close(fig)

os.system("convert -delay 20 -loop 0 Burgers*.png burgers.gif")
