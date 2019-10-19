import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.interpolate import interp1d


def spectral():

    A=np.array([219/256, 55/256, 55/256, 1])	
    B=np.array([230/256, 75/256, 60/256, 1])	
    C=np.array([244/256, 109/256, 67/256, 1])	
    D=np.array([253/256, 174/256, 97/256, 1])	
    E=np.array([254/256, 224/256, 139/256, 1])	
    F=np.array([230/256, 245/256, 152/256, 1])	
    G=np.array([171/256, 221/256, 164/256, 1])	
    H=np.array([102/256, 194/256, 165/256, 1])	
    I=np.array([50/256, 136/256, 189/256, 1])	
    J=np.array([94/256, 79/256, 162/256, 1])
    K=np.array([46/256, 67/256, 92/256, 1])
    L=np.array([38/256, 42/256, 77/256, 1])
    M=np.array([37/256, 43/256, 61/256, 1])
    N=np.array([33/256, 38/256, 51/256, 1])
    
    ncmap = [N, M, L, K, J, I, H, G, F, E, D, C, B, A]
    c_array = []
    
    for i in range(len(ncmap)):
        if i != len(ncmap)-1:
            linfit = interp1d([1,256], np.vstack([ncmap[i], ncmap[i+1]]), axis=0)
            for j in range(255):
                c_array.append(linfit(j+1))
    
    newcmp = ListedColormap(c_array)
    return newcmp

if __name__=="__main__":
    
    newcmp = spectral()
    
    def sinus2d(x, y):
        return np.sin(x) + np.sin(y)
    
    xx, yy = np.meshgrid(np.linspace(0,3*np.pi,512), np.linspace(0,3*np.pi,512))
    z = sinus2d(xx, yy)
    plt.imshow(z, origin='lower', interpolation='none', cmap=newcmp)
    cbar = plt.colorbar()
    plt.show()
