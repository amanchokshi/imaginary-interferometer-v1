import numpy as np

print("E,N")
r = 100
#n = 96
n = 64

for i in range(n):
    x = round(r * np.cos(2*np.pi*i/n))
    y = round(r * np.sin(2*np.pi*i/n))

    print(x,',',y, sep='')
