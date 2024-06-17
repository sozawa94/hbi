import math
import numpy as np

# Constants and parameters
n = 9250
pi = math.pi
ds0 = 1.0
mu0 = 0.6
vref = 1e-6
vs = 3.464
velinit = 1e-9
rigid = 32.04

# Arrays
xcol = np.zeros(n)
ycol = np.zeros(n)
zcol = np.zeros(n)
ang = np.zeros(n)
angd = np.zeros(n)
rake = np.zeros(n)
xr = np.zeros(n)
yr = np.zeros(n)
zr = np.zeros(n)
ds = np.zeros(n)
a = np.zeros(n)
b = np.zeros(n)
dc = np.zeros(n)
f0 = np.zeros(n)
vel = np.zeros(n)
sigma = np.zeros(n)
tau = np.zeros(n)
mu = np.zeros(n)
taudot = np.zeros(n)
sigmadot = np.zeros(n)
xs1 = np.zeros(n)
xs2 = np.zeros(n)
xs3 = np.zeros(n)
ys1 = np.zeros(n)
ys2 = np.zeros(n)
ys3 = np.zeros(n)
zs1 = np.zeros(n)
zs2 = np.zeros(n)
zs3 = np.zeros(n)

# Geometry
with open('bp5t.stl', 'r') as file:
    lines = file.readlines()[1:]

k = 0
for i in range(0, len(lines), 7):
    if i + 6 >= len(lines):
        break
    k += 1
    if k >= n:
        break
    # print(lines[i + 3].split())
    # print(lines[i + 3].split()[-3])
    xs1[k] = float(lines[i + 2].split()[-3])
    xs1[k] = float(lines[i + 2].split()[-3])
    ys1[k] = float(lines[i + 2].split()[-2])
    zs1[k] = float(lines[i + 2].split()[-1])
    xs2[k] = float(lines[i + 3].split()[-3])
    ys2[k] = float(lines[i + 3].split()[-2])
    zs2[k] = float(lines[i + 3].split()[-1])
    xs3[k] = float(lines[i + 4].split()[-3])
    ys3[k] = float(lines[i + 4].split()[-2])
    zs3[k] = float(lines[i + 4].split()[-1])
    
    xcol[k] = (xs1[k] + xs2[k] + xs3[k]) / 3
    ycol[k] = (ys1[k] + ys2[k] + ys3[k]) / 3
    zcol[k] = (zs1[k] + zs2[k] + zs3[k]) / 3
    
    #print(xcol[k], ycol[k], zcol[k])

print(k-1)

# Frictional parameters
a0 = 0.004
b0 = 0.03
dc0 = 0.14
a_max = 0.04

with open('bp5t_param.dat', 'w') as file:
    for i in range(n):
        f0[i] = mu0
        dep = -zcol[i]

        if 4 < dep < 16 and abs(xcol[i]) < 30:
            a[i] = a0
        elif dep < 2 or dep > 18 or abs(xcol[i]) > 32:
            a[i] = a_max
        else:
            r = max(abs(dep - 10) - 6, abs(xcol[i]) - 30) / 2
            a[i] = a0 + r * (a_max - a0)
        
        r = max(abs(dep - 10) - 6, abs(xcol[i]) - 30) / 2
        a[i] = min(a0 + r * (a_max - a0), a_max)
        a[i] = max(a[i], a0)

        b[i] = b0
        dc[i] = dc0

        if abs(xcol[i] + 24) < 6 and abs(dep - 10) < 6:
            dc[i] = 0.13

        sigma[i] = 25.0
        vel[i] = velinit
        dep = -zcol[i]
        if abs(xcol[i] + 24) < 6 and abs(dep - 10) < 6:
            vel[i] = 3e-2
        
        mu[i] = a[i] * math.asinh(0.5 * vel[i] / vref * math.exp((f0[i] + b[i] * math.log(vref / velinit)) / a[i])) + rigid / (2 * vs) * vel[i]
        tau[i] = mu[i] * sigma[i]

        taudot[i] = 0.0
        sigmadot[i] = 0.0

        file.write(f'{0.0:17.8e} {a[i]:17.8e} {b[i]:17.8e} {dc[i]:17.8e} {f0[i]:17.8e} {tau[i]:17.8e} {sigma[i]:17.8e} {vel[i]:17.8e} {taudot[i]:17.8e} {sigmadot[i]:17.8e}\n')

