#parameter file generator for BP5
import numpy as np
import pandas as pd

# Constants
imax = 100
jmax = 40
ncell = imax * jmax
pi = 4.0 * np.arctan(1.0)
ds0 = 1.0
mu0 = 0.6
vref = 1e-6
vs = 3.464
velinit = 1e-9
rigid = 32.04
a0 = 0.004
b0 = 0.03
dc0 = 0.14
a_max = 0.04

# Arrays
rake=np.zeros(ncell);a=np.zeros(ncell);b=np.zeros(ncell);dc=np.zeros(ncell);f0=np.zeros(ncell)
tau=np.zeros(ncell);sigma=np.zeros(ncell);vel=np.zeros(ncell);taudot=np.zeros(ncell);sigmadot=np.zeros(ncell)
xcol=np.zeros(ncell);zcol=np.zeros(ncell)

k=-1
for i in range(imax):
    for j in range(jmax):
        k=k+1
        f0[k] = mu0
        xcol[k]=(i+0.5-imax/2)*ds0
        zcol[k]=-(j+0.5)*ds0
        dep = -zcol[k]

        if (dep > 4.0 and dep < 16.0 and abs(xcol[k]) < 30.0):
            a[k] = a0
        elif (dep < 2.0 or dep > 18.0 or abs(xcol[k]) > 32.0):
            a[k] = a_max
        else:
            r = max(abs(dep - 10.0) - 6.0, abs(xcol[k]) - 30.0) / 2.0
            a[k] = a0 + r * (a_max - a0)

        r = max(abs(dep - 10.0) - 6.0, abs(xcol[k]) - 30.0) / 2.0
        a[k] = min(a0 + r * (a_max - a0), a_max)
        a[k] = max(a[k], a0)

        b[k] = b0
        dc[k] = dc0

        if (abs(xcol[k] + 24.0) < 6.0 and abs(dep - 10.0) < 6.0):
            dc[k] = 0.13

        sigma[k] = 25.0
        vel[k] = velinit
        dep = -zcol[k]
        if (abs(xcol[k] + 24.0) < 6.0 and abs(dep - 10.0) < 6.0):
            vel[k] = 0.03
tau=sigma*a*np.arcsinh(0.5*vel/vref*np.exp((mu0+b0*np.log(vref/velinit))/a))+rigid/(2*vs)*vel

data=(rake,a,b,dc,f0,tau,sigma,vel,taudot,sigmadot)
data=np.transpose(data)
df=pd.DataFrame(data)
print(df)
df.to_csv('bp5rparam.dat', sep='\t', index=False, header=None)
