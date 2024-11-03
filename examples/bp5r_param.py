#!cd /work/hp220105o/i25004/hbi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

imax = 100
jmax = 40
ncell = imax * jmax
ds0=1.0
x=np.zeros(ncell);z=np.zeros(ncell)

# Frictional parameters
a0 = 0.004
b0 = 0.03
dc0 = 0.14
a_max = 0.04
mu0=0.6
vref=1e-6
sigma0=25.0
rigid=32.04
cs=3.464
vel0=1e-9
tau0=sigma0*a_max*np.arcsinh(0.5*vel0/vref*np.exp((mu0+b0*np.log(vref/vel0))/a_max))#+rigid/(2*Vs)*vel(i)
eta=rigid/2/cs

a=np.zeros(ncell);dc=np.zeros(ncell)
tau=np.zeros(ncell);vel=np.zeros(ncell)

k=-1
for i in range(imax):
    for j in range(jmax):
        k=k+1
        f0[k] = mu0
        x[k]=(i+0.5-imax/2)*ds0
        z[k]=-(j+0.5)*ds0
        dep = -z[k]
        dc[k]=dc0
        r = max(abs(dep - 10) - 6, abs(x[k]) - 30) / 2
        a[k] = min(a0 + r * (a_max - a0), a_max)
        a[k] = max(a[k], a0)
        vel[k]=vel0
        if abs(x[k] + 24) < 6 and abs(dep - 10) < 6:
            vel[k] = 3e-2
        if abs(x[k] + 24) < 6 and abs(dep - 10) < 6:
            dc[k] = 0.13

        tau[k] = sigma0*(a[k]*np.arcsinh(0.5*vel[k]/vref*np.exp((mu0+b0*np.log(vref/vel0))/a[k]))+rigid/(2*cs)*vel[k])
        #f.write(0.0),str(a[i]))#,b0,dc0,f0,tau0,sigma0,vel0,0.0,0.0)
plt.scatter(x,z,c=tau,s=10)
plt.colorbar()

data=(a,dc,tau,vel)
data=np.transpose(data)
df=pd.DataFrame(data, columns=['a','dc','tau','vel'])
df.to_csv('../examples/bp5r_param.dat', sep='\t', index=False)
