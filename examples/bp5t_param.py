#parameter distribution (for triangular mesh)
#!cd /work/hp220105o/i25004/hbi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ncell=9250
x=np.zeros(ncell);z=np.zeros(ncell)
with open('/work/hp220105o/i25004/hbi/examples/bp5t.stl') as f: #need to be changed
    lines = f.read()
    d=lines.split()
    print(d[12])
    for k in range(ncell):
        xs1=float(d[4+21*k+8]);xs2=float(d[4+21*k+12]);xs3=float(d[4+21*k+16])
        zs1=float(d[4+21*k+10]);zs2=float(d[4+21*k+14]);zs3=float(d[4+21*k+18])
        x[k]=(xs1+xs2+xs3)/3
        z[k]=(zs1+zs2+zs3)/3

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

for i in range(ncell):
    dep = -z[i]
    dc[i]=dc0
    r = max(abs(dep - 10) - 6, abs(x[i]) - 30) / 2
    a[i] = min(a0 + r * (a_max - a0), a_max)
    a[i] = max(a[i], a0)
    vel[i]=vel0
    if abs(x[i] + 24) < 6 and abs(dep - 10) < 6:
        vel[i] = 3e-2
    if abs(x[i] + 24) < 6 and abs(dep - 10) < 6:
        dc[i] = 0.13
    
    tau[i] = sigma0*(a[i]*np.arcsinh(0.5*vel[i]/vref*np.exp((mu0+b0*np.log(vref/vel0))/a[i]))+rigid/(2*cs)*vel[i])
    #f.write(0.0),str(a[i]))#,b0,dc0,f0,tau0,sigma0,vel0,0.0,0.0)
plt.scatter(x,z,c=vel,s=10)
plt.colorbar()

data=(a,dc,tau,vel)
data=np.transpose(data)
df=pd.DataFrame(data, columns=['a','dc','tau','vel'])
df.to_csv('../examples/bp5t_param.dat', sep='\t', index=False)
