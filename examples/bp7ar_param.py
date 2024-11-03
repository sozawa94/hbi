#parameter distribution
#!cd /work/hp220105o/i25004/hbi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ncell=6400

ds=0.01
a_max=0.016
a0=0.004
b0=0.010
dc0=0.53e-3
mu0=0.6
vref=1e-6
sigma0=25.0
rigid=32.04
cs=3.464
vel0=1e-9
tau0=sigma0*a_max*np.arcsinh(0.5*vel0/vref*np.exp((mu0+b0*np.log(vref/vel0))/a_max))#+rigid/(2*Vs)*vel(i)
eta=rigid/2/cs

a=np.zeros(ncell);tau=np.zeros(ncell)
x=np.zeros(ncell);z=np.zeros(ncell)
for i in range(80):
    for j in range(80):
        k=i+j*80
        x[k]=(i+0.5-40)*ds
        z[k]=(j+0.5-40)*ds
        r=np.sqrt(x[k]**2+z[k]**2)
        if r<0.2:
            a[k]=a0
        else:
            a[k]=a_max
        tau[k]=sigma0*a[k]*np.arcsinh(0.5*vel0/vref*np.exp((mu0+b0*np.log(vref/vel0))/a[k]))+eta*vel0
    #f.write(0.0),str(a[i]))#,b0,dc0,f0,tau0,sigma0,vel0,0.0,0.0)
plt.scatter(x,z,c=tau,s=10)
plt.colorbar()

data=(a,tau)
data=np.transpose(data)
df=pd.DataFrame(data, columns=['a','tau'])
#df=df.set_index()
print(df)
df.to_csv('../examples/bp7ar_param.dat', sep='\t', index=False)
