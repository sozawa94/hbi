#parameter distribution for BP1QD
import numpy as np
import pandas as pd

ncell=1600

ds=0.025
H=15.0
h=3.0
a_max=0.025
a0=0.010
b0=0.015
dc0=0.008
mu0=0.6
vref=1e-6
sigma0=50.0
rigid=32.04
cs=3.464
vel0=1.0e-9
tau0=sigma0*a_max*np.arcsinh(0.5*vel0/vref*np.exp((mu0+b0*np.log(vref/vel0))/a_max))#+rigid/(2*Vs)*vel(i)
eta=rigid/2/cs
ncell=1600
print(rigid*dc0/b0/sigma0)
rake=np.zeros(ncell);a=np.zeros(ncell);b=np.zeros(ncell);dc=np.zeros(ncell);f0=np.zeros(ncell)
tau=np.zeros(ncell);sigma=np.zeros(ncell);vel=np.zeros(ncell);taudot=np.zeros(ncell);sigmadot=np.zeros(ncell)
x=np.zeros(ncell);z=np.zeros(ncell)

for k in range(ncell):
    b[k]=b0
    dc[k]=dc0
    f0[k]=mu0
    sigma[k]=sigma0
    vel[k]=vel0
    x[k]=(k+0.5)*ds
    if x[k]<H:
        a[k]=a0
    else:
        a[k]=min(a_max,a0+(a_max-a0)*(x[k]-H)/h)
        
    tau[k]=sigma[k]*a_max*np.arcsinh(0.5*vel0/vref*np.exp((mu0+b0*np.log(vref/vel0))/a_max))+eta*vel0
    #f.write(0.0),str(a[i]))#,b0,dc0,f0,tau0,sigma0,vel0,0.0,0.0)
#plt.plot(x,a-b)
#plt.xlim([15,20])

data=(rake,a,b,dc,f0,tau,sigma,vel,taudot,sigmadot)
data=np.transpose(data)
df=pd.DataFrame(data)
#df=df.set_index()
print(df)
df.to_csv('../examples/bp1_param.dat', sep='\t', index=False, header=None)
