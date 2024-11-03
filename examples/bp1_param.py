#parameter distribution for BP1QD
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

ncell=1600
ds=0.025
H=15.0
h=3.0
a_max=0.025
a0=0.010
a=np.zeros(ncell)
x=np.zeros(ncell)

for k in range(ncell):
    x[k]=(k+0.5)*ds
    if x[k]<H:
        a[k]=a0
    else:
        a[k]=min(a_max,a0+(a_max-a0)*(x[k]-H)/h)

plt.plot(x,a)
#plt.xlim([15,20])

data=(a)
data=np.transpose(data)
df=pd.DataFrame(data, columns=['a'])
print(df)
df.to_csv('../examples/bp1_param.dat', sep='\t', index=False)
