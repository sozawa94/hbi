import numpy as np
import matplotlib.pyplot as plt

jobid=0
coord=np.loadtxt('xyz'+str(jobid)+'_0.dat')
x=coord[:,0]
ncell=len(x)

fp = open('vel'+ str(jobid) +'_0.dat','rb')
ary = np.fromfile(fp, np.float64, -1)
fp.close()
nt=int(ary.shape[0]/ncell)
vel=ary.reshape(nt,ncell)

plt.plot(x,vel[10,:])
