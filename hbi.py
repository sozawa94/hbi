import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams["font.size"] = 14  
year=365*24*3600
day=24*3600

def load(jobid,field):
    coord=np.loadtxt(f'xyz{jobid}.dat')
    ncell=len(coord)
    fp = open(f'{field}{jobid}.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    d=np.loadtxt(f'time{jobid}.dat')
    time=d[:,1]
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
    return time,d,coord
#fault geoemtry in 2D
def geomplot(jobid,alpha=1,size=10):
    coord=np.loadtxt(f'xyz{jobid}.dat')
    x=coord[:,0]
    y=coord[:,1]
    plt.scatter(x,y,alpha=alpha,s=size)
    
#fault geoemtry in 3D    
def geomplot3d(jobid):
    coord=np.loadtxt(f'xyz{jobid}.dat')
    x=coord[:,0]
    y=coord[:,1]
    z=coord[:,2]
    plt.scatter(x,y,z,s=3)
    
#space time step plot in 2D    
def ssplot(jobid,field,vmin=None,vmax=None):
    coord=np.loadtxt(f'xyz{jobid}.dat')
    dep=coord[:,0]
    ncell=len(dep)
    fp = open(f'{field}{jobid}.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
    if field in ['vel', 'veln']:
        d=np.log10(abs(d))
    x=np.arange(0,n,1)
    
    fig=plt.figure(figsize=(12,4))
    plt.pcolormesh(x,dep,d,cmap=plt.cm.viridis,vmin=vmin,vmax=vmax)
    
    switcher = {
        "vel": "log10 Slip rate (m/s)",
        "slip": "Slip (m)",
        "tau": "Shear stress (MPa)",
        "sigma": "Effective normal stress (MPa)",
        "pf": "Fluid pressure (MPa)",
    }
    plt.colorbar(label=switcher.get(field))
    plt.xlabel('Time step')
    plt.ylabel('x [km]')
    #plt.title('jobid='+str(job))
    #plt.xlim([xmin,xmax])
    #plt.ylim([zmax,0])
    
#space time plot in 2D    
def stplot(jobid, field, vmin=None, vmax=None):
    coord=np.loadtxt(f'xyz{jobid}.dat')
    dep=coord[:,0]
    ncell=len(dep)
    fp = open(field+ str(jobid) +'.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    d=np.loadtxt(f'time{jobid}.dat')
    time=d[:,1]
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
    if field=='vel':
        d=np.log10(d)
    
    fig=plt.figure(figsize=(12,4))

    plt.pcolormesh(time/year,dep,d,cmap=plt.cm.viridis,vmin=vmin,vmax=vmax)
    switcher = {
        "vel": "log10 Slip rate (m/s)",
        "slip": "Slip (m)",
        "tau": "Shear stress (MPa)",
        "sigma": "Effective normal stress (MPa)",
        "pf": "Fluid pressure (MPa)",
    }
    plt.colorbar(label=switcher.get(field))
    plt.xlabel('Time (yr)')
    plt.ylabel('Location (km)')
    #plt.title('jobid='+str(job))
    #plt.xlim([xmin,xmax])
    #plt.ylim([np.max(dep),0])

#for frictional-viscous model
def stvelplot(jobid):
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
    dep=coord[:,1]
    ncell=len(dep)
    fp = open('vel'+ str(jobid) +'.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    fp = open('vc'+ str(jobid) +'.dat','rb')
    ary2 = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    d=np.loadtxt('time'+str(jobid)+'.dat')
    time=d[:,1]
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
    ary2=ary2[:n*ncell]
    d2=ary2.reshape(n,ncell)
    d2=d2.T
    d3=d+d2
    d=np.log10(d3)
    fig=plt.figure(figsize=(12,4))

    plt.pcolormesh(time/year,dep,d,cmap=plt.cm.bwr,vmin=-12,vmax=-6)
    plt.axhline(20)

    plt.colorbar(label='log10 slip rate [m/s]')
    plt.xlabel('Time [yr]')
    plt.ylabel('Depth [km]')
    #plt.title('jobid='+str(job))
    #plt.xlim([xmin,xmax])
    plt.ylim([np.max(dep),0])
    plt.savefig('/home/i25004/.notebook/Visco/'+str(jobid)+'.png',dpi=300,bbox_inches='tight')

#cumulative slip distribution in 2D
def cumslip(jobid):
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
    loc=coord[:,0]
    ncell=len(loc)
    fp = open('slip'+ str(jobid) +'.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
   
    fig=plt.figure(figsize=(6,4))
    plt.plot(d,loc,c='b')
    plt.xlabel('Slip [m]')
    plt.ylabel('x [km]')

#snapshot at a given time step in 2D     
def snap(jobid,field,ts=0,label="none",xindex=False):
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
    loc=coord[:,0]
    ncell=len(loc)
    fp = open(field+ str(jobid) +'.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
    if field=='vel':
        d=np.log10(abs(d))
   
    #fig=plt.figure(figsize=(6,4))
    if xindex:
        plt.plot(d[:,ts],label=label)
    else:
        plt.plot(loc,d[:,ts],label=label)
    #plt.xlabel('Time step')
    #plt.ylabel('x [km]')
    
def sigmasnap(jobid,ts=0):
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
    loc=coord[:,0]
    ncell=len(loc)
    fp = open('pf'+ str(jobid) +'.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    fp = open('sigma'+ str(jobid) +'.dat','rb')
    ary2 = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
    
    ary2=ary2[:n*ncell]
    d2=ary2.reshape(n,ncell)
    d2=d2.T

    #fig=plt.figure(figsize=(6,4))
    plt.plot(loc,d[:,ts]+d2[:,ts])
    plt.xlabel('Time step')
    plt.ylabel('x [km]')

#initial condition in 2D
def initplot(field,jobid):
    coord=np.loadtxt(f'xyz{jobid}.dat')
    loc=coord[:,1]
    ncell=len(loc)
    fp = open(field+ str(jobid) +'.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T   
    fig=plt.figure(figsize=(6,4))
    plt.plot(loc,d)
    plt.xlabel('Time step')
    plt.ylabel('x [km]')

#plot the slip distribution for each earthquake in 2D
def EQslip(jobid,xindex=False):
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
    loc=coord[:,0]
    ncell=len(loc)
    fp = open('EQslip'+ str(jobid) +'.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
   
    fig=plt.figure(figsize=(6,4))
    if xindex:
        plt.plot(d[:,:])
    else:
        plt.plot(loc,d[:,:])
    plt.xlabel('x (km)')
    plt.ylabel('Slip (m)')

#plot the time series from monitorX.dat
def monitorplot(jobid, ycol=2,label='none',time='year',offset=0,alpha=1):
    dd=np.loadtxt(f'monitor{jobid}.dat')
    
    if time=="year":
        plt.plot((dd[:,1]-offset)/year,dd[:,ycol],label=label,alpha=alpha)
    if time=="day":
        plt.plot((dd[:,1]-offset)/day,dd[:,ycol],label=label,alpha=alpha)
    if time=="sec":
        plt.plot(dd[:,1]-offset,dd[:,ycol],label=label,alpha=alpha)
    plt.xlabel('time (yr)')
    switcher = {
        2: "log10 max slip rate (m/s)",
        3: "average slip (m)",
        4: "average friction",
        5: "max sigma_eff (MPa)",
        6: "min sigma_eff (MPa)",
    }
    plt.ylabel(switcher.get(ycol))
    
#same with monitorX.dat but x axis is time step   
def monitortsplot(jobid, ycol=2):
    dd=np.loadtxt(f'monitor{jobid}.dat')
    plt.plot(dd[:,0],dd[:,ycol])
    plt.xlabel('time step')
    switcher = {
        2: "log10 max slip rate (m/s)",
        3: "average slip (m)",
        4: "average friction",
        5: "max sigma_eff (MPa)",
        6: "min sigma_eff (MPa)",
    }
    plt.ylabel(switcher.get(ycol))

#plot time series for local data
def localplot(jobid,loc,ycol,skiprows=0):
    dd=np.loadtxt(f'local{jobid}-{loc}.dat',skiprows=skiprows)
    plt.plot(dd[:,1]/year,dd[:,ycol])
    plt.xlabel('time (yr)')

#3D snapshot (only for single MPI)
def snap3drec(jobid,field,imax, jmax, ts=0, xax='x', yax='z',size=8):
    plt.figure(figsize=(size,size))
    plt.axes().set_aspect('equal')
    ds=0.01
    Lx=ds*imax
    nc=imax*jmax
    fp = open(field+ str(jobid) +'.dat','rb')
    ary2 = np.fromfile(fp, np.float64, -1)
    fp.close()
    nc2=len(ary2)
    nt=int(nc2/nc)
    y = np.linspace(-Lx/2, Lx/2, imax)

    dd=ary2[ts*nc:(ts+1)*nc]
    dd=dd.reshape(imax,jmax)
    dd=dd.transpose()
    if field=='vel':
        dd=np.log10(dd)
    plt.imshow(dd,cmap=plt.cm.viridis)
    plt.colorbar(label=field)

#3D snapshot (both for rectangular and triangular)
def snap3d(jobid,field,ts=0, xax='x', yax='y',size=8, vmin=None, vmax=None, time=True,cmap='viridis',absolute=False):
    plt.figure(figsize=(size,size))
    plt.axes().set_aspect('equal')
    ary1=np.loadtxt(f'xyz{jobid}.dat')
    nc=len(ary1)
    #print(nc)
    switcher = {"x": 0,"y": 1,"z": 2}
    x=ary1[:,switcher.get(xax)]
    y=ary1[:,switcher.get(yax)]
    fp = open(f'{field}{jobid}.dat','rb')
    d = np.fromfile(fp, np.float64, -1)
    fp.close()
    nt=int(len(d)/nc)
    print(nt)
    d=d[ts*nc:(ts+1)*nc]
    if absolute:
        d=abs(d)
    if field=='vel' or field=='veln':
        d=np.log10(abs(d))
    #plt.tripcolor(x, z, d, cmap="viridis")
    plt.scatter(x, y, c=d, s=1000/np.sqrt(nc), cmap=cmap,vmin=vmin, vmax=vmax)
    plt.xlabel(xax + ' (km)')
    plt.ylabel(yax + ' (km)')
    if time:
        d=np.loadtxt(f'time{jobid}.dat')
        time=d[:,1]/year
        plt.title('Time = %f (yr)' %time[ts])
    plt.colorbar(label=field)
    
def locking3d(jobid,ts=0, xax='x', yax='y',cmap='viridis',vpl=1e-9,size=8, time=True,):
    plt.figure(figsize=(size,size))
    plt.axes().set_aspect('equal')
    ary1=np.loadtxt(f'xyz{jobid}.dat')
    nc=len(ary1)
    #print(nc)
    switcher = {"x": 0,"y": 1,"z": 2}
    x=ary1[:,switcher.get(xax)]
    y=ary1[:,switcher.get(yax)]
    fp = open(f'vel{jobid}.dat','rb')
    d = np.fromfile(fp, np.float64, -1)
    fp.close()
    nt=int(len(d)/nc)
    print(nt)
    d=d[ts*nc:(ts+1)*nc]
    #plt.tripcolor(x, z, d, cmap="viridis")
    plt.scatter(x, y, c=(vpl-d)/vpl, s=1000/np.sqrt(nc), cmap=cmap,vmin=-1, vmax=1)
    plt.xlabel(xax + ' (km)')
    plt.ylabel(yax + ' (km)')
    if time:
        d=np.loadtxt(f'time{jobid}.dat')
        time=d[:,1]/year
        plt.title('Time = %f (yr)' %time[ts])
    plt.colorbar(label='locking')
#plot the slip distribution for each earthquake in 2D
def EQslip3d(jobid, nstart=0, ne = 0, xax='x', yax='y', magmin=0, save=False, size=6): 
    ary1=np.loadtxt(f'xyz{jobid}.dat')
    nc=len(ary1)
    switcher = {"x": 0,"y": 1,"z": 2}
    x=ary1[:,switcher.get(xax)]
    y=ary1[:,switcher.get(yax)]
    fp = open(f'EQslip{jobid}.dat','rb')
    ary2 = np.fromfile(fp, np.float64, -1)
    fp.close()

    d=np.loadtxt(f'event{jobid}.dat',skiprows=0)
    hypo=d[:,4]
    net=len(d)
    if ne == 0:
        ne=len(d)

    hypox=np.zeros(net)
    hypoy=np.zeros(net)
    hypoz=np.zeros(net)
    time=d[:,2]/year
    mag=d[:,3]

    for i in range(net):
        n=int(hypo[i])
    #print(n)
        hypox[i]=ary1[n-1,switcher.get(xax)]
        hypoy[i]=ary1[n-1,switcher.get(yax)]

    for k in range(nstart,nstart+ne):
        if d[k,3]>magmin:
            dd=ary2[k*nc:(k+1)*nc]
            # Sort by values so that smaller values are plotted first
            sorted_indices = np.argsort(dd)
            xs = x[sorted_indices]
            ys = y[sorted_indices]
            dds = dd[sorted_indices]
            fig = plt.figure(figsize=(size,size))
            plt.axes().set_aspect('equal')
            plt.scatter(xs,ys,c=dds,s=1,cmap='Blues',vmin=0,vmax=np.max(dds))
            plt.colorbar(label='Slip (m)', shrink=0.6)
            plt.scatter(hypox[k], hypoy[k],c='r')
            plt.title(f'Event ID = {k+1}, Time = {time[k]:.3f} (yr), Mw = {mag[k]:.2f}')
            if save:
                plt.savefig(f'/home/i25004/.notebook/HBI/EQslip{jobid}-{k+1}.png',dpi=300)

#moment duration scaling for a given velocity threshold
def momentvsduration(jobid,velth):
    d=np.loadtxt('monitor'+str(jobid)+'.dat')
    switch=0
    moment=[]
    vpeak=[]
    duration=[]

    for i in range(len(d)):
        if switch ==0 and d[i,2] > velth:
            switch = 1
            ionset=i
            tonset=d[i,1]
            momonset = d[i,3]

        if switch == 1 and d[i,2] < velth:
            switch = 0
            tend = d[i,1]
            momend = d[i,3]
            vpeak = np.append(vpeak, np.max(d[ionset:i,2]))
            duration = np.append(duration, tend-tonset)
            moment = np.append(moment, momend-momonset)
    plt.scatter(moment, duration,c=vpeak)
    plt.colorbar(label="Peak slip rate")
    plt.yscale('log')
    plt.xlabel('Moment')
    plt.ylabel('Duration [s]')
      
def animation_plot(jobid, field, xax='x', yax='z', nt=0):
    d = np.loadtxt(f'time{jobid}.dat')
    time = d[:, 1] / year

    ary1 = np.loadtxt(f'xyz{jobid}.dat')
    nc = len(ary1)
    print("Number of coordinates:", nc)

    switcher = {"x": 0, "y": 1, "z": 2}
    x = ary1[:, switcher.get(xax)]
    y = ary1[:, switcher.get(yax)]

    # Read binary field data
    with open(f'{field}{jobid}.dat', 'rb') as fp:
        data = np.fromfile(fp, np.float64)

    if nt == 0:
        nt = int(len(data) / nc)

    if field in ['vel', 'veln']:
        data = np.log10(np.abs(data) + 1e-20)  # prevent log10(0)

    # Set up figure and axis once
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')
    sc = ax.scatter(x, y, c=data[:nc], vmin=-14, vmax=2, s=3, cmap='viridis')
    plt.colorbar(sc, ax=ax, label="log10 V (m/s)")

    def animate(i):
        dd = data[i * nc:(i + 1) * nc]
        sc.set_array(dd)
        ax.set_title('Time = %.2f yr' % time[i])
        return sc,

    anim = animation.FuncAnimation(fig, animate, interval=100, frames=nt, blit=False)

    writervideo = animation.FFMpegWriter(fps=20)
    anim.save(f'/home/i25004/.notebook/HBI/{field}{jobid}.mp4', writer=writervideo, dpi=300)
