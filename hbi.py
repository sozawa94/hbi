import numpy as np
import matplotlib.pyplot as plt
yr=365*24*3600

#fault geoemtry in 2D
def geomplot(jobid):
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
    x=coord[:,0]
    y=coord[:,1]
    plt.plot(x,y)
    
#fault geoemtry in 3D    
def geomplot3d(jobid):
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
    x=coord[:,0]
    y=coord[:,1]
    z=coord[:,2]
    plt.scatter(x,y,z,s=3)
    
#space time step plot in 2D    
def ssplot(jobid,field):
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
    dep=coord[:,0]
    ncell=len(dep)
    fp = open(field+ str(jobid) +'.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
    if field=='vel':
        d=np.log10(d)
    x=np.arange(0,n,1)
    
    fig=plt.figure(figsize=(12,4))
    if field=='vel':
        plt.pcolormesh(x,dep,d,cmap=plt.cm.viridis,vmin=-12,vmax=0)
    else:
        plt.pcolormesh(x,dep,d,cmap=plt.cm.viridis)
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
def stplot(jobid,field):
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
    dep=coord[:,0]
    ncell=len(dep)
    fp = open(field+ str(jobid) +'.dat','rb')
    ary = np.fromfile(fp, np.float64, -1)
    fp.close()
    
    d=np.loadtxt('time'+str(jobid)+'.dat')
    time=d[:,1]
    
    n=int(ary.shape[0]/ncell)
    print(n)
    ary=ary[:n*ncell]
    d=ary.reshape(n,ncell)
    d=d.T
    if field=='vel':
        d=np.log10(d)
    fig=plt.figure(figsize=(12,4))
    if field=='vel':
        plt.pcolormesh(time/yr,dep,d,cmap=plt.cm.bwr,vmin=-12,vmax=-6)
    else:
        plt.pcolormesh(time/yr,dep,d,cmap=plt.cm.viridis)
    switcher = {
        "vel": "log10 Slip rate (m/s)",
        "slip": "Slip (m)",
        "tau": "Shear stress (MPa)",
        "sigma": "Effective normal stress (MPa)",
        "pf": "Fluid pressure (MPa)",
    }
    plt.colorbar(label=switcher.get(field))
    plt.colorbar()
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

    plt.pcolormesh(time/yr,dep,d,cmap=plt.cm.bwr,vmin=-12,vmax=-6)
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
def snap(jobid,field,ts=0):
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
        d=np.log10(d)
   
    #fig=plt.figure(figsize=(6,4))
    plt.plot(loc,d[:,ts])
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
    coord=np.loadtxt('xyz'+str(jobid)+'.dat')
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
def EQslip(jobid):
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
    plt.plot(loc,d[:,:])
    plt.xlabel('x (km)')
    plt.ylabel('Slip (m)')

#plot the time series from monitorX.dat
def monitorplot(jobid, ycol=2,label='none'):
    dd=np.loadtxt('monitor'+str(jobid)+'.dat')
    plt.plot(dd[:,1]/yr,dd[:,ycol],label=label)
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
    dd=np.loadtxt('monitor'+str(jobid)+'.dat')
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
def localplot(jobid,loc,col):
    dd=np.loadtxt('local'+str(jobid)+'-'+str(loc)+'.dat')
    plt.plot(dd[:,1]/yr,dd[:,col])
    plt.xlabel('time (yr)')

#3D snapshot (only for single MPI)
def snap3drec(jobid,field,imax, jmax, ts=0, xax='x', yax='z'):
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
def snap3d(jobid,field,ts=0, xax='x', yax='y',size=8):
    plt.figure(figsize=(size,size))
    plt.axes().set_aspect('equal')
    ary1=np.loadtxt('xyz'+ str(jobid) +'.dat')
    nc=len(ary1)
    #print(nc)
    x=ary1[:,0]
    y=ary1[:,1]
    z=ary1[:,2]
    fp = open(field+ str(jobid) +'.dat','rb')
    d = np.fromfile(fp, np.float64, -1)
    fp.close()
    nt=int(len(d)/nc)
    print(nt)
    d=d[ts*nc:(ts+1)*nc]
    if field=='vel':
        d=np.log10(d)
    #plt.tripcolor(x, z, d, cmap="viridis")
    plt.scatter(x, y, c=d, s=1000/np.sqrt(nc), cmap="viridis")
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.colorbar(label=field)
    
#plot the slip distribution for each earthquake in 2D
def EQslip3d(jobid): 
    ary1=np.loadtxt('xyz'+ str(jobid) +'.dat')
    nc=len(ary1)
    fp = open('EQslip'+ str(jobid) +'.dat','rb')
    ary2 = np.fromfile(fp, np.float64, -1)
    fp.close()

    d=np.loadtxt('event'+ str(jobid) +'.dat',skiprows=0)
    hypo=d[:,4]
    ne=len(d)

    hypox=np.zeros(ne)
    hypoy=np.zeros(ne)
    hypoz=np.zeros(ne)


    for i in range(ne):
        n=int(hypo[i])
    #print(n)
        hypox[i]=ary1[n-1,0]
        hypoy[i]=ary1[n-1,1]
        hypoz[i]=ary1[n-1,2]

    for k in range(ne):
        dd=ary2[k*nc:(k+1)*nc]
        fig = plt.figure(figsize=(6,3))
        plt.scatter(ary1[:,0],ary1[:,1],c=dd,s=5,cmap='Blues')#,vmin=-12,vmax=0)
        plt.scatter(hypox[k], hypoy[k],c='r')

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
