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
      
def animation_plot(jobid, field, xax='x', yax='z', nt=0, *,
                   output_file=None, tmin=None, tmax=None, vmin=-14, vmax=2,
                   fps=20, dpi=120, label=None, radius=None, show_vmax=False,
                   gif_file=None, preview_file=None):
    """Animate saved fields from the current data directory.

    Time limits are in years. Each frame is one saved snapshot: playback
    does not imply uniform physical-time spacing. Existing positional
    arguments remain supported; outputs default to the current directory.
    """
    from pathlib import Path

    time_data = np.atleast_2d(np.loadtxt(f'time{jobid}.dat'))
    coord = np.atleast_2d(np.loadtxt(f'xyz{jobid}.dat'))
    nc = len(coord)
    raw = np.fromfile(f'{field}{jobid}.dat', np.float64)
    if raw.size % nc or raw.size // nc != len(time_data):
        raise ValueError('Saved field frames and time rows must match exactly')
    data = raw.reshape(-1, nc)
    if not np.all(np.isfinite(data)) or not np.all(np.diff(time_data[:, 1]) > 0):
        raise ValueError('Nonfinite field or nonincreasing snapshot times')
    time = time_data[:, 1] / year
    selected = np.arange(len(time))
    if nt:
        selected = selected[:nt]
    if tmin is not None:
        selected = selected[time[selected] >= tmin]
    if tmax is not None:
        selected = selected[time[selected] <= tmax]
    if not len(selected):
        raise ValueError('No snapshots within the requested window')
    is_velocity = field in ['vel', 'veln']
    values = np.log10(np.maximum(np.abs(data), 1e-20)) if is_velocity else data
    switcher = {'x': 0, 'y': 1, 'z': 2}
    x, z = coord[:, switcher[xax]], coord[:, switcher[yax]]
    if show_vmax:
        if not is_velocity:
            raise ValueError('show_vmax requires a velocity field')
        fig = plt.figure(figsize=(8, 8))
        grid = fig.add_gridspec(2, 1, height_ratios=[3, 1])
        ax, trace = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
        monitor_path = Path(f'monitor{jobid}.dat')
        if field == 'vel' and monitor_path.exists():
            monitor = np.atleast_2d(np.loadtxt(monitor_path))
            trace.plot(monitor[:, 1] / year, monitor[:, 2], color='0.4', lw=1)
        else:
            trace.plot(time, values.max(axis=1), color='0.4', lw=1)
        cursor = trace.axvline(time[selected[0]], color='crimson', lw=1.5)
        point, = trace.plot([], [], 'o', color='crimson', ms=4)
        trace.set(xlim=(time[selected[0]], time[selected[-1]]), ylim=(vmin, vmax),
                  xlabel='Time (years)', ylabel='log10 Vmax')
        trace.grid(alpha=0.2)
    else:
        fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')
    sc = ax.scatter(x, z, c=values[selected[0]], vmin=vmin, vmax=vmax,
                    s=12, marker='s', linewidths=0, cmap='viridis', rasterized=True)
    ax.set(xlabel=f'{xax} (km)', ylabel=f'{yax} (km)')
    if radius is not None:
        angle = np.linspace(0, 2*np.pi, 300)
        ax.plot(radius*np.cos(angle), radius*np.sin(angle), 'k--', lw=0.8)
    fig.colorbar(sc, ax=ax, label='log10 |V| (m/s)' if is_velocity else field)
    fig.suptitle(label or f'{field}, run {jobid}', fontsize=13)
    fig.text(0.5, 0.015, 'One frame per saved snapshot; physical time spacing varies',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, 0.045, 1, 0.94))

    def animate(index):
        sc.set_array(values[index])
        ax.set_title(f'Time = {time[index]:.4f} yr; step {int(time_data[index, 0])}', fontsize=11)
        if show_vmax:
            cursor.set_xdata([time[index], time[index]])
            point.set_data([time[index]], [values[index].max()])
        return (sc,)

    target = Path(output_file or f'{field}{jobid}.mp4')
    target.parent.mkdir(parents=True, exist_ok=True)
    anim = animation.FuncAnimation(fig, animate, interval=1000/fps,
                                   frames=selected, blit=False)
    anim.save(str(target), writer=animation.FFMpegWriter(fps=fps, codec='libx264',
              extra_args=['-pix_fmt', 'yuv420p']), dpi=dpi)
    if gif_file:
        Path(gif_file).parent.mkdir(parents=True, exist_ok=True)
        anim.save(str(gif_file), writer=animation.PillowWriter(fps=fps), dpi=70)
    if preview_file:
        peak_index = selected[np.argmax(np.max(values[selected], axis=1))]
        animate(peak_index)
        fig.savefig(preview_file, dpi=dpi)
    plt.close(fig)
    return {'run': jobid, 'frames': len(selected), 'start_year': float(time[selected[0]]),
            'end_year': float(time[selected[-1]]), 'fps': fps,
            'snapshot_max_velocity': float(np.abs(data[selected]).max()) if is_velocity else None,
            'output': str(target)}
