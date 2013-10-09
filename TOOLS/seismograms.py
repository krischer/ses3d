# -*- coding: iso-8859-1 -*-

import numpy as np
import matplotlib.pylab as plt

#########################################################################
# define seismogram class
#########################################################################

class ses3d_seismogram:

  def __init__(self):

    self.nt=10
    self.dt=1
    self.t=np.linspace(0,10,10)

    self.rx=0.0
    self.ry=0.0
    self.rz=0.0

    self.sx=0.0
    self.sy=0.0
    self.sz=0.0

    self.trace_x=0.0*np.arange(10,dtype=float)
    self.trace_y=0.0*np.arange(10,dtype=float)
    self.trace_z=0.0*np.arange(10,dtype=float)

    self.integrate=False

  #########################################################################
  # read seismograms
  #########################################################################

  def read(self,directory,station_name,integrate=False):
    """ read seismogram
    read(directory,station_name,integrate)

    directory: directory where seismograms are located
    station_name: name of the station, without '.x', '.y' or '.z'
    integrate: integrate original velocity seismograms to displacement seismograms (important for adjoint sources)
    """

    self.integrate=integrate

    # open files ====================================================

    fx=open(directory+station_name+'.x','r')
    fy=open(directory+station_name+'.y','r')
    fz=open(directory+station_name+'.z','r')

    # read content ==================================================

    fx.readline()
    self.nt=int(fx.readline().strip().split('=')[1])
    self.dt=float(fx.readline().strip().split('=')[1])
    fx.readline()
    line=fx.readline().strip().split('=')
    self.rx=float(line[1].split('y')[0])
    self.ry=float(line[2].split('z')[0])
    self.rz=float(line[3])
    print 'receiver: colatitude={} deg, longitude={} deg, depth={} m'.format(self.rx,self.ry,self.rz)
    fx.readline()
    line=fx.readline().strip().split('=')
    self.sx=float(line[1].split('y')[0])
    self.sy=float(line[2].split('z')[0])
    self.sz=float(line[3])
    print 'source: colatitude={} deg, longitude={} deg, depth={} m'.format(self.sx,self.sy,self.sz)

    for k in range(7):
      fy.readline()
      fz.readline()

    self.trace_x=np.empty(self.nt,dtype=np.float64)
    self.trace_y=np.empty(self.nt,dtype=np.float64)
    self.trace_z=np.empty(self.nt,dtype=np.float64)

    for k in range(self.nt):
      self.trace_x[k]=float(fx.readline().strip())
      self.trace_y[k]=float(fy.readline().strip())
      self.trace_z[k]=float(fz.readline().strip())

    self.t=np.linspace(0,self.nt*self.dt,self.nt)

    # integrate to displacement seismograms =========================

    if integrate==True:
      self.trace_x=np.cumsum(self.trace_x)*self.dt
      self.trace_y=np.cumsum(self.trace_y)*self.dt
      self.trace_z=np.cumsum(self.trace_z)*self.dt

    # close files ===================================================

    fx.close()
    fy.close()
    fz.close()

  #########################################################################
  # plot seismograms
  #########################################################################

  def plot(self,scaling=1.0):
    """ plot seismograms
    plot(scaling)
    scaling: scaling factor for the seismograms
    """
    plt.subplot(311)
    plt.plot(self.t,scaling*self.trace_x,'k')
    plt.grid(True)
    plt.xlabel('time [s]')

    if self.integrate==True:
      plt.ylabel(str(scaling)+'*u_theta [m]')
    else:
      plt.ylabel(str(scaling)+'*v_theta [m/s]')

    plt.subplot(312)
    plt.plot(self.t,scaling*self.trace_y,'k')
    plt.grid(True)
    plt.xlabel('time [s]')

    if self.integrate==True:
      plt.ylabel(str(scaling)+'*u_phi [m]')
    else:
      plt.ylabel(str(scaling)+'*v_phi [m/s]')

    plt.subplot(313)
    plt.plot(self.t,scaling*self.trace_z,'k')
    plt.grid(True)
    plt.xlabel('time [s]')

    if self.integrate==True:
      plt.ylabel(str(scaling)+'*u_r [m]')
    else:
      plt.ylabel(str(scaling)+'*v_r [m/s]')

    plt.show()
