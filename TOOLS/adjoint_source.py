# -*- coding: iso-8859-1 -*-

import numpy as np
import matplotlib.pylab as plt
import seismograms as seismo
import prolate_spheroidal as ps

#########################################################################
# define adjoint source class
#########################################################################

class adjoint_source:
  """ class for computing, visualising and storing adjoint sources

  The prodedure is as follows: (1) fetch a seismogram using fetch_seismogram, (2) compute an adjoint source using one of 
  the functions make_adsrc_*, (3) plot the adjoint source with plot, and (4) write it to a file with write.
  """

  def __init__(self):

    self.s=seismo.ses3d_seismogram()

    self.adsrc_x=np.zeros(10,dtype=float)
    self.adsrc_y=np.zeros(10,dtype=float)
    self.adsrc_z=np.zeros(10,dtype=float)

  #########################################################################
  # fetch seismogram
  #########################################################################

  def fetch_seismogram(self,seismogram):
    """ function for fetching a seismogram, i.e. a member of the ses3d_seismogram class
  
    fetch_seismogram(ses3d_seismogram)
    """

    self.s=seismogram

    self.adsrc_x=np.zeros(seismogram.nt,dtype=float)
    self.adsrc_y=np.zeros(seismogram.nt,dtype=float)
    self.adsrc_z=np.zeros(seismogram.nt,dtype=float)

  #########################################################################
  # plot adjoint source
  #########################################################################

  def plot(self):
    """ plot the adjoint source that has been computed previously
    """

    plt.subplot(311)
    plt.plot(self.s.t,self.adsrc_x,'k')
    plt.grid(True)
    plt.xlabel('time [s]')
    plt.ylabel('adsrc_theta')
    plt.title('adjoint source time function (in forward time)')

    plt.subplot(312)
    plt.plot(self.s.t,self.adsrc_y,'k')
    plt.grid(True)
    plt.xlabel('time [s]')
    plt.ylabel('adsrc_phi')

    plt.subplot(313)
    plt.plot(self.s.t,self.adsrc_z,'k')
    plt.grid(True)
    plt.xlabel('time [s]')
    plt.ylabel('adsrc_r')

    plt.show()

  #########################################################################
  # write adjoint source
  #########################################################################

  def write(self,directory,filename):
    """ write a pre-computed adjoint source to a file

    write(directory,filename)
    """

    fid=open(directory+filename,'w')

    fid.write('-- adjoint source ------------------\n')
    fid.write('-- source coordinates (colat,lon,depth)\n')
    fid.write(str(self.s.rx)+' '+str(self.s.ry)+' '+str(self.s.rz)+'\n')
    fid.write('-- source time function (x, y, z) --\n')

    for k in np.arange(self.s.nt-1,-1,-1):

      fid.write(str(self.adsrc_x[k])+' '+str(self.adsrc_y[k])+' '+str(self.adsrc_z[k])+'\n')

    fid.close()

  #########################################################################
  # make adjoint source for multi-taper traveltime measurement
  #########################################################################

  def make_adsrc_mttime(self,component,win_left,win_right,nu,verbose=False):
    """ compute adjoint source for a multi-taper measurement of a time delay at a specific frequency

    make_adsrc_mttime(component,win_left,win_right,nu,verbose=False)

    component: seismogram component ('x','y' or 'z') for (theta, phi or r)
    win_left: left window boundary [s]
    win_right: right window boundary [s]
    nu: frequency [Hz]
    verbose: additional information as screen output if 'True'
    """

    #- transcribe the trace to be considered

    if component=='x':
      d=self.s.trace_x
    elif component=='y':
      d=self.s.trace_y
    elif component=='z':
      d=self.s.trace_z
    else:
      print('enter correct component (x, y or z)')
      d=np.zeros(self.s.nt,dtype=float)

    #- compute start and stop indices (make sure N is an odd number)

    N_0=int(round(win_left/self.s.dt))
    N=int(2*round((round(win_right/self.s.dt)-N_0)/2.0)+1)

    if verbose==True:
      print 'starting index: '+str(N_0)
      print 'length of index series: '+str(N)

    #- set some parameters
    nu=nu*self.s.dt	#- non-dimensional normalised frequency
    k_max=4		#- maximum index of PSWF (number of tapers = k_max+1)
    W=3.0*k_max/(2.0*N)	#- frequency bandwidth

    # REMARK: Choosing the frequency bandwidth sufficiently large ensures a good time-localisation of the tapers. However, it also increases spectral leakage.

    if verbose==True:
      print 'frequency bandwidth of the tapers: '+str(W/self.s.dt)+' Hz'
      print 'non-dimensional frequency bandwidth of the tapers: '+str(W)
      print 'maximum number of PSWFs: '+str(k_max+1)

    #- compute the first k_max prolate spheroidal wave functions and sequences

    eigval,Psi=ps.PSFreq(N,W,k_max,np.arange(0,1))
    eigval,psi=ps.PSIdx(N,W,k_max,np.arange(-(N-1)/2,(N-1)/2+1))

    #- plot the tapers when verbose mode

    if verbose==True:
      for k in np.arange(k_max+1):
	taper=np.zeros(self.s.nt,dtype=float)
	taper[np.arange(N_0,N_0+N)]=np.real(psi[:,k])
	plt.plot(self.s.t,taper,color=[float(k)/float(2*k_max), float(k)/float(2*k_max), float(k)/float(2*k_max)])
	plt.grid(True)
	plt.xlabel('time [s]')
	plt.title('tapers')

      plt.plot((win_left, win_left),(-0.02, 0.02),'r')
      plt.plot((win_right, win_right),(-0.02, 0.02),'r')
      plt.show()

    #- compute the expansion coefficients eta^(k) and the spectral estimate U(nu)

    eta=np.zeros(k_max+1,dtype=complex)
    n=np.arange(N_0,N_0+N)
    u=0.0+0.0j

    for k in np.arange(k_max+1):
      eta[k]=np.sum(d[n]*psi[:,k]*np.exp(-2*np.pi*1j*n*nu)) 

    for k in np.arange(k_max+1):
      u=u+eta[k]*Psi[0,k]/eigval[k]

    if verbose==True:
      print 'expansion coeffcients:'+str(eta)
      print 'eigenvalues:'+str(eigval)
      print 'spectral estimate:'+str(u)

    #- compute adjoint source time function (with forward directed time axis)

    adsrc=np.zeros(self.s.nt,dtype=complex)

    for k in np.arange(k_max+1):
      adsrc[n]=adsrc[n]+Psi[0,k]*psi[:,k]*np.exp(-2*np.pi*1j*n*nu)/eigval[k]
   
    adsrc=adsrc/(2*np.pi*self.s.dt*nu*u)

    #- differentiate adjoint source time function for measurements on velocity seismograms
    #- The previous development assumes displacement seismograms

    if self.s.integrate==False:
      adsrc[N_0:N_0+N-1]=-np.diff(adsrc[N_0:N_0+N])/self.s.dt
      adsrc[N_0+N]=0.0

    #- transcribe to adjoint source components

    if component=='x':
      self.adsrc_x=np.imag(adsrc)
    elif component=='y':
      self.adsrc_y=np.imag(adsrc)
    elif component=='z':
      self.adsrc_z=np.imag(adsrc)


  #########################################################################
  # make adjoint source for multi-taper amplitude measurement
  #########################################################################

  def make_adsrc_mtamp(self,component,win_left,win_right,nu,verbose=False):
    """ compute adjoint source for a multi-taper measurement of a relative amplitude difference at a specific frequency

    make_adsrc_mtamp(component,win_left,win_right,nu,verbose=False)

    component: seismogram component ('x','y' or 'z') for (theta, phi or r)
    win_left: left window boundary [s]
    win_right: right window boundary [s]
    nu: frequency [Hz]
    verbose: additional information as screen output if 'True'
    """

    #- transcribe the trace to be considered

    if component=='x':
      d=self.s.trace_x
    elif component=='y':
      d=self.s.trace_y
    elif component=='z':
      d=self.s.trace_z
    else:
      print('enter correct component (x, y or z)')
      d=np.zeros(self.s.nt,dtype=float)

    #- compute start and stop indices (make sure N is an odd number)

    N_0=int(round(win_left/self.s.dt))
    N=int(2*round((round(win_right/self.s.dt)-N_0)/2.0)+1)

    if verbose==True:
      print 'starting index: '+str(N_0)
      print 'length of index series: '+str(N)

    #- set some parameters
    nu=nu*self.s.dt	#- non-dimensional normalised frequency
    k_max=4		#- maximum index of PSWF (number of tapers = k_max+1)
    W=3.0*k_max/(2.0*N)	#- frequency bandwidth

    # REMARK: Choosing the frequency bandwidth sufficiently large ensures a good time-localisation of the tapers. However, it also increases spectral leakage.

    if verbose==True:
      print 'frequency bandwidth of the tapers: '+str(W/self.s.dt)+' Hz'
      print 'non-dimensional frequency bandwidth of the tapers: '+str(W)
      print 'maximum number of PSWFs: '+str(k_max+1)

    #- compute the first k_max prolate spheroidal wave functions and sequences

    eigval,Psi=ps.PSFreq(N,W,k_max,np.arange(0,1))
    eigval,psi=ps.PSIdx(N,W,k_max,np.arange(-(N-1)/2,(N-1)/2+1))

    #- plot the tapers when verbose mode

    if verbose==True:
      for k in np.arange(k_max+1):
	taper=np.zeros(self.s.nt,dtype=float)
	taper[np.arange(N_0,N_0+N)]=np.real(psi[:,k])
	plt.plot(self.s.t,taper,color=[float(k)/float(2*k_max), float(k)/float(2*k_max), float(k)/float(2*k_max)])
	plt.grid(True)
	plt.xlabel('time [s]')
	plt.title('tapers')

      plt.plot((win_left, win_left),(-0.02, 0.02),'r')
      plt.plot((win_right, win_right),(-0.02, 0.02),'r')
      plt.show()

    #- compute the expansion coefficients eta^(k) and the spectral estimate U(nu)

    eta=np.zeros(k_max+1,dtype=complex)
    n=np.arange(N_0,N_0+N)
    u=0.0+0.0j

    for k in np.arange(k_max+1):
      eta[k]=np.sum(d[n]*psi[:,k]*np.exp(-2*np.pi*1j*n*nu)) 

    for k in np.arange(k_max+1):
      u=u+eta[k]*Psi[0,k]/eigval[k]

    if verbose==True:
      print 'expansion coeffcients:'+str(eta)
      print 'eigenvalues:'+str(eigval)
      print 'spectral estimate:'+str(u)

    #- compute adjoint source time function (with forward directed time axis)

    adsrc=np.zeros(self.s.nt,dtype=complex)

    for k in np.arange(k_max+1):
      adsrc[n]=adsrc[n]+Psi[0,k]*psi[:,k]*np.exp(-2*np.pi*1j*n*nu)/eigval[k]
   
    adsrc=adsrc/(self.s.dt*u)

    #- differentiate adjoint source time function for measurements on velocity seismograms
    #- The previous development assumes displacement seismograms

    if self.s.integrate==False:
      adsrc[N_0:N_0+N-1]=-np.diff(adsrc[N_0:N_0+N])/self.s.dt
      adsrc[N_0+N]=0.0

    #- transcribe to adjoint source components

    if component=='x':
      self.adsrc_x=np.real(adsrc)
    elif component=='y':
      self.adsrc_y=np.real(adsrc)
    elif component=='z':
      self.adsrc_z=np.real(adsrc)


  #########################################################################
  # make adjoint source for the measurement of an L2 amplitude difference
  #########################################################################

  def make_adsrc_L2amp(self,component,win_left,win_right,taper_width,verbose=False):
    """ compute adjoint source for the measurement of an L2 amplitude difference

    make_adsrc_mtamp(component,win_left,win_right,taper_width,verbose=False)

    component: seismogram component ('x','y' or 'z') for (theta, phi or r)
    win_left: left window boundary [s]
    win_right: right window boundary [s]
    taper_width: width of the cosine taper applied to the seismogram [s]
    verbose: additional information as screen output if 'True'
    """

    #- transcribe the trace to be considered

    if component=='x':
      d=self.s.trace_x
    elif component=='y':
      d=self.s.trace_y
    elif component=='z':
      d=self.s.trace_z
    else:
      print('enter correct component (x, y or z)')
      d=np.zeros(self.s.nt,dtype=float)

    #- compute the cosine taper and apply to seismogram

    N_0=int(round(win_left/self.s.dt))
    N=int(round(win_right/self.s.dt)-N_0)

    taper=np.zeros(self.s.nt,dtype=float)
     
    n=(self.s.t>win_left) & (self.s.t<win_right)
    taper[n]=1.0

    n=self.s.t<=taper_width
    h=0.5*(1.0-np.cos(np.pi*self.s.t[n]/taper_width))
    taper[N_0:N_0+len(h)]=h
    h=0.5*(1.0+np.cos(np.pi*self.s.t[n]/taper_width))
    taper[N_0+N-len(h):N_0+N]=h

    if verbose==True:
	plt.plot(self.s.t,taper)
	plt.grid(True)
	plt.xlabel('time [s]')
	plt.title('cosine taper')
	plt.show()

    d=d*taper

    #- compute adjoint source

    adsrc=np.zeros(self.s.nt,dtype=float)
    
    d=np.cumsum(d)*self.s.dt
    adsrc=-d/(np.sum(d*d)*self.s.dt)

    #- differentiate adjoint source time function for measurements on velocity seismograms
    #- The previous development assumes displacement seismograms

    if self.s.integrate==False:
      adsrc[0:len(d)-1]=-np.diff(adsrc)/self.s.dt
      adsrc[len(d)-1]=0.0

    #- transcribe to adjoint source components

    if component=='x':
      self.adsrc_x=adsrc
    elif component=='y':
      self.adsrc_y=adsrc
    elif component=='z':
      self.adsrc_z=adsrc

  #########################################################################
  # make adjoint source for the measurement of a cross-correlation time shift
  #########################################################################

  def make_adsrc_cc(self,component,win_left,win_right,taper_width,verbose=False):
    """ compute adjoint source for the measurement of a cross-correlation time shift

    make_adsrc_mtamp(component,win_left,win_right,taper_width,verbose=False)

    component: seismogram component ('x','y' or 'z') for (theta, phi or r)
    win_left: left window boundary [s]
    win_right: right window boundary [s]
    taper_width: width of the cosine taper applied to the seismogram [s]
    verbose: additional information as screen output if 'True'
    """

    #- transcribe the trace to be considered

    if component=='x':
      d=self.s.trace_x
    elif component=='y':
      d=self.s.trace_y
    elif component=='z':
      d=self.s.trace_z
    else:
      print('enter correct component (x, y or z)')
      d=np.zeros(self.s.nt,dtype=float)

    #- compute the cosine taper and apply to seismogram

    N_0=int(round(win_left/self.s.dt))
    N=int(round(win_right/self.s.dt)-N_0)

    taper=np.zeros(self.s.nt,dtype=float)
     
    n=(self.s.t>win_left) & (self.s.t<win_right)
    taper[n]=1.0

    n=self.s.t<=taper_width
    h=0.5*(1.0-np.cos(np.pi*self.s.t[n]/taper_width))
    taper[N_0:N_0+len(h)]=h
    h=0.5*(1.0+np.cos(np.pi*self.s.t[n]/taper_width))
    taper[N_0+N-len(h):N_0+N]=h

    if verbose==True:
	plt.plot(self.s.t,taper)
	plt.grid(True)
	plt.xlabel('time [s]')
	plt.title('cosine taper')
	plt.show()

    d=d*taper

    #- compute adjoint source

    adsrc=np.zeros(self.s.nt,dtype=float)
    
    adsrc=d/(np.sum(d*d)*self.s.dt)

    #- differentiate adjoint source time function for measurements on velocity seismograms
    #- The previous development assumes displacement seismograms

    if self.s.integrate==False:
      adsrc[0:len(d)-1]=-np.diff(adsrc)/self.s.dt
      adsrc[len(d)-1]=0.0

    #- transcribe to adjoint source components

    if component=='x':
      self.adsrc_x=adsrc
    elif component=='y':
      self.adsrc_y=adsrc
    elif component=='z':
      self.adsrc_z=adsrc