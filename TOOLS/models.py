# -*- coding: iso-8859-1 -*-

import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.basemap import Basemap
import colormaps as cm
import rotation as rot

#########################################################################
#- parameters for plotting
#########################################################################

res='i'   # c, l, i, h f

#########################################################################
#- define submodel model class
#########################################################################

class ses3d_submodel(object):
  """ class defining an ses3d submodel
  """

  def __init__(self):

    self.lat=np.zeros(1)
    self.lon=np.zeros(1)
    self.r=np.zeros(1)

    self.v=np.zeros((1, 1, 1))


#########################################################################
#- define ses3d model class
#########################################################################

class ses3d_model(object):
  """ class for reading, writing, plotting and manipulating and ses3d model
  """

  def __init__(self):
    """ initiate the ses3d_model class
    initiate list of submodels and read rotation_parameters.txt
    """

    self.nsubvol=0
    self.lat_min=0.0
    self.lat_max=0.0
    self.lon_min=0.0
    self.lon_max=0.0
    self.lat_centre=0.0
    self.lon_centre=0.0
    self.global_regional="global"

    self.m=[]

    #- read rotation parameters

    fid=open('rotation_parameters.txt','r')
    fid.readline()
    self.phi=float(fid.readline().strip())
    fid.readline()
    line=fid.readline().strip().split(' ')
    self.n=np.array([float(line[0]),float(line[1]),float(line[2])])
    fid.close()

  #########################################################################
  #- multiplication with a scalar
  #########################################################################

  def __rmul__(self,factor):
    """ override left-multiplication of an ses3d model by a scalar factor
    """

    res=ses3d_model()

    res.nsubvol=self.nsubvol
    res.lat_min=self.lat_min
    res.lat_max=self.lat_max
    res.lon_min=self.lon_min
    res.lon_max=self.lon_max
    res.lat_centre=self.lat_centre
    res.lon_centre=self.lon_centre
    res.phi=self.phi
    res.n=self.n
    res.global_regional=self.global_regional

    for k in np.arange(self.nsubvol):

      subvol=ses3d_submodel()

      subvol.lat=self.m[k].lat
      subvol.lon=self.m[k].lon
      subvol.r=self.m[k].r
      subvol.v=factor*self.m[k].v

      res.m.append(subvol)

    return res

  #########################################################################
  #- adding two models
  #########################################################################

  def __add__(self,other_model):
    """ override addition of two ses3d models
    """

    res=ses3d_model()

    res.nsubvol=self.nsubvol
    res.lat_min=self.lat_min
    res.lat_max=self.lat_max
    res.lon_min=self.lon_min
    res.lon_max=self.lon_max
    res.lat_centre=self.lat_centre
    res.lon_centre=self.lon_centre
    res.phi=self.phi
    res.n=self.n
    res.global_regional=self.global_regional

    for k in np.arange(self.nsubvol):

      subvol=ses3d_submodel()

      subvol.lat=self.m[k].lat
      subvol.lon=self.m[k].lon
      subvol.r=self.m[k].r
      subvol.v=other_model.m[k].v+self.m[k].v

      res.m.append(subvol)

    return res

  #########################################################################
  #- read a 3D model
  #########################################################################

  def read(self,directory,filename,verbose=False):
    """ read an ses3d model from a file

	read(self,directory,filename,verbose=False):
    """

    #- read block files ====================================================

    fid_x=open(directory+'block_x','r')
    fid_y=open(directory+'block_y','r')
    fid_z=open(directory+'block_z','r')

    if verbose==True:
      print 'read block files:'
      print '\t '+directory+'block_x'
      print '\t '+directory+'block_y'
      print '\t '+directory+'block_z'

    dx=np.array(fid_x.read().strip().split('\n'),dtype=float)
    dy=np.array(fid_y.read().strip().split('\n'),dtype=float)
    dz=np.array(fid_z.read().strip().split('\n'),dtype=float)

    fid_x.close()
    fid_y.close()
    fid_z.close()

    #- read coordinate lines ===============================================

    self.nsubvol=int(dx[0])

    if verbose==True:
      print 'number of subvolumes: '+str(self.nsubvol)

    idx=np.zeros(self.nsubvol,dtype=int)+1
    idy=np.zeros(self.nsubvol,dtype=int)+1
    idz=np.zeros(self.nsubvol,dtype=int)+1

    for k in np.arange(1,self.nsubvol,dtype=int):
      idx[k]=int(dx[idx[k-1]])+idx[k-1]+1
      idy[k]=int(dy[idy[k-1]])+idy[k-1]+1
      idz[k]=int(dz[idz[k-1]])+idz[k-1]+1

    for k in np.arange(self.nsubvol,dtype=int):
      subvol=ses3d_submodel()
      subvol.lat=90.0-dx[(idx[k]+1):(idx[k]+1+dx[idx[k]])]
      subvol.lon=dy[(idy[k]+1):(idy[k]+1+dy[idy[k]])]
      subvol.r  =dz[(idz[k]+1):(idz[k]+1+dz[idz[k]])]
      self.m.append(subvol)

    #- read model volume ==================================================

    fid_m=open(directory+filename,'r')

    if verbose==True:
      print 'read model file: '+directory+filename

    v=np.array(fid_m.read().strip().split('\n'),dtype=float)

    fid_m.close()

    #- assign values ======================================================

    idx=1
    for k in np.arange(self.nsubvol):

      n=int(v[idx])
      nx=len(self.m[k].lat)-1
      ny=len(self.m[k].lon)-1
      nz=len(self.m[k].r)-1

      self.m[k].v=v[(idx+1):(idx+1+n)].reshape(nx,ny,nz)

      idx=idx+n+1

    #- decide on global or regional model==================================

    self.lat_min=90.0
    self.lat_max=-90.0
    self.lon_min=180.0
    self.lon_max=-180.0

    for k in np.arange(self.nsubvol):
      if np.min(self.m[k].lat) < self.lat_min: self.lat_min = np.min(self.m[k].lat)
      if np.max(self.m[k].lat) > self.lat_max: self.lat_max = np.max(self.m[k].lat)
      if np.min(self.m[k].lon) < self.lon_min: self.lon_min = np.min(self.m[k].lon)
      if np.max(self.m[k].lon) > self.lon_max: self.lon_max = np.max(self.m[k].lon)

    if ((self.lat_max-self.lat_min) > 30.0 or (self.lon_max-self.lon_min) > 30.0):
      self.global_regional = "global"

      self.lat_centre = (self.lat_max+self.lat_min)/2.0
      self.lon_centre = (self.lon_max+self.lon_min)/2.0

      self.lat_centre,self.lon_centre = rot.rotate_coordinates(self.n,-self.phi,90.0-self.lat_centre,self.lon_centre)
      self.lat_centre = 90.0-self.lat_centre

    else:
      self.global_regional = "regional"

      self.d_lat=5.0
      self.d_lon=5.0


  #########################################################################
  #- write a 3D model to a file
  #########################################################################

  def write(self,directory,filename,verbose=False):
    """ write ses3d model to a file

    write(self,directory,filename,verbose=False):
    """

    fid_m=open(directory+filename,'w')

    if verbose==True:
      print 'write to file '+directory+filename

    fid_m.write(str(self.nsubvol)+'\n')

    for k in np.arange(self.nsubvol):

      nx=len(self.m[k].lat)-1
      ny=len(self.m[k].lon)-1
      nz=len(self.m[k].r)-1

      fid_m.write(str(nx*ny*nz)+'\n')

      for idx in np.arange(nx):
	for idy in np.arange(ny):
	  for idz in np.arange(nz):

	    fid_m.write(str(self.m[k].v[idx,idy,idz])+'\n')

    fid_m.close()

  #########################################################################
  #- CUt Depth LEvel
  #########################################################################

  def cudle(self,r_min,r_max,verbose=False):
    """ cut out a certain radius range, mostly in order to produce smaller models for vtk

    new_model=cudle(self,depth_min,depth_max)
    """

    m_new=ses3d_model()

    #- march through subvolumes

    for n in np.arange(self.nsubvol):

      if (np.min(self.m[n].r)<=r_min) & (np.max(self.m[n].r)>=r_max):

	if verbose==True:
	  print 'subvolume '+str(n)+': r_min='+str(np.min(self.m[n].r))+' km, r_max='+str(np.max(self.m[n].r))+' km\n'

	subvol=ses3d_submodel()
	subvol.lat=self.m[n].lat
	subvol.lon=self.m[n].lon

	idr=(self.m[n].r>=r_min) & (self.m[n].r<=r_max)
	subvol.r=self.m[n].r[idr]

	idr=idr[1:(len(idr))]
	subvol.v=self.m[n].v[:,:,idr]

	m_new.m.append(subvol)
	m_new.nsubvol=m_new.nsubvol+1

    return m_new

  #########################################################################
  #- convert to vtk format
  #########################################################################

  def convert_to_vtk(self,directory,filename,verbose=False):
    """ convert ses3d model to vtk format for plotting with Paraview

    convert_to_vtk(self,directory,filename,verbose=False):
    """

    #- preparatory steps

    nx=np.zeros(self.nsubvol,dtype=int)
    ny=np.zeros(self.nsubvol,dtype=int)
    nz=np.zeros(self.nsubvol,dtype=int)
    N=0

    for n in np.arange(self.nsubvol):
      nx[n]=len(self.m[n].lat)
      ny[n]=len(self.m[n].lon)
      nz[n]=len(self.m[n].r)
      N=N+nx[n]*ny[n]*nz[n]

    #- open file and write header

    fid=open(directory+filename,'w')

    if verbose==True:
      print 'write to file '+directory+filename

    fid.write('# vtk DataFile Version 3.0\n')
    fid.write('vtk output\n')
    fid.write('ASCII\n')
    fid.write('DATASET UNSTRUCTURED_GRID\n')

    #- write grid points

    fid.write('POINTS '+str(N)+' float\n')

    for n in np.arange(self.nsubvol):

      if verbose==True:
	print 'writing grid points for subvolume '+str(n)

      for i in np.arange(nx[n]):
	for j in np.arange(ny[n]):
	  for k in np.arange(nz[n]):

	    theta=90.0-self.m[n].lat[i]
	    phi=self.m[n].lon[j]

	    #- rotate coordinate system

	    if self.phi!=0.0:
	      theta,phi=rot.rotate_coordinates(self.n,-self.phi,theta,phi)

	    #- transform to cartesian coordinates and write to file

	    theta=theta*np.pi/180.0
	    phi=phi*np.pi/180.0

	    r=self.m[n].r[k]
	    x=r*np.sin(theta)*np.cos(phi);
            y=r*np.sin(theta)*np.sin(phi);
            z=r*np.cos(theta);

	    fid.write(str(x)+' '+str(y)+' '+str(z)+'\n')

    #- write connectivity

    n_cells=0

    for n in np.arange(self.nsubvol):
      n_cells=n_cells+(nx[n]-1)*(ny[n]-1)*(nz[n]-1)

    fid.write('\n')
    fid.write('CELLS '+str(n_cells)+' '+str(9*n_cells)+'\n')

    count=0

    for n in np.arange(self.nsubvol):

      if verbose==True:
	print 'writing conectivity for subvolume '+str(n)

      for i in np.arange(1,nx[n]):
	for j in np.arange(1,ny[n]):
	  for k in np.arange(1,nz[n]):
								# i j k
	    a=count+k+(j-1)*nz[n]+(i-1)*ny[n]*nz[n]-1     	# 0 0 0
	    b=count+k+(j-1)*nz[n]+(i-1)*ny[n]*nz[n]       	# 0 0 1
	    c=count+k+(j)*nz[n]+(i-1)*ny[n]*nz[n]-1       	# 0 1 0
	    d=count+k+(j)*nz[n]+(i-1)*ny[n]*nz[n]         	# 0 1 1
	    e=count+k+(j-1)*nz[n]+(i)*ny[n]*nz[n]-1       	# 1 0 0
	    f=count+k+(j-1)*nz[n]+(i)*ny[n]*nz[n]         	# 1 0 1
	    g=count+k+(j)*nz[n]+(i)*ny[n]*nz[n]-1         	# 1 1 0
	    h=count+k+(j)*nz[n]+(i)*ny[n]*nz[n]           	# 1 1 1

	    fid.write('8 '+str(a)+' '+str(b)+' '+str(c)+' '+str(d)+' '+str(e)+' '+str(f)+' '+str(g)+' '+str(h)+'\n')

      count=count+nx[n]*ny[n]*nz[n]

    #- write cell types

    fid.write('\n')
    fid.write('CELL_TYPES '+str(n_cells)+'\n')

    for n in np.arange(self.nsubvol):

      if verbose==True:
	print 'writing cell types for subvolume '+str(n)

      for i in np.arange(nx[n]-1):
	for j in np.arange(ny[n]-1):
	  for k in np.arange(nz[n]-1):

	    fid.write('11\n')

    #- write data

    fid.write('\n')
    fid.write('POINT_DATA '+str(N)+'\n')
    fid.write('SCALARS scalars float\n')
    fid.write('LOOKUP_TABLE mytable\n')

    for n in np.arange(self.nsubvol):

      if verbose==True:
	print 'writing data for subvolume '+str(n)

      idx=np.arange(nx[n])
      idx[nx[n]-1]=nx[n]-2

      idy=np.arange(ny[n])
      idy[ny[n]-1]=ny[n]-2

      idz=np.arange(nz[n])
      idz[nz[n]-1]=nz[n]-2

      for i in idx:
	for j in idy:
	  for k in idz:

	    fid.write(str(self.m[n].v[i,j,k])+'\n')

    #- clean up

    fid.close()

  #########################################################################
  #- plot horizontal slices
  #########################################################################

  def plot_slice(self,depth,min_val_plot,max_val_plot,colormap='tomo',verbose=False):
    """ plot horizontal slices through an ses3d model

    plot_slice(depth,min_val_plot,max_val_plot,verbose=False):

    depth=depth in km of the slice
    min_val_plot, max_val_plot=minimum and maximum values of the colour scale
    colormap='tomo','mono'
    """

    radius=6371.0-depth

    #- set up a map and colourmap

    if self.global_regional=='regional':
      m=Basemap(projection='merc',llcrnrlat=self.lat_min,urcrnrlat=self.lat_max,llcrnrlon=self.lon_min,urcrnrlon=self.lon_max,lat_ts=20,resolution=res)
      m.drawparallels(np.arange(self.lat_min,self.lat_max,self.d_lon),labels=[1,0,0,1])
      m.drawmeridians(np.arange(self.lon_min,self.lon_max,self.d_lat),labels=[1,0,0,1])
    elif self.global_regional=='global':
      m=Basemap(projection='ortho',lon_0=self.lon_centre,lat_0=self.lat_centre,resolution=res)
      m.drawparallels(np.arange(-80.0,80.0,10.0),labels=[1,0,0,1])
      m.drawmeridians(np.arange(-170.0,170.0,10.0),labels=[1,0,0,1])

    m.drawcoastlines()
    m.drawcountries()

    m.drawmapboundary(fill_color=[1.0,1.0,1.0])

    if colormap=='tomo':
      my_colormap=cm.make_colormap({0.0:[0.1,0.0,0.0], 0.2:[0.8,0.0,0.0], 0.3:[1.0,0.7,0.0],0.48:[0.92,0.92,0.92], 0.5:[0.92,0.92,0.92], 0.52:[0.92,0.92,0.92], 0.7:[0.0,0.6,0.7], 0.8:[0.0,0.0,0.8], 1.0:[0.0,0.0,0.1]})
    elif colormap=='mono':
      my_colormap=cm.make_colormap({0.0:[1.0,1.0,1.0], 0.15:[1.0,1.0,1.0], 0.85:[0.0,0.0,0.0], 1.0:[0.0,0.0,0.0]})

    #- loop over subvolumes

    for k in np.arange(self.nsubvol):

      r=self.m[k].r

      #- check if subvolume has values at target depth

      if (max(r)>=radius) & (min(r)<radius):
        
        r=r[0:len(r)-1]
        idz=min(np.where(min(np.abs(r-radius))==np.abs(r-radius))[0])
        if idz==len(r): idz-=idz

        if verbose==True:
          print 'true plotting depth: '+str(6371.0-r[idz])+' km'

        nx=len(self.m[k].lat)
        ny=len(self.m[k].lon)
        nz=len(self.m[k].r)

        lon,lat=np.meshgrid(self.m[k].lon[0:ny],self.m[k].lat[0:nx])

        #- rotate coordinate system if necessary

        if self.phi!=0.0:

          lat_rot=np.zeros(np.shape(lon),dtype=float)
          lon_rot=np.zeros(np.shape(lat),dtype=float)

          for idx in np.arange(nx):
            for idy in np.arange(ny):

              colat=90.0-lat[idx,idy]
              lat_rot[idx,idy],lon_rot[idx,idy]=rot.rotate_coordinates(self.n,-self.phi,colat,lon[idx,idy])
              lat_rot[idx,idy]=90.0-lat_rot[idx,idy]

          lon=lon_rot
          lat=lat_rot

        #- convert to map coordinates and plot

        x,y=m(lon,lat)
        im=m.pcolor(x,y,self.m[k].v[:,:,idz],cmap=my_colormap,vmin=min_val_plot,vmax=max_val_plot)

        m.colorbar(im,"right", size="3%", pad='2%')
        plt.title(str(depth)+' km')
        plt.show()

  #########################################################################
  #- plot depth to a certain threshold value
  #########################################################################

  def plot_threshold(self,val,min_val_plot,max_val_plot,colormap='tomo',verbose=False):
    """ plot depth to a certain threshold value 'val' in an ses3d model

    plot_threshold(val,min_val_plot,max_val_plot,colormap='tomo',verbose=False):
    val=threshold value
    min_val_plot, max_val_plot=minimum and maximum values of the colour scale
    colormap='tomo','mono'
    """

    #- set up a map and colourmap

    if self.global_regional=='regional':
      m=Basemap(projection='merc',llcrnrlat=self.lat_min,urcrnrlat=self.lat_max,llcrnrlon=self.lon_min,urcrnrlon=self.lon_max,lat_ts=20,resolution=res)
      m.drawparallels(np.arange(self.lat_min,self.lat_max,self.d_lon),labels=[1,0,0,1])
      m.drawmeridians(np.arange(self.lon_min,self.lon_max,self.d_lat),labels=[1,0,0,1])
    elif self.global_regional=='global':
      m=Basemap(projection='ortho',lon_0=self.lon_centre,lat_0=self.lat_centre,resolution=res)
      m.drawparallels(np.arange(-80.0,80.0,10.0),labels=[1,0,0,1])
      m.drawmeridians(np.arange(-170.0,170.0,10.0),labels=[1,0,0,1])

    m.drawcoastlines()
    m.drawcountries()

    m.drawmapboundary(fill_color=[1.0,1.0,1.0])

    if colormap=='tomo':
      my_colormap=cm.make_colormap({0.0:[0.1,0.0,0.0], 0.2:[0.8,0.0,0.0], 0.3:[1.0,0.7,0.0],0.48:[0.92,0.92,0.92], 0.5:[0.92,0.92,0.92], 0.52:[0.92,0.92,0.92], 0.7:[0.0,0.6,0.7], 0.8:[0.0,0.0,0.8], 1.0:[0.0,0.0,0.1]})
    elif colormap=='mono':
      my_colormap=cm.make_colormap({0.0:[1.0,1.0,1.0], 0.15:[1.0,1.0,1.0], 0.85:[0.0,0.0,0.0], 1.0:[0.0,0.0,0.0]})

    #- loop over subvolumes

    for k in np.arange(self.nsubvol):

      depth=np.zeros(np.shape(self.m[k].v[:,:,0]))

      nx=len(self.m[k].lat)
      ny=len(self.m[k].lon)

      #- find depth

      r=self.m[k].r
      r=0.5*(r[0:len(r)-1]+r[1:len(r)])

      for idx in np.arange(nx-1):
	for idy in np.arange(ny-1):

	  n=self.m[k].v[idx,idy,:]>=val
	  depth[idx,idy]=6371.0-np.max(r[n])

      #- rotate coordinate system if necessary

      lon,lat=np.meshgrid(self.m[k].lon[0:ny],self.m[k].lat[0:nx])

      if self.phi!=0.0:

	lat_rot=np.zeros(np.shape(lon),dtype=float)
	lon_rot=np.zeros(np.shape(lat),dtype=float)

	for idx in np.arange(nx):
	  for idy in np.arange(ny):

	    colat=90.0-lat[idx,idy]

	    lat_rot[idx,idy],lon_rot[idx,idy]=rot.rotate_coordinates(self.n,-self.phi,colat,lon[idx,idy])
	    lat_rot[idx,idy]=90.0-lat_rot[idx,idy]

	lon=lon_rot
	lat=lat_rot

	#- convert to map coordinates and plot

      x,y=m(lon,lat)
      im=m.pcolor(x,y,depth,cmap=my_colormap,vmin=min_val_plot,vmax=max_val_plot)

    m.colorbar(im,"right", size="3%", pad='2%')
    plt.title('depth to '+str(val)+' km/s [km]')
    plt.show()
