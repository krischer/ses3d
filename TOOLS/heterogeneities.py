import os
import numpy as np
import models as m

#==============================================================================================
#- Class for the generation of test heterogeneities.
#==============================================================================================
class heterogeneities(object):
	"""
	Class for the generation of test heterogeneities.
	"""

	def __init__(self, directory, filename="heterogeneity"):
		"""
		__init__(self, directory, filename="heterogeneity")
		"""

		#- Initialise directory from which to read and filename to write. -------------------------
		self.directory = directory
		self.filename = filename

		#- Read drho from the "directory" and turn it into a zero model. --------------------------
		self.d = m.ses3d_model()
		self.d.read(self.directory,"drho","false")

		for k in range(self.d.nsubvol):
			self.d.m[k].v[:,:,:] = 0.0


	def gaussian(self,theta_0=51.75,phi_0=34.0,r_0=6301.0,amp=0.15,sigma_horizontal=100.0,sigma_vertical=30.0,crop_distance=300.0):
		"""
		Heterogeneity in the form of a Gaussian blob.

		model = gaussian(self,theta_0=51.75,phi_0=34.0,r_0=6301.0,amp=0.15,sigma_horizontal=100.0,sigma_vertical=30.0,crop_distance=300.0):

		Colatitude (deg), longitude (deg) and radius (km) of the Gaussian peak: theta_0, phi_0, r_0.
		Amplitude of the Gaussian: amp.
		Standard deviations in vertical and horizontal directions in km: sigma_horizontal, sigma_horizontal.
		"""

		#- Compute Cartesian x, y, z coordinates of the Gaussian peak on the r_0-sphere. ----------

		theta_0=np.pi*theta_0/180.0
		phi_0=np.pi*phi_0/180.0

		x_0=r_0*np.cos(phi_0)*np.sin(theta_0)
		y_0=r_0*np.sin(phi_0)*np.sin(theta_0)
		z_0=r_0*np.cos(theta_0)

		#- March through the subvolumes. ----------------------------------------------------------

		for n in np.arange(self.d.nsubvol):
			
			for k in range (0, np.size(self.d.m[n].r)-1):
				r=self.d.m[n].r[k]
				for j in range (0, np.size(self.d.m[n].lon)-1):
					phi=np.pi*self.d.m[n].lon[j]/180.0
					for i in range (0, np.size(self.d.m[n].lat)-1):
						theta=np.pi*(90.0-self.d.m[n].lat[i])/180.0

						#- Compute Cartesian x, y, z coordinates of that point on the r_0-sphere. -

						x=r_0*np.cos(phi)*np.sin(theta)
						y=r_0*np.sin(phi)*np.sin(theta)
						z=r_0*np.cos(theta)

						#- Horizontal distance to the Gaussian peak on the sphere through r_0. ----

						dist_horizontal=np.sqrt((x-x_0)**2+(y-y_0)**2+(z-z_0)**2)

						#- Apply Gaussian. ---------------------------------------------------------
				
						if  dist_horizontal> crop_distance : 
							self.d.m[n].v[i,j,k]=0
				
						else: self.d.m[n].v[i,j,k]=amp*np.exp(-dist_horizontal/sigma_horizontal-((r-r_0)**2)/sigma_vertical**2)
				
		#- Return and write. ----------------------------------------------------------------------

		if os.path.exists(os.path.join(self.directory,self.filename)):
			print "File already exists. No output written."
		else:
			self.d.write(self.directory,self.filename)


		return self.d



