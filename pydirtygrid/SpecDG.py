#!/usr/bin/env python
from astropy.table import Table
import numpy as np
from astropy.io import fits
fits.column.ASCII2NUMPY['E'] = fits.column.ASCII2NUMPY['F'] = 'f8'
import matplotlib.pyplot as plt
import string as st
import subprocess
from scipy import interpolate
from scipy.integrate import simps
from PhotDG import PhotDG
import h5py

class SpecDG:
	def __init__(self):
		"""
		Read in the mapping file for all spectra
		
		Returns
        -------
        (self.)plan:	Table
        	table containing the identification of each spectrum, as well as its parameters
        	
        (self.)seds:	list
        	empty list to save the spectra if you want to
		"""
		mapfile = '/astro/dust_kg3/klaw/cloudy2/dirtygrid_db/django/param_table6_combined.tsvx'
		self.plan = Table.read(mapfile, format='ascii')
		self.seds = []
		#self.waves = []
		
	
	def findFile(self,file_id):
		"""
		Returns the full filename for a given GID from the spectrum mapping.
		
		Parameter
        -------
		file_id:	string
			identification number of the spectrum to be read; 'gid' column of the self.plan Table
		
		Returns
        -------
        filepath+filename:	string
        	the absolute path and name of the .fits file
		"""
		alphab = []
		for i in range(26): alphab.append(st.uppercase[i])
		alphab = np.array(alphab)
		ind_id = []
		for i in range(len(file_id)):
			tmp_i = np.where(file_id[i] == alphab)
			if len(tmp_i[0] == 1):
				ind_id.append(tmp_i[0][0])
			else:
				break
		prefix = ''.join(alphab[ind_id])
		filepath = '/astro/dust_kg3/klaw/cloudy2/nasa_fits/'+prefix+'/'+str(file_id[-4:-2])+'/'
		files = []
		proc = subprocess.Popen(['ls', filepath], stdout=subprocess.PIPE)
		for line in proc.stdout.readlines(): files.append(line.rstrip())
		filename = [s for s in files if file_id in s]
		if len(filename) == 0: 
			return 'Void'
		else:
			filename = filename[0]
			return filepath+filename
	

	def findGidFromParam(self, grain, geom, sf_type, metal, age, sfr, tau):
		"""
		Returns the GID given a set of parameters (If you want to plot or something)
		
		Parameters
        -------
        grain:	string
        	type of grain

    	geom:	string
    		geometry

		sf_type: string
			star formation type 

		metal:	float
			metallicity
		
		age:	float
			age of the stellar population
		
		sfr:	float
			star formation rate
			
		tau: 	float
			optical depth
		
		Returns
        -------
        file_gid:	string
        	identification number of the spectrum with the given set of parameters
		"""
		# ToDo: A better way to do this?
		ind = np.where( (self.plan['grain'] == grain) & (self.plan['geom'] == geom) & (self.plan['sf_type'] == sf_type) & (self.plan['metal'] == metal) & (self.plan['age'] == age) & (self.plan['sfr'] == sfr) & (self.plan['tau'] == tau) )
		# Not all spectra available -- Exit if non existant - Until interpolation?
		if len(ind[0]) == 0: return
		# If valid, then go on
		file_gid = self.plan['gid'][ind[0][0]]
		return file_gid
				
	
	def findParamsFromGid(self, thisgid):
		"""
		Returns the parameter values for a given GID (Used later to update the cube)
		 
		Parameters
        -------
        thisgid:	string
        	identification number of the spectrum
        	
        Returns
        -------
        thisgt, thisgm, thissf, thismt, thissa, thissr, thista:		integers
        	indices of the parameter values for the given spectrum, in the multidimensional photometry cube
		"""
		ind_gid = np.where(self.plan['gid'] == thisgid)
		thisgt = self.plan['grain'][ind_gid[0][0]]		# Grain type
		thisgm = self.plan['geom'][ind_gid[0][0]]		# Geometry
		thissf = self.plan['sf_type'][ind_gid[0][0]]	# Star formation type
		thismt = self.plan['metal'][ind_gid[0][0]]		# Metal
		thissa = self.plan['age'][ind_gid[0][0]]		# Stellar age
		thissr = self.plan['sfr'][ind_gid[0][0]]		# Star formation rate
		thista = self.plan['tau'][ind_gid[0][0]]		# Optical depth
		return thisgt, thisgm, thissf, thismt, thissa, thissr, thista
	
	
	def specGet(self, filename):
		"""
		Save a spectrum from the filename
		 
		Parameters
        -------
        filename:	string
        	absolute path to the spectrum you want to read
    	
    	Returns
        -------
        (self.)seds:	list
        	update the list to add the SED
    	(self.)waves:	array
    		the wavelengths to the SED
		"""
		# Read the fits file
		hdulist = fits.open(filename)
		scidata = hdulist[1].data	
		self.seds.append(scidata['Flux'])
		self.waves = scidata['Wavelength']

			
	def specPlot(self, ind=-1):
		"""
		Plot SEDs, either giving a specific index of which if more that one saved, or all of them
		
		Parameters
        -------
        ind:	integer(s)	(Optional)
        	indices of the SED you wish to plot
		"""
		plt.xscale('log')
		plt.yscale('log')
		plt.ylabel('Luminosity  [$erg$ $s^{-1}$ $\mu m^{-1}$]', fontsize=15)
		plt.xlabel('Wavelength  [$\mu m$]', fontsize=15)
		if ind == -1:
			print 'Plotting all SEDs in object'
			for i in range(len(self.seds)): 
				plt.plot(self.waves, self.seds[i])
		else:
			#print 'Plotting just one'
			for i in ind: plt.plot(self.waves, np.array(self.seds[ind]))
		# ToDo: need a legend
		# ToDo: specify an figure/subplot if desired?
		plt.show()

	
	def spec2Phot(self, trans_curve, trans_waves, wave0, energy=1):
		"""
		 Compute the new photometry
		 
		Parameters
        -------
		trans_curve:	array
			transmission curve of the new filter /!\ Alread read, as an array
		
		trans_waves:	array
			wavelengths of the new filter
			
		wave0:	float
			center/effective wavelength of the new filter

		energy:	integer	(Optional)
			specification for integration with energy(=1, default) or photons
			
		Returns
        -------
        newcube:	array(3, 6, 2, 5, 50, 29, 25) (float)
        	new photometry cube
		"""
		phtemp = PhotDG(silent=1)	# Create an instance of the PhotDG class to call functions
		light_speed = 2.998e14		# microns/s
		n_files = len(self.plan['gid'])
		# Create the array to fill
		newcube = np.empty((3,6,2,5,50,29,25))
		# Read in each spectrum
		for i in range(n_files):
			if (i % 1000) == 0: print(i)
			filename = self.findFile(self.plan['gid'][i])
			if filename == 'Void': continue
			hdu = fits.open(filename)
			scidata = hdu[1].data
			spec = scidata['Flux']						# SED
			if i == 0: dgwave = scidata['Wavelength']	# Wavelength (same for all files)
			# Interpolate the spectral response on DG wavelength grid
			# Do it in log space
			func = interpolate.interp1d(trans_waves, np.log10(trans_curve), kind='slinear')
			tr_ind = np.where( (dgwave > np.min(trans_waves)) & (dgwave < np.max(trans_waves)))
			newtran = func(dgwave[tr_ind])
			newtran = 10**newtran
			# Frequencies
			nu0 = light_speed/wave0#[indf]
			nu = light_speed/dgwave[tr_ind]
			# Total response curve, for integration		(Used simps) 
			if energy == 1:
				tot_resp_curve = simps(newtran*nu0/nu, dgwave[tr_ind])
				tmp_phot = simps(spec[tr_ind]*newtran, dgwave[tr_ind]) 
				newphot = tmp_phot / tot_resp_curve
			else:		# In photons. Do the Blackbody too?
				tot_resp_curve = simps(newtran*nu0/nu**2, dgwave[tr_ind])
				tmp_phot = simps(spec[tr_ind]*newtran/nu, dgwave[tr_ind]) 
				newphot = tmp_phot / tot_resp_curve
			# Fill the array 
			curgt, curgm, cursf, curmt, cursa, cursr, curta = self.findParamsFromGid(self.plan['gid'][i])
			indgt, indgm, indsf, indmt, indsa, indsr, indta = phtemp.findIndexFromParams(curgt, curgm, cursf, curmt, cursa, cursr, curta)
			newcube[indgt, indgm, indsf, indmt, indsa, indsr, indta] = newphot
		return newcube
		#ToDo: maybe save sometimes to avoid all recompute if fails somewhere?

##########
	
		
		
		


