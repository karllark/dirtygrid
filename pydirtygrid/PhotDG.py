#!/usr/bin/env python
#import tables as pt
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
import h5py

class PhotDG:
	def __init__(self, silent=0):
		"""
		Read in the DirtyGrid cube as a HDF5

        Returns
        -------
        (self.)seds:	a list that grows with the SED you choose to extract
        (self.)dgrid:	the HDF5 file: 	with attributes corresponding to the parameters values;
        								with each dataset corresponding to each band.
        """
		f = h5py.File("/astro/dust_kg/jchastenet/dirtycube_29mar17.hdf5", "r+")
		# Create the objects
		self.seds = []
		self.dgrid = f
		print 'SEDs in erg/s/micrometer'
		print
		if silent == 0:
			print('-- Discrete parameters')
			print('    - Grain types: ' + str(self.dgrid.attrs['grain']))
			print('    - Geometries: ' + str(self.dgrid.attrs['geom']))
			print('    - Star formation types: ' + str(self.dgrid.attrs['sf_type']))
			print('-- Continuous parameters')
			print('    - Metallicities: ' + str(self.dgrid.attrs['metal'][0]) + ' to ' + str(self.dgrid.attrs['metal'][-1]) + ' (' + str(len(self.dgrid.attrs['metal'])) + ' values)')
			print('    - Stellar ages: ' + str(self.dgrid.attrs['age'][0]) + ' to ' + str(self.dgrid.attrs['age'][-1]) + ' [Myr] (' + str(len(self.dgrid.attrs['age'])) + ' values)')
			print('    - Star formation rates: ' + str(self.dgrid.attrs['sfr'][0]) + ' to ' + str(self.dgrid.attrs['sfr'][-1]) + ' [Solar mass] (' + str(len(self.dgrid.attrs['sfr'])) + ' values)')  
			print(r'       /!\ 5e-11 factor if Star formation type is Continuous [Solar mass / year]')
			print('    - Optical depth: ' + str(self.dgrid.attrs['tau'][0]) + ' to ' + str(self.dgrid.attrs['tau'][-1]) + ' (' + str(len(self.dgrid.attrs['tau'])) + ' values)')
			print(' -- Bands available')
			print(str(self.dgrid.attrs['band']))
			
			
	def photGet(self, grain, geom, sf_type, metal, age, sfr, tau, bands=-1):
		"""
		Save photometry points given a set of parameters. The function allows for :
		 	- values not in the DGrid parameter sampling, and finds the closest point
		 	- specific bands only
		 	
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
    	(self.)seds:	save the corresponding SED
		"""
		
		tmpsed = []
		# grain
		fd_grain = [s for s in self.dgrid.attrs['grain'] if grain in s] #In case the whole name is not known
		w_grain = np.where(fd_grain == self.dgrid.attrs['grain'])
		ind_grain = w_grain[0][0]
		# geometry
		fd_geom = np.where(geom == self.dgrid.attrs['geom'])
		ind_geom = fd_geom[0][0]
		# sf type
		sf_type = sf_type.lower()
		if sf_type == 'burst':
			ind_sf_type = 0
		else:
			ind_sf_type = 1
		# metal
		fd_metal = np.where((min(np.abs(self.dgrid.attrs['metal']-metal)) == np.abs(self.dgrid.attrs['metal']-metal)))
		ind_metal = fd_metal[0][0]
		# age
		fd_age = np.where((min(np.abs(self.dgrid.attrs['age']-age)) == np.abs(self.dgrid.attrs['age']-age)))
		ind_age = fd_age[0][0]
		# sfr
		fd_sfr = np.where((min(np.abs(self.dgrid.attrs['sfr']-sfr)) == np.abs(self.dgrid.attrs['sfr']-sfr)))
		ind_sfr = fd_sfr[0][0]
		# tau
		fd_tau = np.where((min(np.abs(self.dgrid.attrs['tau']-tau)) == np.abs(self.dgrid.attrs['tau']-tau)))
		ind_tau = fd_tau[0][0]
		# bands
		if bands == -1:	bands = self.dgrid.attrs['band']	# Output SED - all bands
		#ToDo: Test if the requested bands exist
		#ToDo: Add the corresponding wavelengths to the extracted SED
		for b in bands: tmpsed.append(self.dgrid[b][ind_grain, ind_geom, ind_sf_type, ind_metal, ind_age, ind_sfr, ind_tau])
		self.seds.append(tmpsed)

	
	def photPlot(self, ind=-1):
		"""
		Plot photometry points, either giving a specific set of bands, or all of them
		                         either giving a specific index of which if more that one saved, or all of them
		Parameters
        -------
        ind:	integer(s)	(Optional)
        	indices of the SED you wish to plot
		"""
		# ToDo: Add colors
		plt.ion()
		# Plot
		plt.xscale('log')
		plt.yscale('log')
		plt.ylabel('Energy  [$erg$ $s^{-1}$ $\mu m^{-1}$]', fontsize=15)
		plt.xlabel('Wavelength  [$\mu m$]', fontsize=15)
		#effwaves= np.array([0.155, 0.2275, 0.3543, 0.4770, 0.6231, 0.7625, 0.9134, 0.365, 0.445, 0.551, 0.658, 0.806, 1.2483, 1.6313, 2.2010, 3.6, 4.5, 5.8, 8.0, 24., 70., 160., 250., 350., 500.,100])
		if ind == -1:
			print 'Plotting all SEDs in object'
			for i in range(len(self.seds)): plt.scatter(self.dgrid.attrs['effwaves'], np.array(self.seds[i]))
		else:
			#print 'Plotting just one'
			for i in ind: plt.scatter(self.dgrid.attrs['effwaves'], np.array(self.seds[i]))

	
	def findIndexFromParams(self, grain, geom, sf_type, metal, age, sfr, tau):
		"""
		Give a set of parameters and find where they will be stored in the DirtyGrid cube
		
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
    	indgt, indgm, indst, indmt, indsa, indsr, indta:	integers
    		the corresponding indices
        """
        
		indgt = np.where(self.dgrid.attrs['grain'] == grain)
		indgm = np.where(self.dgrid.attrs['geom'] == geom)
		indst = np.where(self.dgrid.attrs['sf_type'] == sf_type)
		indmt = np.where(self.dgrid.attrs['metal'] == metal)
		indsa = np.where(self.dgrid.attrs['age'] == age)
		indsr = np.where(self.dgrid.attrs['sfr'] == sfr)
		indta = np.where(self.dgrid.attrs['tau'] == tau)
		return indgt, indgm, indst, indmt, indsa, indsr, indta
		
	
	def photAddNew(self, wave0, band_name, newcube):
		"""
		Create the new dataset and add it to the file
		
		Parameters
        -------
        wave0:	float
        	central wavelength of the new filter
    	
    	band_name:	string
    		name of the new filter
		
		newcube:	array(3, 6, 2, 5, 50, 29, 25) (float)
			new photometry to be added to the file as a new dataset
		"""
		fnew = h5py.File("/astro/dust_kg/jchastenet/dirtycube_29mar17.hdf5", "r+")
		# Update the nominal wavelengths
		oldwaves = f.attrs['effwaves'].tolist()
		oldwaves.append(wave0)
		fnew.attrs['effwaves'] = oldwaves
		# Update the attribute
		oldbands = fnew.attrs['band'].tolist()
		oldbands.append(band_name)
		fnew.attrs['band'] = oldbands
		# Create new dataset
		fnew[band_name] = newcube
		fnew.flush()
		fnew.close()










##########

