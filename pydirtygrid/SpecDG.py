#!/usr/bin/env python
import tables as pt
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table

class PhotDG:
	def __init__(self):
		# Prevent a pytable buffer size warning
		pt.parameters.BUFFER_TIMES = 1048576*2 # 2GB
		# Load data
		pt_in_filename = "/astro/dust_kg/klaw/projects/dg_interpol3/interpol6.h5"
		pt_in = pt.open_file(pt_in_filename)
		dimensions = pt_in.root.dirty_cube.attrs.dimensions
		assert dimensions == ['grain', 'geom', 'sf_type', 'metal', 'age', 'sfr', 'tau', 'band']
		param_values = {d: pt_in.root.dirty_cube.attrs[d + 's'] for d in dimensions}
		dirty_cube = pt_in.root.dirty_cube
		# Create the objects
		self.seds = []
		self.grid = dirty_cube
		self.param_values = param_values
		print 'SEDs in erg/s/micrometer'		
				
	def photGet(self, grain, geom, sf_type, metal, age, sfr, tau, bands=-1):
		# grain
		fd_grain = [s for s in self.param_values['grain'] if grain in s]
		w_grain = np.where(fd_grain == self.param_values['grain'])
		ind_grain = int(w_grain[0])
		# geometry
		fd_geom = np.where(geom == self.param_values['geom'])
		ind_geom = int(fd_geom[0])
		# sf type
		sf_type = sf_type.lower()
		if sf_type == 'burst':
			ind_sf_type = 0
		else:
			ind_sf_type = 1
		# metal
		fd_metal = np.where((min(np.abs(self.param_values['metal']-metal)) == np.abs(self.param_values['metal']-metal)))
		ind_metal = int(fd_metal[0])
		# age
		fd_age = np.where((min(np.abs(self.param_values['age']-age)) == np.abs(self.param_values['age']-age)))
		ind_age = int(fd_age[0])
		# sfr
		fd_sfr = np.where((min(np.abs(self.param_values['sfr']-sfr)) == np.abs(self.param_values['sfr']-sfr)))
		ind_sfr = int(fd_sfr[0])
		# tau
		fd_tau = np.where((min(np.abs(self.param_values['tau']-tau)) == np.abs(self.param_values['tau']-tau)))
		ind_tau = int(fd_tau[0])
		# bands
		if bands == -1:
		# Output SED - all bands
			self.seds.append(self.grid[ind_grain, ind_geom, ind_sf_type, ind_metal, ind_age, ind_sfr, ind_tau,:])
		else:
		# Output SED - selective bands
		# ToDo: Test if the requested bands exist
			ind_bands=[]
			for i in range(len(bands)): 
				tmp_bd = np.where(self.param_values['band'] == bands[i])
				if len(tmp_bd[0] != 0):
					ind_bands.append(tmp_bd[0][0])
				else:
					print 'Band ' + bands[i] + ' is missing in the DIRTY Grid - Cancelling that last call'
					return
			self.seds.append(self.grid[ind_grain, ind_geom, ind_sf_type, ind_metal, ind_age, ind_sfr, ind_tau, ind_bands])
		
	def photPlot(self, bands =-1, ind=-1):
		plt.ion()
		# Wavelengths
		# FUV NUV u g r i z U B V R I J H K ...
		waves = [0.155, 0.2275, 0.3543, 0.4770, 0.6231, 0.7625, 0.9134, 0.365, 0.445, 0.551, 0.658, 0.806, 1.2483, 1.6313, 2.2010, 3.6, 4.5, 5.8, 8.0, 24, 70, 160, 250, 350, 500]
		# Plot
		plt.xscale('log')
		plt.yscale('log')
		plt.ylabel('Energy  [$erg$ $s^{-1}$ $\mu m^{-1}$]', fontsize=15)
		plt.xlabel('Wavelength  [$\mu m$]', fontsize=15)
		waves=np.array([0.3543, 0.4770, 0.6231, 0.7625, 3.6, 4.5, 5.8, 8.0, 24, 70])
		if ind == -1:
			print 'Plotting all SEDs in object'
			for i in range(len(self.seds)):
				plt.plot(waves, np.array(self.seds[i]))
		else:
			#print 'Plotting just one'
			plt.plot(waves, np.array(self.seds[ind]))

	def photAddNew(self, wave0, band_name, newphot):
		# Here, to write the code to add a new photometry value to DIRTY Grid


##########

