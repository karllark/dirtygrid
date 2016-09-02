#!/usr/bin/env python
from astropy.table import Table
import numpy as np
from astropy.io import fits
fits.column.ASCII2NUMPY['E'] = fits.column.ASCII2NUMPY['F8'] = 'f8'
import matplotlib.pyplot as plt
import string as st
import subprocess

class SpecDG:
	def __init__(self):
		mapfile = '/astro/dust_kg3/klaw/cloudy2/dirtygrid_db/django/param_table6_combined.tsvx'
		self.map = Table.read(mapfile, format='ascii')
		self.seds = []
		self.waves = []
		
	def SpecGet(self, grain, geom, sf_type, metal, age, sfr, tau, bands=-1):
		# Find the GID for the given parameter
		# ToDo: A better way to do this?
		ind = np.where( (self.map['grain'] == grain) & (self.map['geom'] == geom) & (self.map['sf_type'] == sf_type) & (self.map['metal'] == metal) & (self.map['age'] == age) & (self.map['sfr'] == sfr) & (self.map['tau'] == tau) )
		file_id = self.map['gid'][ind[0][0]]
		# Find the fits file
		alphab = []
		for i in range(26): alphab.append(st.uppercase[i])
		alphab=np.array(alphab)
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
		filename = [s for s in files if file_id in s][0]
		#filepath = '/astro/dust_kg3/klaw/cloudy2/nasa_fits/BG/01/'
		#filename='dg_BG0195-s3.16228e+06-t10.0-m0.0001-a18-shell3_global_lum.table.fits'
		# Read the fits file
		hdulist = fits.open(filepath+filename)
		scidata = hdulist[1].data	
		self.seds.append(scidata['Flux'])
		self.waves = scidata['Wavelength']
			
	def SpecPlot(self, ind=-1):
		plt.ion()
		# Wavelengths
		# Plot
		plt.xscale('log')
		plt.yscale('log')
		plt.ylabel('Energy  [$erg$ $s^{-1}$ $\mu m^{-1}$]', fontsize=15)
		plt.xlabel('Wavelength  [$\mu m$]', fontsize=15)
		if ind == -1:
			print 'Plotting all SEDs in object'
			for i in range(len(self.seds)):
				plt.plot(self.waves, self.seds[i])
		#else:
			#print 'Plotting just one'
			#plt.plot(self.waves, np.array(self.seds[ind])

