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

class SpecDG:
	def __init__(self):
		mapfile = '/astro/dust_kg3/klaw/cloudy2/dirtygrid_db/django/param_table6_combined.tsvx'
		self.map = Table.read(mapfile, format='ascii')
		self.seds = []
		self.waves = []
		
	def findFile(file_id):
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
		filename = [s for s in files if file_id in s][0]
		return filepath+filename
	
	def specGet(self, grain, geom, sf_type, metal, age, sfr, tau, bands=-1):
		# Find the GID for the given parameter
		# ToDo: A better way to do this?
		ind = np.where( (self.map['grain'] == grain) & (self.map['geom'] == geom) & (self.map['sf_type'] == sf_type) & (self.map['metal'] == metal) & (self.map['age'] == age) & (self.map['sfr'] == sfr) & (self.map['tau'] == tau) )
		# Not all spectra available -- Exit if non existant - Until interpolation?
		if len(ind[0]) == 0: return
		# If valid, then go on
		file_id = self.map['gid'][ind[0][0]]
		filename = findFile(file_id)
		# Read the fits file
		hdulist = fits.open(filename)
		scidata = hdulist[1].data	
		self.seds.append(scidata['Flux'])
		self.waves = scidata['Wavelength']
			
	def specPlot(self, ind=-1):
		plt.ion()
		# Wavelengths
		# Plot
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
			plt.plot(self.waves, np.array(self.seds[ind]))

	
	def spec2Phot(self, band_names, resp_files, cen_wave, energy=1):
		light_speed = 2.998e14
		n_files = len(self.map['gid'])
		# Read the first fits to get Wavelength grid
#		filename = findFile(self.map['gid'][0])
#		hdulist = fits.open(filename)
#		scidata = hdulist[1].data
#		dgspec = scidata['Flux']
		# Read in each spectrum
		for i in range(n_files):
			filename = findFile(self.map['gid'][i])
			hdulist = fits.open(filename)
			scidata = hdulist[1].data
			spec = scidata['Flux']
			if i == 0: dgwave = scidata['Wavelength']
			for fresp in resp_files:
				indf = 0
				# Interpolate the spectral response on DG wavelength grid
				# Q.: How to know to right way to read the file?
				spec_resp = Table.read(fresp, format='ascii')
				# Do it in log space
				func = interpolate.interp1d(spec_resp['col1'], np.log10(spec_resp['col2']), kind='slinear')							## Prob with column names
				tr_ind = np.where( (dgwave > np.min(spec_resp['col1'])) & (dgwave < np.max(spec_resp['col1']) )							## Prob with column names
				newtran = func(dgwave[tr_ind])
				newtran = 10**newtran
				# Frequencies
				nu0 = light_speed/cen_wave[indf]
				nu = light_speed/dgwave[tr_ind]
				# Total response curve, for integration		(Used simps) 
				if energy == 1:
					tot_resp_curve = simps(newtran*nu0/nu, dgwave[tr_ind])
					tmp_phot = simps(spec[tr_ind]*newtran, dgwave[tr_ind]) 
					newphot = tmp_phot / tot_resp_curves
				## /!\ What about integration in photon?!
				#else:
					#tmp_phot = np.sum(0.5*(spec[0:n_waves-2] + spec[1:n_waves-1]) * 0.5*(new_transm[0:n_waves-2]+new_transm[1:n_waves-1]) *dfreq/afreq) /tot_resp_curves
				
				indf++
##########
## ToDo: Add the new photometry to the DGrid.
	
		
		
		


