#!/usr/bin/env python
#def functionname( parameters ):
#   "function_docstring"
#   function_suite
#   return [expression]
def	plot_and_print_dirtygrid( grain, geom, sf_type, metal, age, sfr, tau, inp='indices'):
	"Plots the SED for the given parameters or indices of the DIRTYGrid, and prints the parameters"
   
	import tables as pt
	import numpy as np
	from astropy.io import fits
	import matplotlib.pyplot as plt

	# you may need this to prevent a pytable buffer size warning
	pt.parameters.BUFFER_TIMES = 1048576*2 # 2GB

	# load data
	pt_in_filename = "/astro/dust_kg/klaw/projects/dg_interpol3/interpol6.h5"
	pt_in = pt.open_file(pt_in_filename)
	dimensions = pt_in.root.dirty_cube.attrs.dimensions
	assert dimensions == ['grain', 'geom', 'sf_type', 'metal', 'age', 'sfr', 'tau', 'band']
	param_values = {d: pt_in.root.dirty_cube.attrs[d + 's'] for d in dimensions}
	
	dirty_cube = pt_in.root.dirty_cube

	# Choice possible between inputs as indices or values
	# If indices: just pick the DIRTYGrid
	if inp == 'indices':
		ind_grain = grain
		ind_geom = geom
		ind_sf_type = sf_type
		ind_metal = metal
		ind_age = age
		ind_sfr = sfr
		ind_tau = sfr
		# Output SED
		sed = dirty_cube[ind_grain, ind_geom, ind_sf_type, ind_metal, ind_age, ind_sfr, ind_tau,:]
	# If it's values: get the closer one if the DIRTYGrid space
	else:
		# grain
		fd_grain = [s for s in param_values['grain'] if grain in s]
		w_grain = np.where(fd_grain == param_values['grain'])
		ind_grain = int(w_grain[0])
		# geometry
		fd_geom = np.where(geom == param_values['geom'])
		ind_geom = int(fd_geom[0])
		# sf type
		sf_type = sf_type.lower()
		if sf_type == 'burst':
			ind_sf_type = 0
		else:
			ind_sf_type = 1
		# metal
		fd_metal = np.where((min(np.abs(param_values['metal']-metal)) == np.abs(param_values['metal']-metal)))
		ind_metal = int(fd_metal[0])
		# age
		fd_age = np.where((min(np.abs(param_values['age']-age)) == np.abs(param_values['age']-age)))
		ind_age = int(fd_age[0])
		# sfr
		fd_sfr = np.where((min(np.abs(param_values['sfr']-sfr)) == np.abs(param_values['sfr']-sfr)))
		ind_sfr = int(fd_sfr[0])
		# tau
		fd_tau = np.where((min(np.abs(param_values['tau']-tau)) == np.abs(param_values['tau']-tau)))
		ind_tau = int(fd_tau[0])
		# Output SED
		sed = dirty_cube[ind_grain, ind_geom, ind_sf_type, ind_metal, ind_age, ind_sfr, ind_tau,:]

	# Print out
	print 'Grain type: '+str(param_values['grain'][ind_grain])
	print 'Geometry: ' + param_values['geom'][ind_geom]
	print 'Star formation type: '+str(param_values['sf_type'][ind_sf_type])
	print 'Metallicity value in the DIRTYGrid: '+str(param_values['metal'][ind_metal])
	print 'Stellar age value in the DIRTYGrid: {0:.3e}'.format(param_values['age'][ind_age])
	if ind_sf_type == 0:
		sigma_m = param_values['sfr'][ind_sfr] * param_values['age'][ind_age]
		print 'Stellar mass surface density value: {0:.3e}'.format(sigma_m)
	else:
		print 'Star formation rate (surface density) value in the DIRTYGrid: {0:.3e}'.format(param_values['sfr'][ind_sfr])
	print 'Optical depth value in the DIRTYGrid: {0:.4e}'.format(param_values['tau'][ind_tau])

	# wavelgth = np.array([0.364,0.442,0.540,0.658,0.806,0.900,1.020,1.220,1.630,2.190,3.6,4.5,5.8,8.0,24,70,160,250,350,500])
	
	# Plot
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_yscale('log')
	ax.set_xscale('log')
	plt.ylabel('Energy  [$erg$ $s^{-1}$ $\mu m^{-1}$]', fontsize=15)
	plt.xlabel('Wavelength  [$\mu m$]', fontsize=15)
	plt.plot([3.6,4.5,5.8,8.0,24,70,160,250,350,500], sed[15:])
	plt.show()

	return;

# TestCall
plot_and_print_dirtygrid(0,1,1,1,0,1,0)
#plot_and_print_dirtygrid(grain='SMC', geom='shell3', sf_type='burst', metal=0.005, age=11e3, sfr=4.e11, tau=0.5, inp='bob')




