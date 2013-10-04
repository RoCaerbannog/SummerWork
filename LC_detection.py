import asciitable
import numpy
import matplotlib as mpt
from matplotlib import pyplot
from VOT import *
import timeit

class take_data:

	def __init__(self, Burst = '130427A', redshift = 0.34, eblModel = 'Dominguez', instrument = 'VERITAS', zenith = 20):	

		self.redshift = redshift
		self.eblModel = eblModel
		self.instrument = instrument
		self.zenith = zenith
		

	# this finds the GRB file and reads it 
		GRB = 'GRBs/'+ Burst +'.csv'
		burst = asciitable.read(GRB, delimiter = ',')
	
	# this makes the arrays to which we will add values from the GRB table
		start_time = []
		stop_time = []
		Flux_erg = []
		Photon_index = []
		EF_erg = []
		EFM_erg = []
	
	# this takes the values found under start time, stop time, etc. and places them in the arrays we made before
		for x_1 in burst['start time']:
			start_time.append(x_1)
		for x_2 in burst['stop time']:
			stop_time.append(x_2)
		for y in burst['Energy flux (erg/cm2/s)']:
			Flux_erg.append(y)
		for p in burst['Photon index']:
			Photon_index.append(p)
		for a in burst['Energy flux plus']:
			EF_erg.append(a)
		for b in burst['Energy flux minus']:
			EFM_erg.append(b)
			
	# these make new arrays. One is the average time and the other is the time between the start time and the average
		time,time_err = zip(*[(((y-x)/2)+x,(y-x)/2) for x,y in zip(start_time, stop_time)])
	# this changes ergs to GeV
		Flux_GeV =  [624.150934 * x for x in Flux_erg]
		self.Flux_GeV = Flux_GeV
	 
 		dNdE_new = np.array([])
	 	int_flux = []
		for i, v in zip(Flux_GeV, Photon_index,):
			myVOT = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 10000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			
			# makes a boolean array (Trues and Falses) of all the energy bins depending on whether E > 200 GeV
			mymask = myVOT.VS.EBins > 200
			#print np.shape(mymask),type(myVOT.VS.EBins),type(myVOT.VS.dNdE)
			#print mask
			
			# makes a new array of energies E >= 200 GeV
			E_new = [200] + myVOT.VS.EBins[mymask]
			# makes a new differential flux with an interpolated value at 200 GeV and all values above (not interpolated)
			dNdE_new = [np.interp(200, myVOT.VS.EBins, myVOT.VS.dNdE)] + np.array(myVOT.VS.dNdE)[mymask]
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			int_flux.append(np.trapz(dNdE_new, x = E_new))
		
		Flux_GeV_error =  [624.150934 * x for x in EF_erg]
		Flux_error = np.array(Flux_GeV_error) - np.array(Flux_GeV)
		self.Flux_error = Flux_error
	 # this takes care of the flux error
 		dNdE_new_error = np.array([])
	 	int_flux_error = []
		for i, v in zip(Flux_GeV_error, Photon_index,):
			myVOT_error = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 10000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			
			# makes a boolean array (Trues and Falses) of all the energy bins depending on whether E > 200 GeV
			mymask = myVOT_error.VS.EBins > 200
			#print np.shape(mymask),type(myVOT.VS.EBins),type(myVOT.VS.dNdE)
			#print mask
			
			# makes a new array of energies E >= 200 GeV
			E_new_error = [200] + myVOT_error.VS.EBins[mymask]
			# makes a new differential flux with an interpolated value at 200 GeV and all values above (not interpolated)
			dNdE_new_error = [np.interp(200, myVOT.VS.EBins, myVOT_error.VS.dNdE)] + np.array(myVOT_error.VS.dNdE)[mymask]
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			int_flux_error.append(np.trapz(dNdE_new_error, x = E_new_error))
		
			# 
		error = np.array(int_flux_error)-np.array(int_flux) 
		# assigns instances(attributes) to the take_data class so that it can be called in the notebook
		self.time = time
		self.int_flux = int_flux
		self.error = error
		self.time_err = time_err
		self.E_new = E_new
		self.myVOT = myVOT
		
	
		#print int_flux	
		#fig2 = pyplot.figure(figsize=(16,8))
		#fig2ax1 = fig2.add_subplot(111)
		#fig2ax1.set_ylabel(r'Integral Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]')
		#fig2ax1.set_xlabel('Time [s]')
		#fig2ax1.errorbar(time, int_flux, xerr=time_err, fmt='o')
		#fig2ax1.loglog(time, int_flux)
		#legend()
		#pyplot.show()