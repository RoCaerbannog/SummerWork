import asciitable
import numpy
import matplotlib as mpt
import ephem
from matplotlib import pyplot
from VOT import *
import timeit


class take_data:

	def __init__(self, Burst = '130427A', redshift = 0.34, eblModel = 'Dominguez', instrument = 'VERITAS', date = '2013/4/25 08:05:12', GC = [173.1370800,27.6990200]):	

		self.redshift = redshift
		self.eblModel = eblModel
		self.instrument = instrument
		self.date = date
		self.GC = GC

	# this finds the GRB file and reads it 
		GRB = 'GRBs/'+ Burst +'.csv'
		burst = asciitable.read(GRB, delimiter = ',')
	
	# this makes the arrays to which we will add values from the GRB table
		start_time = []
		stop_time = []
		Photon_index = []
		EF_erg = []
		EFM_erg = []
		Flux_erg = []
		One_GeV = []
		Four_GeV = []
		One_TeV = []
	# this takes the values found under start time, stop time, etc. and places them in the arrays we made before
		for x_1 in burst['start time']:
			start_time.append(x_1)
		for x_2 in burst['stop time']:
			stop_time.append(x_2)
		for y in burst['Energy flux (erg/cm2/s)']:
			Flux_erg.append(y)
		for y_1 in burst['Energy flux error']:
			EF_erg.append(y_1)
		for p in burst['Photon index']:
			Photon_index.append(p)
		
	# these make new arrays. One is the average time and the other is the time between the start time and the average
		time,time_err = zip(*[(((y-x)/2)+x,(y-x)/2) for x,y in zip(start_time, stop_time)])
		"""
		Here we write in the locations of the various telescopes. This uses some tricks of the dictionary capability of python.
		For now I'm just using the TeV Cat locations. Elevations I got from the CTA paper/HAWC website
		Locations:- 110.952158','31.675058
		Veritas: http://tevcat.uchicago.edu/
		H.E.S.S.: http://tevcat.uchicago.edu/
		HAWC: http://www.ice.phys.psu.edu/~deyoung/Home/HAWC.html
		CTA:    http://dx.doi.org/10.1016/j.bbr.2011.03.031
		http://dx.doi.org/10.1016/j.bbr.2011.03.031
		"""
		
		# sets veritas as a telescope
		veritas = ephem.Observer()
		veritas.lon, veritas.lat = '-110.952158','31.675058'
		veritas.elevation = 1275.
		
		#sets HESS as a telescope
		hess = ephem.Observer()
		hess.lon, hess.lat = '16.5','-23.16'
		hess.elevation = 1800.
		
		# sets CTA as a telescope
		cta = ephem.Observer()
		cta.lon, cta.lat = '-25.', '30.'
		cta.elevation = 2000.
		
		#sets HAWC as a telescope
		hawc = ephem.Observer()
		hawc.lon, hawc.lat = '-97.307', '18.995'
		hawc.elevation = 4100.
		
		# The dictionary 
		instruments = { 'VERITAS': veritas, 'HAWC': hawc, 'HESS': hess, 'CTA': cta}
		# Sets a Fixed Astronomical Body
		grb = ephem.FixedBody(GC[0],GC[1])
		
		# Takes the difference between the time bins in order to create steps 
		tdiff = [time[n]-time[n-1] for n in range(1,len(time))]
		
		myinstrument = instruments[instrument]
		
		dates = []
		zenith = []
		moon = []
		myinstrument.date = date
		
		for t,t_diff in zip(time, tdiff):  
   			grb.compute(veritas)
    		m = ephem.Moon(veritas)
    		#print veritas.date,grb.alt
    		dates = np.append(dates,myinstrument.date)
    		zenith = np.append(zenith, 90 + np.degrees(grb.alt))
    		moon = np.append(moon, np.degrees(m.alt))
    		myinstrument.date += t_diff/(24.*60.)
		
		self.zenith = zenith
		self.dates = dates
	# this changes ergs to GeV
		Flux_GeV =  [624.150934 * x for x in Flux_erg]
		self.Flux_GeV = Flux_GeV
	 	detTimes =[]
 		dNdE_new = np.array([])
	 	int_flux = []
	 	crab_flux =[]
	 	self.crab_flux = crab_flux
		for i, v, t in zip(Flux_GeV, Photon_index, zenith):
			myVOT = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = t , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			
			# makes a boolean array (Trues and Falses) of all the energy bins depending on whether 100 GeV > E > 0.1 GeV
			mymask = (myVOT.VS.EBins > 200) 
			#print np.shape(mymask),type(myVOT.VS.EBins),type(myVOT.VS.dNdE)
			#print mask
			
			# makes a new array of energies above 200 GeV
			E_new = [200] + myVOT.VS.EBins[mymask]
			# makes a new differential flux with an interpolated value at 200 GeV and all values above (not interpolated)
			dNdE_new = [np.interp(200, myVOT.VS.EBins, myVOT.VS.dNdE)] + np.array(myVOT.VS.dNdE)[mymask] 
			# interpolates the differential flux at 1 GeV, 400 GeV and 1 TeV 
			One_GeV.append(np.interp(1,myVOT.VS.EBins, myVOT.VS.dNdE))
			Four_GeV.append(np.interp(400,myVOT.VS.EBins, myVOT.VS.dNdE))
			One_TeV.append(np.interp(1000,myVOT.VS.EBins, myVOT.VS.dNdE))
			self.One_GeV = One_GeV
			self.Four_GeV = Four_GeV
			self.One_TeV = One_TeV
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			int_flux.append(np.trapz(dNdE_new, x = E_new))
			
			#checks how long detection takes
			
			crabFlux = 100*myVOT.rate*60./myVOT.VR.crabRate
			crab_flux.append(myVOT.VR.crabRate)
			detTime = np.interp([crabFlux*0.01], myVOT.VR.SensCurve[:,0], myVOT.VR.SensCurve[:,1])
			detTimes.append(detTime[0])
			self.detTimes = detTimes
		print detTimes
		
		Flux_GeV_error =  [624.150934 * x for x in EF_erg]
		EFP_GeV_error = [] 
		for x,y in zip(Flux_GeV_error,Flux_GeV):
			EFP_GeV_error.append(x+y) 
		Flux_error = np.array(Flux_GeV_error) - np.array(Flux_GeV)
		self.Flux_error = Flux_error
	 # this takes care of the flux error
 		dNdE_new_error = np.array([])
	 	int_flux_error = []
		for i, v, t in zip(EFP_GeV_error, Photon_index, zenith):
			myVOT_error = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = t , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			
			# makes a boolean array (Trues and Falses) of all the energy bins depending on whether E > 200 GeV
			mymask = (myVOT_error.VS.EBins > 200)  
			#print np.shape(mymask),type(myVOT.VS.EBins),type(myVOT.VS.dNdE)
			#print mask
			
			# makes a new array of energies 100 GeV > E > 0.1 GeV
			E_new_error = [200] + myVOT_error.VS.EBins[mymask] 
			# makes a new differential flux with an interpolated value at 200 GeV and all values above (not interpolated)
			dNdE_new_error = [np.interp(200, myVOT.VS.EBins, myVOT_error.VS.dNdE)] + np.array(myVOT_error.VS.dNdE)[mymask] 
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			int_flux_error.append(np.trapz(dNdE_new_error, x = E_new))
		
			# 
		error = np.array(int_flux_error)-np.array(int_flux) 
		# assigns instances(attributes) to the take_data class so that it can be called in the notebook
		self.time = time
		self.int_flux = int_flux
		self.int_flux_error = int_flux_error
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