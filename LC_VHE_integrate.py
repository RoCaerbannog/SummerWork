import asciitable
import numpy
import matplotlib as mpt
from matplotlib import pyplot
from prettytable import *
from VOT import *
import timeit

"""
This version is designed to run Jeremy's extrapolation program and give out an integrated flux, errors, and the detection
"""

class flux:

	def __init__(self, redshift = 0.34, eblModel = 'Dominguez', instrument = 'VERITAS', zenith = 20, index = -1.2, index_err = 0.17, N0 = 1000, start_time = 0, stop_time = 1000):	

		self.redshift = redshift
		self.eblModel = eblModel
		self.instrument = instrument
		self.zenith = zenith
		self.index = index
		self.index_err = index_err
		self.N0 = N0
		self.start_time = start_time
		self.stop_time = stop_time

	# this finds the GRB file and reads it 
		#GRB = 'GRBs/'+ Burst +'.csv'
		#burst = asciitable.read(GRB, delimiter = ',')
	
	# this makes the arrays to which we will add values from the GRB table
		#start_time = []
		#stop_time = []
		#Photon_index = []
		#EF_erg = []
		#EFM_erg = []
		#Flux_erg = []
		#One_GeV = []
		#Four_GeV = []
		#One_TeV = []
	# this takes the values found under start time, stop time, etc. and places them in the arrays we made before
		
			
	# these make new arrays. One is the average time and the other is the time between the start time and the average
		#time,time_err = zip(*[(((y-x)/2)+x,(y-x)/2) for x,y in zip(start_time, stop_time)])
	# this changes ergs to GeV
		#Flux_GeV =  [624.150934 * x for x in Flux_erg]
		#self.Flux_GeV = Flux_GeV
	 	#detTimes =[]
 		#dNdE_new = np.array([])
	 	#int_flux = []
		myVOT = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = self.N0, index = self.index, E0 = 1)
		myVOT_P = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = self.N0, index = self.index + self.index_err, E0 = 1)	
		myVOT_M = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = self.N0, index = self.index - self.index_err, E0 = 1)	    
		
		self.myVOT = myVOT
		self.myVOT_P = myVOT_P
		self.myVOT_M = myVOT_M
		
		mymask = (myVOT.VS.EBins > 200) 
		mymask_p = (myVOT_P.VS.EBins > 200)
		mymask_m =  (myVOT_M.VS.EBins > 200)
		
		
		E_new = [200] + myVOT.VS.EBins[mymask]
		E_new_p = [200] + myVOT_P.VS.EBins[mymask_p]
		E_new_m = [200] + myVOT_M.VS.EBins[mymask_m]
		dNdE_new = [np.interp(200, myVOT.VS.EBins, myVOT.VS.dNdE)] + np.array(myVOT.VS.dNdE)[mymask] 
		dNdE_new_p = [np.interp(200, myVOT_P.VS.EBins, myVOT_P.VS.dNdE)] + np.array(myVOT_P.VS.dNdE)[mymask_p]
		dNdE_new_m = [np.interp(200, myVOT_M.VS.EBins, myVOT_M.VS.dNdE)] + np.array(myVOT_M.VS.dNdE)[mymask_m]
		
		# percent error and integrated flux
		int_flux = np.trapz(dNdE_new, x = E_new)
		int_flux_p = np.trapz(dNdE_new_p, x = E_new_p)
		int_flux_m = np.trapz(dNdE_new_m, x = E_new_m)
		pp_error = (np.absolute(int_flux_p - int_flux)/int_flux)*100
		pm_error = (np.absolute(int_flux_m - int_flux)/int_flux)*100
		self.int_flux = int_flux
		self.int_flux_p = int_flux_p
		self.int_flux_m = int_flux_m
		self.pp_error = pp_error
		self.pm_error = pm_error
		
		
		time_bin = stop_time - start_time
		
		#checks how long detection takes
		crabFlux = 100*myVOT.rate*60./myVOT.VR.crabRate
		self.crabFlux = crabFlux
		crabFlux_p = 100*myVOT_P.rate*60./myVOT_P.VR.crabRate
		crabFlux_m = 100*myVOT_M.rate*60./myVOT_M.VR.crabRate
		detTime = np.interp([crabFlux*0.01], myVOT.VR.SensCurve[:,0], myVOT.VR.SensCurve[:,1])*60
		detTime_p = np.interp([crabFlux_p*0.01], myVOT_P.VR.SensCurve[:,0], myVOT_P.VR.SensCurve[:,1])*60
		detTime_m = np.interp([crabFlux_m*0.01], myVOT_M.VR.SensCurve[:,0], myVOT_M.VR.SensCurve[:,1])*60
		self.detTime = detTime
		self.detTime_p = detTime_p
		self.detTime_m = detTime_m
		# gives an output of observed or not observed depending on whether the detection time is less than or greater than the time bin (exposure length)
		if time_bin >= detTime[0]:
			detection = 'Observed'
		else:
			detection = 'Not Observed'
		if time_bin >= detTime_p[0]:
			detection_p = 'Observed'
		else:
			detection_p = 'Not Observed'
		if time_bin >= detTime_m[0]:
			detection_m = 'Observed'
		else:
			detection_m = 'Not Observed'
		self.detection = detection	
		self.detection_p = detection_p
		self.detection_m = detection_m
		
		
	def make_table(self):
		
		t = PrettyTable(["Photon Index","Integrated Flux","Percent Error", "detection"])
		#t.align["Photon Index"] = "l" #left align photon indicies  
		t.add_row([self.index, self.int_flux, 'N/A', self.detection])
		t.add_row([self.index + self.index_err, self.int_flux_p, str(self.pp_error)+'%', self.detection_p])
		t.add_row([self.index - self.index_err, self.int_flux_m, str(self.pm_error)+'%', self.detection_m])
		print t
	
	def PI_plot(self):
		fig2 = pyplot.figure(figsize=(16,8))
		fig2ax1 = fig2.add_subplot(111)
		fig2ax1.set_ylabel(r'Differential Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]')
		fig2ax1.set_xlabel('Energy [GeV]')
		l1 = fig2ax1.loglog(self.myVOT.VS.EBins, self.myVOT.VS.dNdE, 'b')
		l2 = fig2ax1.loglog(self.myVOT_P.VS.EBins, self.myVOT_P.VS.dNdE, 'g')
		l3 = fig2ax1.loglog(self.myVOT_M.VS.EBins, self.myVOT_M.VS.dNdE, 'r')
		fig2ax1.set_ylim((1e-16,1e-6))
		fig2ax1.set_xlim([0,1000])
		fig2ax2 = fig2ax1.twinx()
		l0 = fig2ax2.loglog((10**self.myVOT.VR.EACurve[0:,0])[::2],(self.myVOT.VR.EACurve[0:,1])[::2],'ko')
		fig2ax2.set_ylabel(r'Effective Area [cm$^2$]')
		fig2ax2.axvline(10**self.myVOT.VR.EASummary['minSafeE'],color='k')
		fig2ax2.axvline(10**self.myVOT.VR.EASummary['maxSafeE'],color='k')
		fig2.legend([l1[0], l2[0], l3[0], l0[0]], ["$\Gamma =$" + str(self.index),"$\Gamma =$" + str(self.index + self.index_err),"$\Gamma =$" + str(self.index - self.index_err), 'EA'], 1)
		pyplot.show()
	
	def N_plot(self):
		fig2 = pyplot.figure(figsize=(16,8))
		fig2ax1 = fig2.add_subplot(111)
		fig2ax1.set_ylabel(r'Differential Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]')
		fig2ax1.set_xlabel('Energy [GeV]')
		l1 = fig2ax1.loglog(self.myVOT.VS.EBins, self.myVOT.VS.dNdE_absorbed, 'b')
		l2 = fig2ax1.loglog(self.myVOT_P.VS.EBins, self.myVOT_P.VS.dNdE_absorbed, 'g')
		l3 = fig2ax1.loglog(self.myVOT_M.VS.EBins, self.myVOT_M.VS.dNdE_absorbed, 'r')
		fig2ax1.set_ylim((1e-16,1e-6))
		fig2ax1.set_xlim([0,1000])
		fig2ax2 = fig2ax1.twinx()
		l0 = fig2ax2.loglog((10**self.myVOT.VR.EACurve[0:,0])[::2],(self.myVOT.VR.EACurve[0:,1])[::2],'ko')
		fig2ax2.set_ylabel(r'Effective Area [cm$^2$]')
		fig2ax2.axvline(10**self.myVOT.VR.EASummary['minSafeE'],color='k')
		fig2ax2.axvline(10**self.myVOT.VR.EASummary['maxSafeE'],color='k')
		fig2.legend([l1[0], l2[0], l3[0], l0[0]], ["$\Gamma =$" + str(self.index),"$\Gamma =$" + str(self.index + self.index_err),"$\Gamma =$" + str(self.index - self.index_err), 'EA'], 1)
		pyplot.show()
		