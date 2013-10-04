import asciitable
import numpy
import matplotlib as mpt
from matplotlib import pyplot
from VOT import *
import timeit
import operator

'''
Ian's NASA Summer 2013 Work

To begin, I just learned python this summer so I may be spelling out some things that should
would be plainly obvious to any programmer. This will help me to make better sense of what it
is that I have helped to create. 
'''
'''
This is the meat and bones of the tool. It uses what is called object-oriented programming
to define an object. What this means is that the object or class take_data has certain properties 
or iterables. These properties have an effect on the processes
'''
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
		PI = []
		PI_err = []
		EF_erg = []
		PF = []
		PF_err = []
		EF_erg_err = []
		Flux_erg = []
		Flux_erg_error = []
		One_GeV = []
		One_GeV_P = []
		One_GeV_M = []
		Four_GeV = []
		Four_GeV_P = []
		Four_GeV_M = []
		One_TeV = []
		One_TeV_P = []
		One_TeV_M = []
	# this takes the values found under start time, stop time, etc. and places them in the arrays we made before
		for x_1 in burst['start time']:
			start_time.append(x_1)
		for x_2 in burst['stop time']:
			stop_time.append(x_2)
		for y in burst['photon flux']:
			PF.append(y)
		for y_1 in burst['photon flux err']:
			PF_err.append(y_1)
		for p in burst['Photon index']:
			PI.append(p)
		for p_1 in burst['Photon index error']:
			PI_err.append(p_1)
		for E_1 in burst['Energy flux (erg/cm2/s)']:
			Flux_erg.append(E_1)
		for E_e in burst['Energy flux error']:
			Flux_erg_error.append(E_e)
		
		self.PI = PI
		self.PI_err = PI_err
	# these make new arrays. One is the average time and the other is the time between the start time and the average
		time,time_err = zip(*[(((y-x)/2)+x,(y-x)/2) for x,y in zip(start_time, stop_time)])
		self.time = time
		self.time_err = time_err
	# this changes ergs to GeV
		Flux_GeV =  [624.150934 * x for x in Flux_erg]
		self.Flux_GeV = Flux_GeV
		Flux_GeV_error =  [624.150934 * x for x in Flux_erg_error] 
		self.Flux_GeV_error = Flux_GeV_error
		#The energy range for the LAT goes from 0.1 GeV to 100 GeV
		E2 =  100 # GeV
		E1 = 0.1 # GeV
		Norm = []
		Norm_P = []
		Norm_M = []
		for F, F_err, I in zip(PF, PF_err, PI):
			Norm.append((F/(E2-E1))*((E2**(1-I)-E1**(1-I))/(1-I)))
			Norm_P.append(((F+F_err)/(E2-E1))*((E2**(1-I)-E1**(1-I))/(1-I)))
			Norm_M.append(((F-F_err)/(E2-E1))*((E2**(1-I)-E1**(1-I))/(1-I)))
		self.Norm = Norm
		self.Norm_P = Norm_P
		self.Norm_M = Norm_M
	 	detTimes =[]
	 	detTimes_p = []
	 	detTimes_m = []
 		dNdE_new = np.array([])
	 	int_flux = []
	 	self.int_flux = int_flux
	 	int_flux_p = []
	 	self.int_flux_p = int_flux_p
	 	int_flux_m = []
	 	self.int_flux_m = int_flux_m
	 	crab_flux = []
	 	crab_flux_p = []
	 	crab_flux_m = []
	 	self.crab_flux = crab_flux
	 	self.crab_flux_p = crab_flux_p
	 	self.crab_flux_m = crab_flux_m
		for i, v in zip(Norm, PI):
			myVOT = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			
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
		
		for i, v in zip(Norm_P, PI):
			myVOT_PFP = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			
			# makes a boolean array (Trues and Falses) of all the energy bins depending on whether 100 GeV > E > 0.1 GeV
			mymask_p = (myVOT_PFP.VS.EBins > 200)
		
			
			# makes a new array of energies above 200 GeV
			E_new_p = [200] + myVOT_PFP.VS.EBins[mymask_p]
			# makes a new differential flux with an interpolated value at 200 GeV and all values above (not interpolated)
			dNdE_new_p = [np.interp(200, myVOT_PFP.VS.EBins, myVOT_PFP.VS.dNdE)] + np.array(myVOT_PFP.VS.dNdE)[mymask_p]
			# interpolates the differential flux at 1 GeV, 400 GeV and 1 TeV 
			One_GeV_P.append(np.interp(1,myVOT_PFP.VS.EBins, myVOT_PFP.VS.dNdE))
			Four_GeV_P.append(np.interp(400,myVOT_PFP.VS.EBins, myVOT_PFP.VS.dNdE))
			One_TeV_P.append(np.interp(1000,myVOT_PFP.VS.EBins, myVOT_PFP.VS.dNdE))
			self.One_GeV_P = One_GeV_P
			self.Four_GeV_P = Four_GeV_P
			self.One_TeV_P = One_TeV_P
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			int_flux_p.append(np.trapz(dNdE_new_p, x = E_new_p))
			
			
			#checks how long detection takes
			
			crabFlux_p = 100*myVOT_PFP.rate*60./myVOT_PFP.VR.crabRate
			crab_flux_p.append(myVOT_PFP.VR.crabRate)
			detTime_P = np.interp([crabFlux_p*0.01], myVOT_PFP.VR.SensCurve[:,0], myVOT_PFP.VR.SensCurve[:,1])
			detTimes_p.append(detTime_P[0])
			self.detTimes_p = detTimes_p
		for i, v in zip(Norm_M, PI):
			myVOT_PFM = VOT("custom", eMin = 0.001, emax = 100000, Nbins= 100000, redshift = self.redshift, eblModel = self.eblModel, instrument = self.instrument,
				    zenith = self.zenith , spectralModel = 'PowerLaw', N0 = i, index = v, E0 = 1)
			
			# makes a boolean array (Trues and Falses) of all the energy bins depending on whether 100 GeV > E > 0.1 GeV
			mymask_m =  (myVOT_PFM.VS.EBins > 200)
			#print np.shape(mymask),type(myVOT.VS.EBins),type(myVOT.VS.dNdE)
			#print mask
			
			# makes a new array of energies above 200 GeV
			E_new_m = [200] + myVOT_PFM.VS.EBins[mymask_m]
			# makes a new differential flux with an interpolated value at 200 GeV and all values above (not interpolated)
			dNdE_new_m = [np.interp(200, myVOT_PFM.VS.EBins, myVOT_PFM.VS.dNdE)] + np.array(myVOT_PFM.VS.dNdE)[mymask_m]		
			# interpolates the differential flux at 1 GeV, 400 GeV and 1 TeV 
			One_GeV_M.append(np.interp(1,myVOT_PFM.VS.EBins, myVOT_PFM.VS.dNdE))
			Four_GeV_M.append(np.interp(400,myVOT_PFM.VS.EBins, myVOT_PFM.VS.dNdE))
			One_TeV_M.append(np.interp(1000,myVOT_PFM.VS.EBins, myVOT_PFM.VS.dNdE))
			self.One_GeV_M = One_GeV_M
			self.Four_GeV_M = Four_GeV_M
			self.One_TeV_M = One_TeV_M
			
			# Integrates differential flux with respect to energy bins using the trapezoid rule
			int_flux_m.append(np.trapz(dNdE_new_m, x = E_new_m))
			
			#checks how long detection takes
			
			crabFlux_m = 100*myVOT_PFM.rate*60./myVOT_PFM.VR.crabRate
			crab_flux_m.append(myVOT_PFM.VR.crabRate)
			detTime_m = np.interp([crabFlux_m*0.01], myVOT_PFM.VR.SensCurve[:,0], myVOT_PFM.VR.SensCurve[:,1])
			detTimes_m.append(detTime_m[0])
			self.detTimes = detTimes_m
		perror = []
		self.perror= perror
		for err, perr in zip(int_flux, int_flux_p):
			perror.append(perr - err)
		merror = []
		self.merror = merror
		for err, merr in zip(int_flux, int_flux_m):
			merror.append(err-merr)
		
		
		One_GeV_M_error,Four_GeV_M_error, One_TeV_M_error = zip(*[((x-s),(y-t),(z-u)) for s,t,u,x,y,z in zip(One_GeV_M, Four_GeV_M, One_TeV_M, One_GeV, Four_GeV, One_TeV)])
		One_GeV_P_error,Four_GeV_P_error, One_TeV_P_error = zip(*[((s-x),(t-y),(u-z)) for s,t,u,x,y,z in zip(One_GeV_M, Four_GeV_M, One_TeV_M, One_GeV, Four_GeV, One_TeV)])
		self.One_GeV_M_error = One_GeV_M_error
		self.Four_GeV_M_error = Four_GeV_M_error
		self.One_TeV_M_error = One_TeV_M_error
		self.One_GeV_P_error = One_GeV_P_error
		self.Four_GeV_P_error = Four_GeV_P_error
		self.One_TeV_P_error = One_TeV_P_error
	def intplot(self):
		fig2 = pyplot.figure(figsize=(16,8))
		fig2ax1 = fig2.add_subplot(111)
		fig2ax2 = fig2.add_subplot(112)
		fig2ax1.set_ylabel(r'Integral Flux [cm$^{-2}$ s$^{-1}$ ]')
		fig2ax1.set_xlabel('Time [s]')
		fig2ax1.set_xscale('log')
		fig2ax1.set_yscale('log')
		msk1 = 2*np.array(self.time_err) > np.array(self.detTimes)
		msk2 = 2*np.array(self.time_err) < np.array(self.detTimes)
		l1 = fig2ax1.errorbar(np.array(self.time)[msk1], np.array(self.int_flux)[msk1], xerr = np.array(self.time_err)[msk1],
			 yerr = [np.array(self.merror)[msk1], np.array(self.perror)[msk1]], fmt='+', color = 'g')
		l11 = fig2ax1.errorbar(np.array(self.time)[msk2], np.array(self.int_flux)[msk2], xerr = np.array(self.time_err)[msk2],
			 yerr = [np.array(self.merror)[msk2], np.array(self.perror)[msk2]], fmt='+', color = 'r')
		#l2 = fig2ax1.scatter(np.array(self.time)[msk1], np.array(self.int_flux)[msk1], color = 'g')
		#l22 = fig2ax1.scatter(np.array(self.time)[msk2], np.array(self.int_flux)[msk2], color = 'r')
		fig2ax2 = fig2ax1.twinx()
		fig2.legend([l1, l11], ['Detectedable', 'Undetectedable'], 1)
		#fig2ax1.set_xlim(1,1e6)
		#fig2ax1.set_ylim(1e-12,1e-2)
		fig2ax1.axvline(100, color = 'k')
		pyplot.show()
	def dnplot(self):
		fig2 = pyplot.figure(figsize=(16,8))
		fig2ax1 = fig2.add_subplot(111)
		fig2ax1.set_ylabel(r'Differential Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]')
		fig2ax1.set_xlabel('Time [s]')
		fig2ax1.set_xscale('log')
		fig2ax1.set_yscale('log')
		timer = np.arange(1, 10000, 1)
		x = 1 # 1 GeV
		y = 400 # 400 GeV
		z = 1000 # 1000 GeV
		crab_flux_1 = 2.83e-11*(x/1000.)**(-2.62) 
		crab_flux_2 = 2.83e-11*(y/1000.)**(-2.62) 
		crab_flux_3 = 2.83e-11*(z/1000.)**(-2.62) 
		cflux1 = []
		cflux2 = []
		cflux3 = []
		for t in timer:
   			cflux1.append(crab_flux_1)
    		cflux2.append(crab_flux_2)
    		cflux3.append(crab_flux_3)
		
		l20 = fig2ax1.errorbar(np.array(self.time), np.array(self.One_GeV), xerr = self.time_err, yerr = [np.array(self.One_GeV_M_error), np.array(self.One_GeV_P_error)], fmt = '+', color = 'b')
		l21 = fig2ax1.errorbar(np.array(self.time), np.array(self.Four_GeV), xerr = self.time_err, yerr = [np.array(self.Four_GeV_M_error), np.array(self.Four_GeV_P_error)], fmt = '+', color = 'g')
		l22 = fig2ax1.errorbar(np.array(self.time), np.array(self.One_TeV), xerr = self.time_err, yerr = [np.array(self.One_TeV_M_error), np.array(self.One_TeV_P_error)], fmt = '+', color = 'r')
		l31 = fig2ax1.plot(timer, cflux1)
		l32 = fig2ax1.plot(timer, cflux2)
		l33 = fig2ax1.plot(timer, cflux3)
		
	
	
		fig2.legend([l20[0],l21[0],l22[0]], ['1 GeV', '400 GeV','1 TeV '], 1)
		
		fig2ax1.set_xlim(1e-2,1e6)
		fig2ax1.set_ylim(1e-15,1e4)
		fig2ax1.axvline(100, color = 'k')
		pyplot.show()
	

	
		