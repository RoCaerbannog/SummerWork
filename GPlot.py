import matplotlib as mpt
import numpy as np
from matplotlib import pyplot
from matplotlib import animation
def plotG(n):

	fig2 = pyplot.figure(figsize=(16,8))
	fig2ax1 = fig2.add_subplot(111)
	fig2ax1.set_ylabel(r'Differential Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]')
	fig2ax1.set_xlabel('Energy [GeV]')
	fig2ax1.loglog(n.VS.EBins, n.VS.dNdE)
	l1 = fig2ax1.loglog(n.VS.EBins, n.VS.dNdE_absorbed, color = 'g')
	fig2ax1.set_ylim((1e-16,1e-6))
	fig2ax1.set_xlim([1,1000])
	fig2ax2 = fig2ax1.twinx()
	l2 = fig2ax2.loglog((10**n.VR.EACurve[0:,0])[::2],(n.VR.EACurve[0:,1])[::2],'ko')
	fig2ax2.set_ylabel(r'Effective Area [cm$^2$]')
	fig2ax2.axvline(10**n.VR.EASummary['minSafeE'],color='k')
	fig2ax2.axvline(10**n.VR.EASummary['maxSafeE'],color='k')
	fig2.legend([l1[0], l2[0]], ['z = 0.34', 'EA'], 1)
	pyplot.show()

def plotG_all(n):

	fig2 = pyplot.figure(figsize=(16,8))
	fig2ax1 = fig2.add_subplot(111)
	fig2ax1.set_ylabel(r'Differential Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]')
	fig2ax1.set_xlabel('Energy [GeV]')
	fig2ax1.loglog(n.VS.EBins, n.VS.dNdE)
	l1 = fig2ax1.loglog(n.VS.EBins, n.VS.dNdE_absorbed, color = 'g')
	fig2ax1.set_ylim((1e-16,1e-6))
	fig2ax1.set_xlim([1,1000])
	fig2ax2 = fig2ax1.twinx()
	l2 = fig2ax2.loglog((10**n.VR.EACurve[0:,0])[::2],(n.VR.EACurve[0:,1])[::2],'ko')
	fig2ax2.set_ylabel(r'Effective Area [cm$^2$]')
	fig2ax2.axvline(10**n.VR.EASummary['minSafeE'],color='k')
	fig2ax2.axvline(10**n.VR.EASummary['maxSafeE'],color='k')
	fig2.legend([l1[0], l2[0]], ['z = 0.34', 'EA'], 1)


	
	