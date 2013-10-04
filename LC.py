import asciitable
import numpy
import matplotlib as mpt
from matplotlib import pyplot

def take_data(Burst = '130427A'):	
	import asciitable
	GRB = 'GRBs/'+ Burst +'.csv'
	burst = asciitable.read(GRB, delimiter = ',')
	start_time = []
	stop_time = []
	Flux_erg = []
	for x_1 in burst['start time']:
		  start_time.append(x_1)
	for x_2 in burst['stop time']:
		  stop_time.append(x_2)
	for y in burst['Energy flux (erg/cm2/s)']:
		Flux_erg.append(y)
	time,time_err = zip(*[(((y-x)/2)+x,(y-x)/2) for x,y in zip(start_time, stop_time)])
	Flux_GeV =  [624.150934 * x for x in Flux_erg]
	fig2 = pyplot.figure(figsize=(16,8))
	fig2ax1 = fig2.add_subplot(111)
	fig2ax1.set_ylabel(r'Differential Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]')
	fig2ax1.set_xlabel('Time [s]')
	fig2ax1.errorbar(time, Flux_GeV, xerr=time_err, fmt='o')
	fig2ax1.loglog(time, Flux_GeV)
	pyplot.show()
	
def Plot():
	fig2 = pyplot.figure(figsize=(16,8))	
	ax1 = fig1.add_subplot(111)
	ax1.set_xlabel('Time [s]')
	ax1.set_ylabel(r'Differential Flux [cm$^{-2}$ s$^{-1}$ GeV$^{-1}$]')
	ax1.loglog(burst['start time'], burst['Energy flux(GeV/cm2/s)'])