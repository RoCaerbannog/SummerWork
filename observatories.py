import ephem
import numpy as np

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
class VHEtelescope:

	def __init__(self, instrument = 'VERITAS', time = range(5000), GC = [173.1370800,27.6990200], date = '2013/4/25 08:05:12'):		
		
		self.time = time
		self.GC = GC
		self.instrument = instrument
		self.date = date
		
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
		self.instruments = instruments
		# Sets a Fixed Astronomical Body
		grb = ephem.FixedBody(GC[0],GC[1])
		self.grb = grb
		
		# Takes the difference between the time bins in order to create steps 
		#tdiff = [time[n]-time[n-1] for n in range(1,len(time))]
		
		myinstrument = instruments[self.instrument]
		self.myinstrument = myinstrument
		
		dates = np.array([])
		zenith = np.array([])
		moon = np.array([])
		self.myinstrument.date = date 
		self.date = date
		for t in self.time: 
   			self.grb.compute(self.myinstrument)
    		m = ephem.Moon(self.myinstrument)
    		#print veritas.date,grb.alt
    		dates = np.append(dates, self.date)
    		zenith = np.append(zenith, 90 + np.degrees(self.grb.alt))
    		moon = np.append(moon, 90 + np.degrees(m.alt))
    		myinstrument.date += 1./(24.*60.)
		
		self.zenith = zenith
		self.dates = dates
		self.moon = moon
		