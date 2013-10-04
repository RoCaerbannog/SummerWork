import glob
import os
import numpy as np
from scipy import integrate
from VOTSpectrum import VOTSpectrum
from VOTResponse import VOTResponse

'''
This file contains calculations for the methods used to calculate the significance of detection for HAWC
in the paper 1108.6034.
'''

'''
 The equation for significance is 
 $S = \sqrt{\frac{\Delta t}{FN_{pmt}R_{pmt}}}\int^{E_{max}_E_{min}{dE \frac{dN}{dE}A^{scaler}_{eff}(\theta}}}$ 

'''
def HAWCdet(source = 'custom', redshift = 2, zenith = 20, eblModel =):	
	self.Rate = self.VR.convolveSpectrum(self.VS.EBins, 
                                             self.VS.dNdE_absorbed, 
                                             self.EACurve_interpolated,
                                             10**self.VR.EASummary['minSafeE'], 
                                             10**self.VR.EASummary['maxSafeE'])
	self.DetTime =  25*R_pmt*N_pmt/(self.rate)^2
