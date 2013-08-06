import os
import re
import csv
import glob
import numpy as np
import fileinput

def new_converter():
	
	#scans the directory
	#files = os.listdir(path)
	
	# files that end in .txt are acted upon
	files = list(glob.glob('*.txt'))
	filename = ''.join(files)
	for fh in files:
		txt_file = fh
		#csv_file_new = fh[:21]+'_new.csv'
		csv_file = fh[5:14]+".csv"
		with open(txt_file, 'rb') as txt_in:	
			in_txt = csv.reader(((line.replace('<','') for line in txt_in)),delimiter = ' ', )
			with open(csv_file, 'wb') as csv_out:
				out_csv = csv.writer(csv_out, dialect = 'excel', delimiter = ',' )
			# simplifies functions that are used
		#txt_in = open(txt_file, 'rb')
		#csv_out = open(csv_file, 'wb')
	
		 
		
		
		
			# writes headers
				header = ["start time","stop time","TS","photon flux","photon flux err", 
									"Photon index","Photon index error", "Energy flux (erg/cm2/s)","Energy flux error"]	
				out_csv.writerow(header)
				
				out_csv.writerows(in_txt)
				#txt_in.close()
				#csv_out.close()