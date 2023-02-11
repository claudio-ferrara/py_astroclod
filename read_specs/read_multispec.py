import numpy as np

import re
import datetime
import os
from time import gmtime, strftime

import astropy.io.fits as fits

class spettro:

	"""
	---Class spettro---
	---This class is created to store a spectrum read from a "multispec" fits file format---
	---It stores the header, wavelengths and flux---
	"""
	
	def __init__(self, h = [], w = [], f = []):
		self.h = h		
		self.w = w
		self.f = f


class soluzione:

	"""
	---Class soluzione---
	---This class is created to store a wavelength solution read from a "multispec" fits file format---
	---It is used to calculate wavelength arrays for multispec format spectra---
	---INSERT DESCRIPTION OF THE ATTRIBUTES---
	"""

	def __init__(self, ap = [], beam = [], dtype = [], w1 = [], dw = [], nw = [], z = [], aplow = [], aphigh = []):
		self.ap = ap
		self.beam = beam
		self.dtype = dtype
		self.w1 = w1
		self.dw = dw
		self.nw = nw
		self.z = z
		self.aplow = aplow
		self.aphigh = aphigh
	
	def dim(self):	

		"""
		---Method dim---
		---returns the number of orders of the spectrum---
		"""
		return len(self.ap)

def getcaosspectrum(filename):

	"""
	---Function getcaosspectrum---
	---This function reads a spectrum from a fits file and it---
	---returns the spectrum in a "spettro" object---
	---The spectrum has to be written in the multispec (IRAF) format---
	"""

	data = fits.getdata(filename)
	header = fits.getheader(filename)
	#-----reading the solution from the header-----
	num = 0
	solution = ''
	for card in header.cards:
		if header.cards[num][0][:4] == 'WAT2':
			solution+=header.cards[num][1].ljust(68)
		num+=1
	
	#-----initializing a "soluzione" object-----
	caos_sol = soluzione([],[],[],[],[],[],[],[],[])

	counter = 0
	salva = False
	APERTURA = ''
	
	for i in range(len(solution)):
		if solution[i]=='"' and salva == False:
			salva = True
			counter = counter + 1
		elif solution[i]=='"' and salva == True:
			salva = False
			if counter > 0:
				APERTURA = APERTURA.strip()
				columns = APERTURA[1:].split()
				caos_sol.ap.append(int(columns[0]))
				caos_sol.beam.append(int(columns[1]))
				caos_sol.dtype.append(int(columns[2]))
				caos_sol.w1.append(float(columns[3]))
				caos_sol.dw.append(float(columns[4]))
				caos_sol.nw.append(int(columns[5]))
				caos_sol.z.append(float(columns[6]))
				caos_sol.aplow.append(float(columns[7]))
				caos_sol.aphigh.append(float(columns[8]))
				
				APERTURA = ''
		if salva == True:
			APERTURA = APERTURA + solution[i]



	spectrum = spettro([],[],[])

	spectrum.h = header
	#-----iterating all over the orders to construct wavelength arrays-----
	for i in range(caos_sol.dim()):
		spectrum.w.append(np.linspace(caos_sol.w1[i],caos_sol.w1[i] + caos_sol.dw[i]*caos_sol.nw[i],caos_sol.nw[i]))
		spectrum.f.append(data[i])
	
	
	return spectrum

