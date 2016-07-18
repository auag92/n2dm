"""Holds functions to make plotting of coordinations easier."""
import pickle
import os
import sys
import numpy as np
import warnings
class MC_Result:

	def __init__(self,av_coordi=None,av_energy=None,tgas=None,pgas=None,tamc=None,
	natoms=None,diameter=None,filename=None):
		self.version = 1.0
		self.av_coordi = av_coordi
		self.av_energy = av_energy
		self.natoms = natoms
		self.tgas = tgas
		self.pgas = pgas
		self.tamc = tamc
		self.diameter = diameter
		self.filename = filename
		
		
	def PickleIt(self):
		"""Saves object to pickle, parameters should have been set before"""
		
		if self.filename ==None:
			op = "./result0_N"+str(self.natoms)+".mcresult"
			i=1
			while os.path.exists(op):
				op = "./result"+str(i)+"_N"+str(self.natoms)+".mcresult"
			
			self.filename = op
			
		f = open(self.filename,"w")
		pickle.dump(self,f)
		
		f.close()
		
		
	def LoadPickle(self,filename):
		
		f = open(filename,"r")
		p = pickle.load(f)
		f.close()
		if p.version != self.version:
			warnings.warn("Version used does not match pickled version")
		self.av_coordi = p.av_coordi
		self.av_energy = p.av_energy
		self.natoms = p.natoms
		self.tgas = p.tgas
		self.pgas = p.pgas
		self.diameter = p.diameter
		self.tamc = p.tamc
		self.filename = p.filename
	

		
		
		
