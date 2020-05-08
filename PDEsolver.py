"""
PDE solver for the monodomian 
equation

"""

from dolfin import *

class PDEsolver():

	def __init__(self,domain, time, Mi, v_=None,params=None):
		# inptu store
		self._domain = domain 
		self._Mi = Mi 
		self._time = time 

		if params != None
			self.params = params 
		else:
			self.params = self.default_parameters

		# set up function space
		# need to implement change order (k)
		k = 1
		V = FunctionSpace(self._domain,"CG",k) 
		self.V = V

		if v_ != None:
			self.v_ = Function(self.V,name='v_')
		else
			self.v_ = v_
		self.v = Function(V,name="v")

	def step(self,interval):
		"""
		solution at one step
		
		"""
	def default_parameters():
	def solve():