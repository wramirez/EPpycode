
"""
ODE solver 

"""
from dolfin import *

class ODESolver():
	def __init__(self,domain,time,model,I_s=None,params=None):
		# store input
		self._domain = domain
		self._time = time 
		self._model = model #cell model

		# Info from cell model
		self._F = self._model.F 
		self._Iion = self._model.I 
		self._num_states = self._model.num_states()

		# Stimulus
		self._Is = Stimulus(Is)

		if params != None:
			self.params = params 
		else:
			self.params = self.default_parameters()

		# create function space
		self.VS = VectorFunctionSpace(domain,"CG",1,
								dim=self._num_states+1)

		self.vs_ = Function(self.VS,name="vs_")
		self.vs = Function(self.VS,name="vs")

	def time(self):
		return self._time
		
	def step(self,interval):
		pass
	def solution_fields(self):
		pass
	def default_parameters(self):
		pass
	def solve():
		pass