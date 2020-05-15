"""
Base class for cell models

"""

from collections import OrderedDict
from dolfin import Expression

class CellModel():
	def __init__(self,params=None,initial_conditions=None):
		self._params = self.default_parameters()
		self._initial_conditions = self.defaul_initial_conditions()


	def default_parameters():
		return OrderedDict()
	def defaul_initial_conditions():
		return OrderedDict()
	def set_parameters(self,**params):
		for param_name,param_value in params.items():
			self._params[param_name] = param_value

	def set_initial_conditions(self,**init):
		for init_name,init_value in init.items():
			self._initial_conditions[init_name] = init_value

	def initial_conditions(self):
		return Expression(list(self._initial_conditions.keys())
					,degree=1,**self._initial_conditions)
	def parameters(self):
		return self._params

	def F(self,v,s,time=None):
		"""
			Return right hand side for state
			variable evolution
		"""
		error("Must define F=F(v,s)")

	def I(self,v,s,time=None):
		"""
			Return the ionic current
		"""
		error("Must define I = I(v,s)")

	def num_states(self):
		"""
		return number of state variables
		"""
		error("Must overload num_states")
	def __str__(self):
		
		return "Some cardiac cell model"