"""
Base class for cell models

"""

from collections import OrderedDict

class CellModel():
	def __init__(self):
		self._params = self.default_parameters()
		self._initial_conditions = self.defaul_initial_conditions()


	def default_parameters():
		return 
	def defaul_initial_conditions():
		pass 
	def set_parameters(self):
		pass 
	def set_initial_conditions(self):
		pass 
	def initial_conditions(self):
		pass 
	def paramters(self):
		pass 
	def F(self,v,s,time=None):
		pass
	def I(self,v,s,time=None):
		pass
	def num_states(self):
		pass