"""
Cardiac Model class
to handle heart EP inputs

"""


class CardiacModel():

	def  __init__(self,domain,time,Mi,
					cell_model,stimulus):

		self._domain = domain 
		self._time = time 
		self.intracellular_conductivity = Mi
		self._cell_model = cell_model
		self._stimulus = stimulus

	def conductivity(self):
		return self._intracellular_conductivity
	def time(self):
		return self._time
	def domain(self):
		return self._domain
	def cell_model(self):
		return self._cell_model
