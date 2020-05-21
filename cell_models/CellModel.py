"""
Base class for cell models

"""

from collections import OrderedDict
from dolfin import UserExpression

class initial_conditions_expression(UserExpression):
	def __init__(self,initial_conditions,num_states,**kwargs):
		self.initial_conditions = initial_conditions
		self.num_states = num_states
		super().__init__(self,**kwargs)
	def eval(self,value,x):
		for i,(key,values) in enumerate(self.initial_conditions.items()):
			value[i] = values
	def value_shape(self):
		return (self.num_states+1,)

class CellModel(object):
	def __init__(self,params=None,initial_conditions=None):

		self._params = self.default_parameters()
		self._initial_conditions = self.default_initial_conditions()

	@staticmethod
	def default_parameters():
		return OrderedDict()
	
	@staticmethod
	def default_initial_conditions():
		return OrderedDict()
		
	def set_parameters(self,**params):
		for param_name,param_value in params.items():
			self._params[param_name] = param_value

	def set_initial_conditions(self,**init):
		for init_name,init_value in init.items():
			self._initial_conditions[init_name] = init_value

	def initial_conditions(self):

		"""
		This function returns initial conditions
		as an UserExpression
		"""
		return initial_conditions_expression(self._initial_conditions,
				self.num_states(),degree=1)

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

	def __dir__(self):
		return ['default_parameters','defaul_initial_conditions',
			'F','I','num_states','set_parameters','set_initial_conditions',
			'parameters','initial_conditions']