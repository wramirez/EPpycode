"""
heart EP
base class for operator splitting solver
of the EP problems

"""
from dolfin import *
from ODESolver import ODESolver
from PDESolver import PDESolver
import numpy as np

class SplittingSolver:
	def __init__ (self,model,params=None):
		"""
		model = class of the heart
			with the info

		"""
		#initialize by taking the model (heart) info
		self._model = model
		self._domain = model.domain()
		self._time = model.time() 

		if params == None:
			self._parameters = self.default_parameters() 
		else:
			self._parameters = params
		

		# Create ODE solver and extract solution fields
		self.ode_solver = self._create_ode_solver()
		(self.vs_, self.vs) = self.ode_solver.solution_fields()
		self.VS = self.vs.function_space()

		# Create PDE solver and extract solution fields
		self.pde_solver = self._create_pde_solver()
		(self.v_, self.v) = self.pde_solver.solution_fields()

		V = self.v.function_space()
		self.merger = FunctionAssigner(self.VS.sub(0),V)

		
        
	def _create_ode_solver(self):
		"""
		creates cell solver
		"""
		stimulus = self._model.stimulus()

		cell_model = self._model.cell_model()

		# Extract ode solver parameters
		params = self._parameters["ODESolver"]

		solver = ODESolver(self._domain, self._time, cell_model,
		                               I_s=stimulus,
		                               params=params)

		return solver


	def _create_pde_solver(self):
		"""
		creates diffusion solver
		"""

		Mi = self._model.conductivity
		solver = PDESolver(self._domain,self._time,Mi,v_=self.vs[0])

		return solver
	@staticmethod
	def time_stepper(interval,dt):
		(t0,t1) = interval
		t =[(t0,t0+dt)]
		ind = 0
		
		while t[ind][1]<=t1:
			t.append((t[ind][1],t[ind][1]+dt))
		return t 

	def solve(self,interval,dt):
		"""
		solves the problem
		"""

		# implement time stepper
		# return an iterable	
		# time_stepper = self.time_stepper(interval,dt)
		
		t = np.arange(interval[0],interval[1],dt)

		for ti in t:
			t0 = ti 
			t1 = ti + dt
			self.step((t0,t1))
			# yield (t0,t1), self.solution_fields()
			print (t0,t1)
			# update previous solution
			self.vs_.assing(self.vs)

	def step(self,interval):
		"""
		very important: manages the
		steps of the OP method.
		"""
		#parameter theta
		theta = self._parameters["theta"]

		#time domain
		(t0,t1) = interval
		dt = (t1-t0)
		t = t0+theta*dt 

		# solve step 1 of the splitting method
		self.ode_solver.step((t0,t))
		# solve step 2 of the splitiing method
		self.pde_solver.step((t0,t1))
		# update variables
		self.vs_.assing(self.vs)
		# set the voltage part of vs to the actual voltage given 
		# by the solution of the pde
		self.merge(self.vs_) 

		# solve last step (step 3) 
		self.ode_solver.step((t,t1))

	def merge(self,solution):
		"""
		merge solution from ode solver
		and pde solver
		
		"""
		self.merger.assing(solution.sub(0),self.v)

	@staticmethod
	def default_parameters():

		params = Parameters("SplittingSolver")
		params.add("theta",0.5)
		ode_solver_parameters = ODESolver.default_parameters()
		pde_solver_parameters = PDESolver.default_parameters()
		params.add(ode_solver_parameters)
		params.add(pde_solver_parameters)

		return params

	def solution_fields(self):

		return self.vs_, self.vs, self.v

