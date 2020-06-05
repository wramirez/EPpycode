"""
heart EP
base class for operator splitting solver
of the EP problems

"""
from dolfin import *
from ODESolver import ODESolver, PointODESolver
from PDESolver import PDESolver
from Utilities import TimeStepper
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
		

		cell_model = self._model.cell_model()

		# Extract ode solver parameters
		params = self._parameters["ODESolver"]

		solver = ODESolver(self._domain, self._time, cell_model,
		                               I_s=None,
		                               params=params)

		return solver


	def _create_pde_solver(self):
		"""
		creates diffusion solver
		"""
		stimulus = self._model.stimulus()

		params = self._parameters["PDESolver"]
		Mi = self._model.conductivity()
		solver = PDESolver(self._domain,self._time,Mi,
				Is=stimulus,v_=self.vs[0],params=params)

		return solver

	def solve(self,interval,dt):
		"""
		solves the problem
		"""

		time_stepper = TimeStepper(interval,dt)

		for t0,t1 in time_stepper:
			self.step((t0,t1))
			yield (t0,t1), self.solution_fields()
			# update previous solution
			self.vs_.assign(self.vs)
			

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
		# solve step 2 of the splitting method
		self.pde_solver.step((t0,t1))
		# update variables
		self.vs_.assign(self.vs)
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
		self.merger.assign(solution.sub(0),self.v)

	@staticmethod
	def default_parameters():

		params = Parameters("SplittingSolver")
		params.add("theta",0.5)
		ode_solver_parameters = ODESolver.default_parameters()
		ode_solver_parameters["V_polynomial_degree"] = 1
		ode_solver_parameters["V_polynomial_family"] = "CG"
		pde_solver_parameters = PDESolver.default_parameters()
		params.add(ode_solver_parameters)
		params.add(pde_solver_parameters)

		return params

	def solution_fields(self):

		return self.vs_, self.vs, self.v


class PointSplittingSolver(SplittingSolver):
		def __init__(self, model, params= None):
			SplittingSolver.__init__(self,model,params)

		def _create_ode_solver(self):
			"""
			creates cell solver
			"""

			cell_model = self._model.cell_model()

			# Extract ode solver parameters
			Solver = eval(self._parameters["ode_solver_choice"])
			params = self._parameters[Solver.__name__]
			solver = Solver(self._domain, self._time, cell_model,
											I_s=None,params=params)

			return solver
		@staticmethod
		def default_parameters():
			params = Parameters("SplittingSolver")
			params.add("theta",0.5,0,1)
			try:
				params.add("ode_solver_choice","PointODESolver",set(["PointODESolver","ODESolver"]))
			except:
				params.add("ode_solver_choice","PointODESolver",["PointODESolver","ODESolver"])
				pass
			point_ode_solver_params = PointODESolver.default_parameters()
			point_ode_solver_params["scheme"] = "RL1"
			params.add(point_ode_solver_params)

			ode_solver_params = ODESolver.default_parameters()
			ode_solver_params["V_polynomial_degree"] =1
			ode_solver_params["V_polynomial_family"] = "CG"
			params.add(ode_solver_params)

			pde_solver_parameters = PDESolver.default_parameters()
			params.add(pde_solver_parameters)

			return params

