"""
heart EP
base class for operator splitting solver
of the EP problems

"""

class SplittingSolver():
	def __init__ (self,model,params=None):
		"""
		model = class of the heart
			with the info

		"""
		#initialize by taking the model (heart) info
		self._model = model
		self._domain = model.domain
		self._time = model.time #Do not know how this works

		if params = None
			self._parameters = self.default_parameters() 
		else
			self._parameters = params

		# Create ODE solver and extract solution fields
        self.ode_solver = self._create_ode_solver()
        (self.vs_, self.vs) = self.ode_solver.solution_fields()
        self.VS = self.vs.function_space()

        # Create PDE solver and extract solution fields
        self.pde_solver = self._create_pde_solver()
        (self.v_, self.v) = self.pde_solver.solution_fields()


	def _create_ode_solver(self):
		"""
		creates cell solver
		"""
		stimulus = self._model.stimulus

		cell_model = self._model.cell_model

		# Extract ode solver parameters
        params = self.parameters["ODESolver"]

        solver = ODESolver(self._domain, self._time, cell_model,
                                       I_s=stimulus,
                                       params=params)

        return solver


	def _create_pde_solver(self):
		"""
		creates diffusion solver (I guess?)
		"""

		Mi = self._model.conductivity
		solver = PDESolver(self._domain,self._time,Mi,v_=self.vs[0])

		return solver

	def solve(self):
		"""
		solves the problem
		"""


	def step(self):
		"""
		very important: manages the
		steps of the OP method.
		"""
	def merge(self,solution):
		"""
		merge solution from ode solver
		and pde solver
		
		"""
	def default_parameters():

	def solution_fields():

		return self.vs_, self.vs, self.v

