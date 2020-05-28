"""
PDE solver for the monodomian 
equation

"""

from dolfin import *
from Utilities import TimeStepper

class PDESolver():

	def __init__(self,domain, time,Mi , Is=None, v_=None,params=None):
		# inptu store
		self._domain = domain 
		self._Mi = Mi 
		self._time = time 
		self._Is = Is

		if params != None:
			self.params = params 
		else:
			self.params = self.default_parameters

		# set up function space
		k = self.params["polynomial_degree"]
		V = FunctionSpace(self._domain,"CG",k) 
		self.V = V

		if v_ == None:
			self.v_ = Function(self.V,name='v_')
		else:
			self.v_ = v_
		self.v = Function(V,name="v")

	def step(self,interval):
		"""
		solution at one step

		"""
		(t0,t1) = interval
		dt = (t1-t0)
		Mi = self._Mi 
		theta = self.params["theta"]

		# set time
		t = t0 + theta*(t1-t0)
		self._time.assign(t)

		#variational formulation
		v = TrialFunction(self.V)
		w = TestFunction(self.V)
		v_mid = theta*v + (1.0-theta)*self.v_

		(dz,rhs) = (dx,0.0) if self._Is==None \
		 	else self._Is.rhs(self._domain,w)
		G = ((v-self.v_)/Constant(dt))*w*dz() \
			+ inner(Mi*grad(v_mid),grad(w))*dz() \
			- rhs

		a,L = system(G)
		pde = LinearVariationalProblem(a,L,self.v)

		solver = LinearVariationalSolver(pde)
		solver.solve()

	def solution_fields(self):

		return self.v_,self.v 

	@staticmethod
	def default_parameters():
		params = Parameters("PDESolver")
		params.add("theta",0.5)
		params.add("polynomial_degree",1)

		return params

	def solve(self,interval,dt):
		"""
		solves the problem
		"""
			
		time_stepper = TimeStepper(interval,dt)

		for t0,t1 in time_stepper:
			self.step((t0,t1))
			yield (t0,t1), self.solution_fields()
			# update previous solution
			self.v_.assign(self.v)