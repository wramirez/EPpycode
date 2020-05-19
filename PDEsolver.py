"""
PDE solver for the monodomian 
equation

"""

from dolfin import *

class PDESolver():

	def __init__(self,domain, time, Mi, v_=None,params=None):
		# inptu store
		self._domain = domain 
		self._Mi = Mi 
		self._time = time 

		if params != None:
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

		# set time
		t = t0 + theta*(t1-t0)
		self._time.assing(t)

		#variational formulation
		u = TrialFunction(self.V)
		w = TestFunction(self.V)
		v_mid = theta*v + (1.0-theta)*self.v_
		G = ((v-self.v_)/Constant(dt))*w*dx() \
			+ inner(Mi*nabla_grad(v_mid),nabla_grad(w))*dx()

		a,L = system(G)
		pde = LinearVariationalProblem(a,L,self.v)

		solver = LinearVariationalSolver(pde)
		solver.solve()

	def solution_fields(self):

		return self.v_,self.v 

	@staticmethod
	def default_parameters():
		params = Parameters("MonodomainSolver")
		params.add("theta",0.5)

		return params

	def solve():
		"""
		I don't if this is needed at this 
		moment
		"""
		pass