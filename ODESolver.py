
"""
ODE solver 

"""
from dolfin import *

class ODESolver(object):
	def __init__(self,domain,time,model,I_s,params=None):
		# store input
		self._domain = domain
		self._time = time 
		self._model = model #cell model

		self._F = self._model.F
		self._Iion = self._model.I 
		self._num_states = self._model.num_states()

		# Stimulus
		self._Is = I_s

		if params != None:
			self.params = params 
		else:
			self.params = self.default_parameters()

		# create function space
		self.VS = VectorFunctionSpace(domain,"CG",1,
								self._num_states+1)

		self.vs_ = Function(self.VS,name="vs_")
		self.vs = Function(self.VS,name="vs")

	def time(self):
		return self._time

	def step(self,interval):
		"""

		"""

		(t0,t1) = interval
		dt = t1 - t0

		# test function
		wtest = TestFunction(self.VS)

		# extract previous solution
		dim = self._num_states+1
		v_ = self.vs_[0]
		s_ = as_vector([self.vs_[i] for i in range(1,dim)])

		# current variables
		self.vs.assign(self.vs_)
		v = self.vs[0]
		s = as_vector([self.vs[i] for i in range(1,dim)])

		# test functions for voltage and internal variables
		w = wtest[0]
		r = as_vector([wtest[i] for i in range(1,dim)])

		theta = self.params["theta"]

		# time for time dependen variables
		t = t0 + theta*dt 
		self._time.assign(t)

		# designate v_mid and s_mid for the evaluation of F(v,s) and 
		# I(v,s)
		v_mid = theta*v + (1.0-theta)*v_ 
		s_mid = theta*s + (1.0-theta)*s_

		F_theta = self._F(v_mid,s_mid,time=self._time)
		I_theta = - self._Iion(v_mid,s_mid,time=self._time)


		(dz,rhs) = self._Is.rhs(self._domain,w)

		lhs = ((v-v_)/dt-I_theta)*w*dz() + inner((s-s_)/dt-F_theta,r)*dz()

		# set up system of equations
		G = lhs-rhs

		# solve the system
		pde = NonlinearVariationalProblem(G, self.vs, J=derivative(G, self.vs))
		solver = NonlinearVariationalSolver(pde)
		solver_params = self.params["nonlinear_variational_solver"]
		solver.parameters.update(solver_params)
		solver.solve()

	def solution_fields(self):
		
		return self.vs_,self.vs 

	@staticmethod
	def default_parameters():
		
		params = Parameters("ODESolver")
		params.add("theta",0.5)

		# Use iterative solver as default.
		params.add(NonlinearVariationalSolver.default_parameters())
		params["nonlinear_variational_solver"]["newton_solver"]["linear_solver"] = "gmres"


		return params

	def solve():
		pass