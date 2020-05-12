
"""
ODE solver 

"""
from dolfin import *

class ODESolver():
	def __init__(self,domain,time,model,I_s=None,params=None):
		# store input
		self._domain = domain
		self._time = time 
		self._model = model #cell model

		# Info from cell model
		self._F = self._model.F 
		self._Iion = self._model.I 
		self._num_states = self._model.num_states()

		# Stimulus
		self._Is = Stimulus(Is)

		if params != None:
			self.params = params 
		else:
			self.params = self.default_parameters()

		# create function space
		self.VS = VectorFunctionSpace(domain,"CG",1,
								dim=self._num_states+1)

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
		v_ = self.vs[0]
		s_ = as_vector([self.vs[i] for i in range(1,dim)])

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


		(dz,rhs) = rhs_with_stimulus(self._Is,self._domain,w)

		lhs = ((v-v_)/dt-Itheta)*w*dz + inner((s-s_)/dt-F_theta,r)*dz

		# set up system of equations
		G = lhs,rhs


		# solve the system
		pde = NonlinearVariationalProblem(G, self.vs, J=derivative(G, self.vs))
        solver = NonlinearVariationalSolver(pde)
        solver.solve()

	def solution_fields(self):
		
		return self.vs_,self.vs 
		
	def default_parameters(self):
		pass
	def solve():
		pass