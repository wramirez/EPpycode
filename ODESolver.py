
"""
ODE solver 

"""
from dolfin import *
from Stimulus import Stimulus
from Utilities import TimeStepper,state_space,splat
import numpy as np
import os.path

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

		# Create (mixed) function space for potential + states
		v_family = self.params["V_polynomial_family"]
		v_degree = self.params["V_polynomial_degree"]
		s_family = self.params["S_polynomial_family"]
		s_degree = self.params["S_polynomial_degree"]

        #if functions spaces are the same 
		if (v_family == s_family and s_degree == v_degree):
			self.VS = VectorFunctionSpace(self._domain,v_family,v_degree,
								dim=self._num_states+1)
		else: 
			V = FunctionSpace(self._domain, v_family, v_degree)
			S = state_space(self._domain, self._num_states, s_family, s_degree)
			Mx = MixedElement(V.ufl_element(), S.ufl_element())
			self.VS = FunctionSpace(self._domain, Mx)

		self.vs_ = Function(self.VS,name="vs_")
		self.vs = Function(self.VS,name="vs")

		if os.path.exists("outputs_newton_raphson.out"):
			os.remove("outputs_newton_raphson.out")
			self.nrfile = open("outputs_newton_raphson.out","a+")
			self.nrfile.write("time,  ninter,  fnorm,  unorm, convergence\n")
		else:
			self.nrfile = open("outputs_newton_raphson.out","a+")
			self.nrfile.write("time,  ninter,  fnorm,  unorm, convergence\n")

	def time(self):
		return self._time

	def step(self,interval):
		"""

		"""

		(t0,t1) = interval
		dt = Constant(t1 - t0)

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

		F_theta = self._F(v_mid,s_mid,time=self.time())
		I_theta = -self._Iion(v_mid,s_mid,time=self.time())


		(dz,rhs) = self._Is.rhs(self._domain,w) if self._Is !=None else \
			(dx,0.0)

		lhs = (((v-v_)/dt)-I_theta)*w*dz() \
			+ inner(((s-s_)/dt)-F_theta,r)*dz()

		# set up system of equations
		G = lhs-rhs
		Jac = derivative(G,self.vs)
		# solve the system
		if False:
			pde = NonlinearVariationalProblem(G, self.vs, J=Jac)
			solver = NonlinearVariationalSolver(pde)
			solver_params = self.params["nonlinear_variational_solver"]
			solver.parameters.update(solver_params)
			solver.solve()
		else:
			a_tol, r_tol = 1e-7, 1e-10
			U_inc = Function(self.VS)
			nIter = 0
			eps = 1

			while eps > 1e-10 and nIter <= 20:              # Newton iterations
				nIter += 1
				A, b = assemble_system(Jac, -G)
				solve(A, U_inc.vector(), b)     # Determine step direction
				eps = np.linalg.norm(U_inc.vector().get_local(), ord = 2)

				a = assemble(G)
				fnorm = np.linalg.norm(a.get_local(), ord = 2)
				lmbda = 1.0     # step size, initially 1

				self.vs.vector()[:] += lmbda*U_inc.vector()    # New u vector

				print ('      {0:2d}  {1:3.2E}  {2:5e}'.format(nIter, eps, fnorm))
			
			if(nIter >= 20):
				self.nrfile.write("%.2f, %i, %.2e, %.2e, False\n"%(float(self.time()),nIter,fnorm,eps))
			else:
				self.nrfile.write("%.2f, %i, %.2e, %.2e, True\n"%(float(self.time()),nIter,fnorm,eps))
			self.nrfile.flush()
	def solution_fields(self):
		
		return self.vs_,self.vs 

	@staticmethod
	def default_parameters():
		
		params = Parameters("ODESolver")
		params.add("theta",0.5)
		params.add("V_polynomial_degree", 0)
		params.add("V_polynomial_family", "DG")
		params.add("S_polynomial_degree", 0)
		params.add("S_polynomial_family", "DG")

		# Use iterative solver as default.
		params.add(NonlinearVariationalSolver.default_parameters())
		params["nonlinear_variational_solver"]["newton_solver"]["linear_solver"] = "gmres"


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
			self.vs_.assign(self.vs)


class SingleCellSolver(ODESolver):
	def __init__(self,model,time,I_s,params=None):
		self._model=model
		mesh = UnitIntervalMesh(1)
		markers = MeshFunction("size_t",mesh,1)
		markers.set_all(1)
		stimulus = Stimulus((I_s,),(1,),markers)
		ODESolver.__init__(self,mesh,time,model,
				I_s=stimulus,params=params)