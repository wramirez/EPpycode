"""
Test for Splitting solver 

"""

from dolfin import *
from SplittingSolver import SplittingSolver
from Stimulus import Stimulus
from CardiacModel import CardiacModel
import numpy as np
from cell_models.Fenton_Karma_BR_altered import Fenton_Karma_1998_BR_altered

class MyExpression():
	def __init__(self,period,duration,start,Iamp,t):
		self.period = period
		self.duration = duration
		self.start = start
		self.Iamp = Iamp
		self.t = t

	def eval(self,value,x):
		t = self.t
		I =(self.Iamp if t-int(t/self.period)*self.period  <= self.duration+self.start else 0.0) \
			if (t-int(t/self.period)*self.period >= self.start) else 0.0*self.Iamp

def create_mesh(dx,Lx,Ly,Lz,wstim):
	N = lambda v:int(np.int(v))
	mesh = Mesh()
	mesh = BoxMesh(mesh.mpi_comm(),
		Point(0.0, 0.0, 0.0),
		Point(Lx, Ly, Lz),
		N(Lx/dx), N(Ly/dx), N(Lz/dx))

	stimcells = MeshFunction("size_t",mesh,3)
	stimcells.set_all(0)

	V = FunctionSpace(mesh,'DG',0)
	ncells = mesh.num_cells()
	W = VectorFunctionSpace(mesh,'DG',0,3)
	fiber = Function(W)
	fiber_vals = np.zeros(ncells*3)
	fiber_vals[::3] = 1.0
	fiber.vector()[:] = fiber_vals

	for cell in cells(mesh):
		xx = cell.midpoint().x()
		yy = cell.midpoint().y()
		zz = cell.midpoint().z()
		if (xx <= wstim):
			stimcells[cell]= 1

	return mesh,fiber,stimcells

def conductivity(fiber,s_t=None,s_l=None):
	if s_t == None and s_l == None:
		s_t = 1.0
		s_l = 1.0

	t11 = (s_t+(s_l-s_t)*fiber[0]*fiber[0])
	t22 = (s_t+(s_l-s_t)*fiber[1]*fiber[1])
	t33 = (s_t+(s_l-s_t)*fiber[2]*fiber[2])
	t12 = ((s_l-s_t)*fiber[0]*fiber[1])
	t23 = ((s_l-s_t)*fiber[1]*fiber[2])
	t31 = ((s_l-s_t)*fiber[2]*fiber[0])

	M = as_tensor(((t11,t12,t31),(t12,t22,t23),(t31,t23,t33)))

	return M

def setup_cardiac_model():
	time = Constant(0.0)
	mesh,fibers,stimcells = create_mesh(1.0,5.0,5.0,1.0,wstim=1.0)
	M = conductivity(fibers)

	# Surface to volume ratio
	chi = 140.0     # mm^{-1}
	# Membrane capacitance
	C_m = 0.01 # mu F / mm^2
	duration = 2. # ms
	A = 50000. # mu A/cm^3
	cm2mm = 10.
	factor = 1.0/(chi*C_m) # NB: cbcbeat convention
	amplitude = factor*A*(1./cm2mm)**3 # mV/ms
	start = 2.0
	period = 300.0

	I_s = MyExpression(period,duration,start,amplitude,time)

	
	cell_model = Fenton_Karma_1998_BR_altered()
	stim = Stimulus((I_s,),(1,),stimcells)
	heart = CardiacModel(mesh,time,M,cell_model,stim)

	return heart

def solve_cardiac_model():
	cardiac_model = setup_cardiac_model()
	solver = SplittingSolver(cardiac_model)
	dt = 0.1
	interval = (0.0,10.0)

	solver.solve(interval,dt)
	# ts = solver.time_stepper(interval,dt)
	# print(ts)
	# for (timestep,fields) in solver.solve(interval,dt):
		# print("!hola!!")
		# print("(t_0, t_1) = (%g, %g)", timestep)

solve_cardiac_model()




