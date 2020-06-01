"""
Test for Splitting solver 

"""

from dolfin import *
from SplittingSolver import SplittingSolver
from Stimulus import Stimulus
from CardiacModel import CardiacModel
import numpy as np
from cell_models.Fenton_Karma_BR_altered import Fenton_Karma_1998_BR_altered
import matplotlib.pyplot as plt 

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

class MyExpression(UserExpression):
	def __init__(self,period,duration,start,Iamp,w,t,**kwargs):
		self.period = period
		self.duration = duration
		self.start = start
		self.Iamp = Iamp
		self.t = t
		self.w = w
		super().__init__(self,**kwargs)
	def eval(self,value,x):
		t = float(self.t)
		w = self.w
		# duration = self.duration
		# #
		# f = 0.2
		# c = 0.7
		# Iamp = self.Iamp*np.exp(-(((t-self.start)**2)/(f*c**2)))
		Iamp = self.Iamp
		I =(Iamp if t-int(t/self.period)*self.period  <= self.duration+self.start else 0.0) \
			 if (t-int(t/self.period)*self.period >= self.start) else 0.0*Iamp
		
		
		value[0] = I

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
	fiber_vals[:] = 1.0/np.sqrt(3.0)
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
		s_t = 0.1
		s_l = 0.1

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
	wstim = 0.5
	mesh,fibers,stimcells = create_mesh(0.1,10.0,0.1,0.1,wstim=wstim)
	M = conductivity(fibers)

	# Surface to volume ratio
	chi = 140.0     # mm^{-1}
	# Membrane capacitance
	C_m = 0.01 # mu F / mm^2
	duration = 2. # ms
	A = 60000. # mu A/cm^3
	cm2mm = 10.
	factor = 1.0/(chi*C_m) # NB: cbcbeat convention
	amplitude = factor*A*(1./cm2mm)**3 # mV/ms
	print ("amplitude: ",amplitude)
	start = 2.0
	period = 300.0

	I_s = MyExpression(period,duration,start,amplitude,wstim,time)
		
	cell_model = Fenton_Karma_1998_BR_altered()
	stim = Stimulus((I_s,),(1,),stimcells)
	heart = CardiacModel(mesh,time,M,cell_model,stim)

	return heart,cell_model,I_s,mesh

def solve_cardiac_model():
	cardiac_model,cell_model,I_s,mesh = setup_cardiac_model()
	solver = SplittingSolver(cardiac_model)
	(vs_,vs,v) = solver.solution_fields()
	vs_.assign(interpolate(cell_model.initial_conditions(),vs_.function_space()))
	
	dt = 0.05
	interval = (0.0,200.0)

	pvdfile = File("outputs/voltage.pvd")
	printstep = 20.0
	count = 0

	isfile = File("outputs/I_s.pvd")
	sfile = File("outputs/s_0.pvd")
	V = FunctionSpace(mesh,'CG',1)

	for (timestep,fields) in solver.solve(interval,dt):

		(vs_,vs,v) =fields
		(t0,t) = timestep
		(volt,s0,s1) = vs.split(True) 
		if float(count)%printstep == 0:
			pvdfile << (volt,t0)
			sfile << (s1,t0)
			Isp = project(I_s,V)
			isfile << (Isp,t0) 
		count += 1
		
		print("(t_0, t_1) = (%g, %g),dt:%g"%(t0,t,t-t0))
solve_cardiac_model()




