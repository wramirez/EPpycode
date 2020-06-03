"""
test Single Cell

"""

import matplotlib.pyplot as plt
from dolfin import *
from ODESolver import SingleCellSolver, PointSingleCellSolver
from cell_models.Fenton_Karma_BR_altered import Fenton_Karma_1998_BR_altered
import numpy as np

# Turn on FFC/FEniCS optimizations
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
flags = ["-O3", "-ffast-math", "-march=native"]
parameters["form_compiler"]["cpp_optimize_flags"] = " ".join(flags)
parameters["form_compiler"]["quadrature_degree"] = 3

class MyExpression(UserExpression):
	def __init__(self,period,duration,start,Iamp,t,**kwargs):
		self.period = period
		self.duration = duration
		self.start = start
		self.Iamp = Iamp
		self.t = t
		super().__init__(self,**kwargs)
	def eval(self,value,x):
		t = float(self.t)
		I =(self.Iamp if t-int(t/self.period)*self.period  <= self.duration+self.start else 0.0) \
			if (t-int(t/self.period)*self.period >= self.start) else 0.0*self.Iamp
		value[0] = I


time = Constant(0.0)

# Surface to volume ratio
chi = 140.0     # mm^{-1}
# Membrane capacitance
C_m = 0.01 # mu F / mm^2
duration = 5. # ms
A = 25000. # mu A/cm^3
cm2mm = 10.
factor = 1.0/(chi*C_m) # NB: cbcbeat convention
amplitude = factor*A*(1./cm2mm)**3 # mV/ms
start = 2.0
period = 400.0

I_s = MyExpression(period,duration,start,amplitude,time)


cell_model = Fenton_Karma_1998_BR_altered()
cell_solver = PointSingleCellSolver(cell_model,time,I_s)
(vs_,vs) = cell_solver.solution_fields()
vs_.assign(interpolate(cell_model.initial_conditions(),vs_.function_space()))


dt = 0.02
interval = (0.0,400.0)
voltage = []
s0 = []
s1 = []
time = []
printstep = 10.0
count = 0

# if True:
if False:

	for (timestep,fields) in cell_solver.solve(interval,dt):

			(vs_,vs) =fields
			(t0,t) = timestep
			if float(count)%printstep == 0:
				voltage.append(vs.vector().get_local()[0])
				s0.append(vs.vector().get_local()[1])
				s1.append(vs.vector().get_local()[2])
				time.append(t)
			count += 1
			np.save('outputs/voltage.npy',np.array(voltage))
			np.save('outputs/v.npy',np.array(s0))
			np.save('outputs/w.npy',np.array(s1))
			np.save('outputs/time.npy',np.array(time))
			print("(t_0, t_1) = (%g, %g)"%timestep)
else:
	voltage  = np.load('outputs/voltage.npy')
	s0  = np.load('outputs/v.npy')
	s1  = np.load('outputs/w.npy')
	time = np.load('outputs/time.npy')
	plt.figure()
	plt.plot(time,voltage,lw=2,label='V')
	# plt.plot(time,s0,lw=2,label='v')
	# plt.plot(time,s1,lw=2,label='w')
	plt.show()