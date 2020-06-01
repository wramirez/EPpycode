"""
test Single Cell

"""

import matplotlib.pyplot as plt
from dolfin import *
from ODESolver import SingleCellSolver
from cell_models.Fenton_Karma_BR_altered import Fenton_Karma_1998_BR_altered

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
period = 200.0

I_s = MyExpression(period,duration,start,amplitude,time)


cell_model = Fenton_Karma_1998_BR_altered()
cell_solver = SingleCellSolver(cell_model,time,I_s)

dt = 0.1
interval = (0.0,400.0)
voltage = []
time = []
printstep = 10.0
count = 0

for (timestep,fields) in cell_solver.solve(interval,dt):

		(vs_,vs) =fields
		(t0,t) = timestep
		if float(count)%printstep == 0:
			voltage.append(vs.vector().get_local()[2])
			time.append(t)
		count += 1
		
		print("(t_0, t_1) = (%g, %g)"%timestep)

plt.figure()
plt.plot(time,voltage,lw=2)
plt.show()