"""
Test for Splitting solver 

"""

from dolfin import *
import SplittingSolver
import Stimulus
import CardiacModel


def create_mesh(dx,lx,ly,lz,wstim):
	N = lambda v:int(np.int(v))
	mesh = BoxMesh(mpi_comm_world(),
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

def setup_cardiac_model(cell_model,mesh)
	time = Constant(0.0)
	
	M = conductivity(fiber)





