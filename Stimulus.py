"""
Class to handle stimulus

"""
from dolfin import Measure

class Stimulus():
	def __init__(self,objects,keys,markers):
		self._objects = dict(zip(keys,objects))
		self._markers = markers

	def keys(self):
		return self._objects.keys()

	def values(self):
		return self._objects.values()

	def rhs(self,mesh,v):
		markers = self._markers
		dz = Measure("dx", domain=mesh, subdomain_data=markers)
		rhs = sum([g*v*dz(i) for (i, g) in zip(self.keys(), self.values())])
				
		return (dz, rhs)	

