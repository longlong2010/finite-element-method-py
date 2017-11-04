from enum import Enum;
from geometry import *;

class Load(Enum):
	X = Dof.X;
	Y = Dof.Y;
	Z = Dof.Z;

	def __init__(self, dof):
		self.dof = dof;
	
	def getDof(self):
		return self.dof;

class Constraint(Enum):
	X = Dof.X;
	Y = Dof.Y;
	Z = Dof.Z;

	def __init__(self, dof):
		self.dof = dof;

	def getDof():
		return self.dof;
