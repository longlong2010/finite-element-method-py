import numpy;
from enum import Enum;
from load import *;

Dof = Enum('Dof', 'X Y Z');

class Node:
	def __init__(self, x, y, z):
		self.x = x;
		self.y = y;
		self.z = z;

		self.dofs = set();
		self.loads = dict();
		self.constraints = dict();
		self.values = dict();

	def addLoad(self, load, v):
		if self.loads.__contains__(load):
			self.loads[load] += v;
		else:
			self.loads[load] = v;
		
	def addConstraint(self, constraint, v = 0):
		self.constraints[constraint] = v;
	
	def getLoads(self):
		return self.loads;

	def getConstraints():
		return self.constraints;

	def getDofs():
		return self.dofs;

	def getDofNum():
		return len(self.dofs);

class StructuralNode(Node):
	def __init__(self, x, y, z):
		super(StructuralNode, self).__init__(x, y, z);
		self.dofs.add(Dof.X);
		self.dofs.add(Dof.Y);
		self.dofs.add(Dof.Z);

class Element:
	def __init__(self):
		self.nodes = [];
	
	def getNodeNum(self):
		return len(self.nodes);

	def getDofNum(self):
		ndof = 0;
		for node in self.nodes:
			ndof += node.getDofNum();
		return ndof;
