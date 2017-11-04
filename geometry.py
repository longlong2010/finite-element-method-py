import abc;
import numpy;
from enum import Enum;
from load import *;
from property import *;

Dof = Enum('Dof', 'X Y Z');

class Node:
	def __init__(self, x, y, z):
		self.x = x;
		self.y = y;
		self.z = z;

		self.dofs = set();
		self.dofs.add(Dof.X);
		self.dofs.add(Dof.Y);
		self.dofs.add(Dof.Z);
		
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

	def getConstraints(self):
		return self.constraints;

	def getDofs(self):
		return self.dofs;

	def getDofNum(self):
		return len(self.dofs);

class Element(metaclass = abc.ABCMeta):
	def __init__(self, nodes, material, points):
		self.nodes = nodes;
		self.material = material;
		self.points = points;
	
	def getNodeNum(self):
		return len(self.nodes);

	def getDofNum(self):
		ndof = 0;
		for node in self.nodes:
			ndof += node.getDofNum();
		return ndof;

	def getCoord(self):
		nnode = self.getNodeNum();
		coord = numpy.zeros((nnode, 3));
		i = 0;
		for n in self.nodes:
			coord[i] = numpy.array([n.x, n.y, n.z]);
			i += 1;
		return coord;

	@abc.abstractmethod	
	def getShapeDer(self, p):
		pass;
	
	@abc.abstractmethod
	def getShapeFunc(self, p):
		pass;

	def getJacobi(self, p):
		der = self.getShapeDer(p);
		coord = self.getCoord();
		return der.dot(coord);
	
	def getStiffMatrix(self):
		ndof = self.getDofNum();
		Ke = numpy.zeros((ndof, ndof));
		for p in self.points:
			N = self.getShapeFunc(p);
			J = self.getJacobi(p);
		return Ke;
	
	def getMassMatrix(self):
		ndof = self.getDofNum();
		Me = numpy.zeros((ndof, ndof));
		rho = self.material.getProperty(MaterialProperty.rho);
		for p in self.points:
			N = self.getShapeFunc(p);
			J = self.getJacobi(p);
			Me += N.T.dot(N) * numpy.linalg.det(J);
		return Me;

class Tet4Element(Element):
	def __init__(self, n1, n2, n3, n4, material):
		super(Tet4Element, self).__init__([n1, n2, n3, n4], material, [[0.25, 0.25, 0.25, 0.25, 1]]);

	def getShapeFunc(self, p):
		[l1, l2, l3, l4, w] = p;
		return numpy.array([
			[l1],
			[l2],
			[l3],
			[l4]
		]);

	def getShapeDer(self, p):
		return numpy.array([
			[1, 0, 0, -1],
			[0, 1, 0, -1],
			[0, 0, 1, -1]
		]);
