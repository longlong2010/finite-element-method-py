import abc;
import numpy;
from enum import Enum;
Dof = Enum('Dof', 'X Y Z');
from load import *;
from property import *;

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

	def setValue(self, dof, v):
		if dof in self.dofs:
			self.values[dof] = v;

	def getConstraints(self):
		return self.constraints;

	def getDofs(self):
		dofs = [];
		for dof in Dof:
			if dof in self.dofs:
				dofs.append(dof);
		return dofs;

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
	
	def getNodes(self):
		return self.nodes;

	def getCoord(self):
		nnode = self.getNodeNum();
		coord = numpy.zeros((nnode, 3));
		i = 0;
		for n in self.nodes:
			coord[i] = numpy.array([n.x, n.y, n.z]);
			i += 1;
		return coord;

	@abc.abstractmethod	
	def getShapeDerMatrix(self, p):
		pass;
	
	@abc.abstractmethod
	def getShapeMatrix(self, p):
		pass;
	
	@abc.abstractmethod
	def getStrainMatrix(self, p):
		pass;

	def getJacobi(self, p):
		der = self.getShapeDerMatrix(p);
		coord = self.getCoord();
		return der.dot(coord);
	
	def getStiffMatrix(self):
		ndof = self.getDofNum();
		Ke = numpy.zeros((ndof, ndof));
		D = self.getStressStrainMatrix();
		for p in self.points:
			B = self.getStrainMatrix(p);
			J = self.getJacobi(p);
			Ke += B.T.dot(D).dot(B) * abs(numpy.linalg.det(J)) * p[4];
		return Ke;
	
	def getMassMatrix(self):
		ndof = self.getDofNum();
		Me = numpy.zeros((ndof, ndof));
		rho = self.material.getProperty(MaterialProperty.rho);
		for p in self.points:
			N = self.getShapeMatrix(p);
			J = self.getJacobi(p);
			Me += rho * N.T.dot(N) * abs(numpy.linalg.det(J)) * p[4];
		return Me;

	def getStressStrainMatrix(self):
		D = numpy.zeros((6, 6));
		E = self.material.getProperty(MaterialProperty.E);
		nu = self.material.getProperty(MaterialProperty.nu);
		D[0][0] = 1 - nu;
		D[0][1] = nu;
		D[0][2] = nu;
		
		D[1][0] = nu;
		D[1][1] = 1 - nu;
		D[1][2] = nu;
		
		D[2][0] = nu;
		D[2][1] = nu;
		D[2][2] = 1 - nu;

		D[3][3] = (1 - 2 * nu) / 2;
		D[4][4] = (1 - 2 * nu) / 2;
		D[5][5] = (1 - 2 * nu) / 2;

		D *= E / ((1 + nu) * (1 - 2 * nu));
		return D;


class Tet4Element(Element):
	def __init__(self, n1, n2, n3, n4, material):
		super(Tet4Element, self).__init__([n1, n2, n3, n4], material, [[0.25, 0.25, 0.25, 0.25, 1 / 6]]);

	def getShapeMatrix(self, p):
		[l1, l2, l3, l4, w] = p;
		ndof = self.getDofNum();
		N = numpy.zeros((3, ndof));
		for i in range(0, ndof):
			N[i % 3][i] = p[i % 3];
		return N;

	def getShapeDerMatrix(self, p):
		return numpy.array([
			[1, 0, 0, -1],
			[0, 1, 0, -1],
			[0, 0, 1, -1]
		]);

	def getStrainMatrix(self, p):
		der = self.getShapeDerMatrix(p);
		J = self.getJacobi(p);
		Der = numpy.linalg.inv(J).dot(der);
		ndof = self.getDofNum();
		nnode = self.getNodeNum();
		B = numpy.zeros((6, ndof));
		for i in range(0, nnode):
			k = i * 3;
			for j in range(0, 3):
				B[j][i * 3 + j] = Der[j][i];
			B[3][k] = Der[1][i];
			B[3][k + 1] = Der[0][i];
			
			B[4][k + 1] = Der[2][i];
			B[4][k + 2] = Der[1][i];
			
			B[5][k] = Der[2][i];
			B[5][k + 2] = Der[0][i];
		return B;
	
	def getStress(self):
		ndof = self.getDofNum();
		u = numpy.zeros((ndof, 1));
		k = 0;
		for n in self.nodes:
			dofs = n.getDofs();
			ndof = 0;
			for d in dofs:
				u[k + ndof][0] = n.values[d];
				ndof += 1;
			k += n.getDofNum();
		D = self.getStressStrainMatrix();
		B = self.getStrainMatrix(self.points[0]);
		return D.dot(B).dot(u);
