import abc;
import numpy;
import math;
from enum import Enum;
Dof = Enum('Dof', 'X Y Z');
from load import *;
from property import *;

class Node:
	def __init__(self, nid, x, y, z):
		self.nid = nid;
		self.fid = 0;
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

	def getNid(self):
		return self.nid;

	def getFid(self):
		return self.fid;

	def setFid(self, fid):
		self.fid = fid;

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
	def __init__(self, eid, nodes, property, points):
		self.eid = eid;
		self.nodes = nodes;
		self.property = property;
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
	
	@abc.abstractmethod	
	def getShapeDerMatrix(self, p):
		pass;
	
	@abc.abstractmethod
	def getShapeMatrix(self, p):
		pass;

	@abc.abstractmethod
	def getStrainMatrix(self, p):
		pass;
	
	@abc.abstractmethod
	def getStiffMatrix(self, p):
		pass;

	@abc.abstractmethod
	def getStressStrainMatrix(self):
		pass;

class Element3D(Element, metaclass = abc.ABCMeta):
	def getCoord(self):
		nnode = self.getNodeNum();
		coord = numpy.zeros((nnode, 3));
		i = 0;
		for n in self.nodes:
			coord[i] = numpy.array([n.x, n.y, n.z]);
			i += 1;
		return coord;

	@abc.abstractmethod
	def addPload(self, load, v, n1, n2):
		pass;

	def getJacobi(self, p):
		der = self.getShapeDerMatrix(p);
		coord = self.getCoord();
		return der.dot(coord);
	
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
	
	def getStiffMatrix(self):
		ndof = self.getDofNum();
		Ke = numpy.zeros((ndof, ndof));
		D = self.getStressStrainMatrix();
		for p in self.points:
			B = self.getStrainMatrix(p);
			J = self.getJacobi(p);
			Ke += B.T.dot(D).dot(B) * abs(numpy.linalg.det(J)) * p[-1];
		return Ke;
	
	def getMassMatrix(self):
		ndof = self.getDofNum();
		Me = numpy.zeros((ndof, ndof));
		rho = self.property.material.getProperty(MaterialProperty.rho);
		for p in self.points:
			N = self.getShapeMatrix(p);
			J = self.getJacobi(p);
			Me += rho * N.T.dot(N) * abs(numpy.linalg.det(J)) * p[-1];
		return Me;

	def getStressStrainMatrix(self):
		D = numpy.zeros((6, 6));
		E = self.property.material.getProperty(MaterialProperty.E);
		nu = self.property.material.getProperty(MaterialProperty.nu);
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

class TetElement(Element3D, metaclass = abc.ABCMeta):
	@staticmethod
	def getTriangleArea(nodes):
		coord = numpy.zeros((len(nodes), 3));
		k = 0;
		for n in nodes:
			coord[k] = numpy.array([n.x, n.y, n.z]);
			k += 1;
		A = 0;
		M = numpy.ones((3, 3));
		for i in range(0, 3):
			j = (i + 1) % 3;
			for k in range(0, 3):
				M[0, k] = coord[k, i];
				M[1, k] = coord[k, j];
			A += numpy.linalg.det(M) ** 2;
		return 0.5 * math.sqrt(A);

class Tet4Element(TetElement):
	def __init__(self, eid, nodes, property):
		super(Tet4Element, self).__init__(eid, nodes, property, [[0.25, 0.25, 0.25, 0.25, 1 / 6]]);

	def getShapeMatrix(self, p):
		[l1, l2, l3, l4, w] = p;
		ndof = self.getDofNum();
		nnode = self.getNodeNum();
		N = numpy.zeros((3, ndof));
		for i in range(0, nnode):
			for j in range(0, 3):
				N[j][i * 3 + j] = p[i];
		return N;

	def getShapeDerMatrix(self, p):
		return numpy.array([
			[1, 0, 0, -1],
			[0, 1, 0, -1],
			[0, 0, 1, -1]
		]);
	
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

	def addPload(self, load, v, n1, n2):
		k = -1;
		numbers = [];
		nodes = [];
		for i in range(0, 4):
			n = self.nodes[i];
			if n == nid:
				k = i;
			else:
				nodes.append(n);
		assert(k >= 0);
		if k == 0:
			numbers += [1, 2, 3];
		elif k == 1:
			numbers += [0, 2, 3];
		elif k == 2:
			numbers += [0, 1, 3];
		elif k == 3:
			numbers += [0, 1, 2];

		A = TetElement.getTriangleArea(nodes);
		for no in numbers:
			n = self.nodes[no];
			dofs = n.getDofs();
			for dof in dofs:
				if dof == load.getDof():
					n.addLoad(load, v / len(numbers) * A);

class Tet10Element(TetElement):
	def __init__(self, eid, nodes, property):
		alpha = 0.58541020;
		beta = 0.13819660;
		super(Tet10Element, self).__init__(eid, nodes, property, [
				[alpha, beta, beta, beta, 1 / 24], 
				[beta, alpha, beta, beta, 1 / 24], 
				[beta, beta, alpha, beta, 1 / 24], 
				[beta, beta, beta, alpha, 1 / 24]]);

	def getShapeMatrix(self, p):
		[l1, l2, l3, l4, w] = p;
		nnode = self.getNodeNum();
		ndof = self.getDofNum();
		Fun = numpy.zeros(nnode);
		Fun[0] = l1 * (2 * l1 - 1);
		Fun[1] = l2 * (2 * l2 - 1);
		Fun[2] = l3 * (2 * l3 - 1);
		Fun[3] = l4 * (2 * l4 - 1);
		Fun[4] = 4 * l1 * l2;
		Fun[5] = 4 * l2 * l3;
		Fun[6] = 4 * l1 * l3;
		Fun[7] = 4 * l1 * l4;
		Fun[8] = 4 * l2 * l4;
		Fun[9] = 4 * l3 * l4;
		N = numpy.zeros((3, ndof));
		for i in range(0, nnode):
			for j in range(0, 3):
				N[j][i * 3 + j] = Fun[i];
		return N;

	def getShapeDerMatrix(self, p):
		[l1, l2, l3, l4, w] = p;
		nnode = self.getNodeNum();
		Der = numpy.zeros((3, nnode));
		for i in range(0, 3):
			 Der[i][i] = 4 * p[i] - 1;
			 Der[i][3] = 1 - 4 * p[3];

		Der[0][4] = 4 * l2;
		Der[1][4] = 4 * l1;

		Der[1][5] = 4 * l3;
		Der[2][5] = 4 * l2;

		Der[0][6] = 4 * l3;
		Der[2][6] = 4 * l1;

		Der[0][7] = 4 * (l4 - l1);
		Der[1][7] = -4 * l1;
		Der[2][7] = -4 * l1;

		Der[0][8] = -4 * l2;
		Der[1][8] = 4 * (l4 - l2);
		Der[2][8] = -4 * l2;

		Der[0][9] = -4 * l3;
		Der[1][9] = -4 * l3;
		Der[2][9] = 4 * (l4 - l3);
		return Der;
	
	def getStress(self):
		ndof = self.getDofNum();
		u = numpy.zeros((ndof, 1));
		B = numpy.zeros((6, ndof));
		k = 0;
		for n in self.nodes:
			dofs = n.getDofs();
			ndof = 0;
			for d in dofs:
				u[k + ndof][0] = n.values[d];
				ndof += 1;
			k += n.getDofNum();
		D = self.getStressStrainMatrix();
		points = [
			[1, 0, 0, 0, 0],
			[0, 1, 0, 0, 0],
			[0, 0, 1, 0, 0],
			[0, 0, 0, 1, 0]
		];
		stress = numpy.zeros((6, 4));
		k = 0;
		for p in points:
			B = self.getStrainMatrix(p);
			stress[:,k] = D.dot(B).dot(u).reshape((6,));
			k += 1;
		return stress;

	def addPload(self, load, v, n1, n2):
		numbers = [];
		nodes = [];
		k = -1;
		for i in range(0, 4):
			n = self.nodes[i];
			if n == n2:
				k = i;
			else:
				nodes.append(n);
		assert(k >= 0);
		if k == 0:
			numbers += [5, 8, 9];
		elif k == 1:
			numbers += [6, 7 ,9];
		elif k == 2:
			numbers += [4, 7, 8];
		elif k == 3:
			numbers += [4, 5, 6];

		A = TetElement.getTriangleArea(nodes);
		for no in numbers:
			n = self.nodes[no];
			dofs = n.getDofs();
			for dof in dofs:
				if dof == load.getDof():
					n.addLoad(load, v / len(numbers) * A);

class Element1D(Element, metaclass = abc.ABCMeta):
	def getLength(self):
		n1 = self.nodes[0];
		n2 = self.nodes[-1];
		return math.sqrt((n2.x - n1.x) ** 2 + (n2.y - n1.y) ** 2 + (n2.z - n1.z) ** 2);

class Truss(Element1D, metaclass = abc.ABCMeta):
	def getStressStrainMatrix(self):
		D = numpy.zeros((1, 1));
		D[0, 0] = self.property.material.getProperty(MaterialProperty.E);
		return D;
	
	def getStrainMatrix(self, p):
		L = self.getLength();
		Der = self.getShapeDerMatrix(p);
		ndof = self.getDofNum();
		nnode = self.getNodeNum();
		B = numpy.zeros((1, ndof));
		for i in range(0, nnode):
			B[0, i * 3] = Der[0, i] * 2 / L;
		return B;
	
	def getStiffMatrix(self):
		nnode = self.getNodeNum();
		ndof = self.getDofNum();
		Ke = numpy.zeros((ndof, ndof));
		D = self.getStressStrainMatrix();
		A = self.property.A;
		n1 = self.nodes[0];
		n2 = self.nodes[1];
		n = numpy.array([n2.x - n1.x, n2.y - n1.y, n2.z - n1.z]);
		L = numpy.linalg.norm(n);
		n = n / L;
		for p in self.points:
			B = self.getStrainMatrix(p);
			Ke += B.T.dot(D).dot(B) * L * A * p[-1];
		T = numpy.zeros(((ndof, ndof)));
		for i in range(0, nnode):
			T[i * 3, i * 3 : i * 3 + 3] = n;
		return T.T.dot(Ke).dot(T);

class Truss2(Truss):
	def __init__(self, eid, nodes, material):
		super(Truss2, self).__init__(eid, nodes, material, [[0, 1]]);

	def getShapeMatrix(self, p):
		[xi, w] = p;
		ndof = self.getDofNum();
		N = numpy.zeros((3, ndof));
		N[0, 0] = (1 - xi) / 2;
		N[0, 3] = (1 + xi) / 2;
		return N;

	def getShapeDerMatrix(self, p):
		nnode = self.getNodeNum();
		Der = numpy.zeros((1, nnode));
		Der[0, 0] = - 1 / 2;
		Der[0, 1] = 1 / 2;
		return Der;

if __name__ == '__main__':
	n1 = Node(1, 0, 0, 0);
	n2 = Node(2, 1, 1, 0);
	n3 = Node(3, 1, 0, 0);
	E = 1;
	rho = 7.9e3;
	nu = 0.3;
	A = 1e-4;

	m = Material();
	m.setProperty(MaterialProperty.E, E);
	m.setProperty(MaterialProperty.rho, rho);
	m.setProperty(MaterialProperty.nu, nu);
	
	p = Property1D(m, A);


	t = Truss2(1, [n1, n2], p);
	print(t.getStiffMatrix());
