from geometry import *;
import scipy.sparse;
from scipy.sparse import linalg;
import time;

class Model:
	def __init__(self):
		self.nodes = [];
		self.elements = [];

	def addElement(self, e):
		nodes = e.getNodes();
		for n in nodes:
			if n not in self.nodes:
				if len(self.nodes) > 0:
					_n = self.nodes[-1];
					n.setFid( _n.getFid() + _n.getDofNum());
				self.nodes.append(n);
		if e not in self.elements:
			self.elements.append(e);

	def getDofNum(self):
		ndof = 0;
		for node in self.nodes:
			ndof += node.getDofNum();
		return ndof;

	def getElements(self):
		return self.elements;

	def init(self):
		ndof = self.getDofNum();
		self.K = scipy.sparse.dok_matrix((ndof, ndof));
		self.R = numpy.zeros((ndof, 1));

	def integrate(self):
		ndof = self.getDofNum();
		for e in self.elements:
			nodes = e.getNodes();
			size = e.getDofNum();
		
			m = dict();
			l = 0;
			for n in nodes:
				k = n.getFid();
				for i in range(0, n.getDofNum()):
					m[l] = k + i;
					l += 1;

			Ke = e.getStiffMatrix();
			for i in range(0, size):
				for j in range(0, size):
					mi = m[i];
					mj = m[j];
					self.K[mi, mj] += Ke[i][j];

	def integrateLoad(self):
		for n in self.nodes:
			dofs = n.getDofs();
			dofn = 0;
			k = n.getFid();
			for d in dofs:
				loads = n.getLoads();
				for l, v in loads.items():
					if l.getDof() == d:
						self.R[k + dofn][0] = v;
				dofn += 1;

	def addConstraint(self):
		ndof = self.getDofNum();
		for n in self.nodes:
			constraints = n.getConstraints();
			dofs = n.getDofs();
			dofn = 0;
			k = n.getFid();
			for d in dofs:
				for c, v in constraints.items():
					if c.getDof() == d:
						self.K[k + dofn, k + dofn] += 1e20;
						self.R[k + dofn][0] = v * self.K[k + dofn, k + dofn];
				dofn += 1;
	
	def solveEquations(self):
		u = linalg.cg(self.K.tocsr(), self.R)[0];
		k = 0;
		for n in self.nodes:
			dofs = n.getDofs();
			dofn = 0;
			for d in dofs:
				n.setValue(d, u[k + dofn]);
				dofn += 1;
			k += n.getDofNum();

	def outputResult(self):
		print('View "Disp" {');
		for e in self.elements:
			stress = e.getStress();
			print('SS (', end='');
			l1 = [];
			l2 = [];
			for i in range(0, 4):
				l1.append(e.nodes[i].x);
				l1.append(e.nodes[i].y);
				l1.append(e.nodes[i].z);
				#l2.append(stress[2][i]);
				l2.append(e.nodes[i].values[Dof.X]);
			print(','.join(map(str, l1)), end='');
			print(') {', end='');
			print(','.join(map(str, l2)), end='');
			print('};');
		print('};');

	def solve(self):
		self.init();
		self.integrate();
		self.integrateLoad();
		self.addConstraint();
		self.solveEquations();
		self.outputResult();	
