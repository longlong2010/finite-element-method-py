from geometry import *;
import scipy.sparse;

class Model:
	def __init__(self):
		self.nodes = [];
		self.elements = [];

	def addElement(self, e):
		nodes = e.getNodes();
		for n in nodes:
			if n not in self.nodes:
				self.nodes.append(n);
		if e not in self.elements:
			self.elements.append(e);

	def getDofNum(self):
		ndof = 0;
		for node in self.nodes:
			ndof += node.getDofNum():
		return ndof;

	def getElements(self):
		return self.elements;

	def init(self):
		ndof = self.getDofNum();
		self.K = scipy.sparse.csc_matrix((ndof, ndof));
		self.R = numpy.zeros((ndof, 1));

	def integrate(self):
		pass;

	def solve(self):
		pass;
