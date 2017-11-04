from enum import Enum;

MaterialProperty = Enum('MaterialProperty', 'E nu rho k alpha T0');

class Material:
	def __init__(self):
		self.property = dict();

	def getProperty(self, p):
		return self.property[p];

	def setProperty(self, p, v):
		self.property[p] = v;
