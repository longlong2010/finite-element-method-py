from enum import Enum;
import abc;

MaterialProperty = Enum('MaterialProperty', 'E nu rho k alpha T0');

class Material:
	def __init__(self):
		self.property = dict();

	def getProperty(self, p):
		return self.property[p];

	def setProperty(self, p, v):
		self.property[p] = v;

class Property(metaclass = abc.ABCMeta):
	def __init__(self, material):
		self.material = material;

class Property3D(Property):
	def __init__(self, material):
		super(Property3D, self).__init__(material);

class Property1D(Property):
	def __init__(self, material, A):
		super(Property1D, self).__init__(material);
		self.A = A;
