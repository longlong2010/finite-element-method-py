from load import *;
from geometry import *;
from property import *;

if __name__ == '__main__':
	E = 200e9;
	rho = 7.9e3;
	nu = 0.3;

	n1 = Node(1, 0, 0);
	n2 = Node(0, 1, 0);
	n3 = Node(1, 1, 0);
	n4 = Node(0, 0, 1);
	
	m = Material();
	m.setProperty(MaterialProperty.E, E);
	m.setProperty(MaterialProperty.rho, rho);
	m.setProperty(MaterialProperty.nu, nu);

	e = Tet4Element(n1, n2, n3, n4, m);
	print(e.getMassMatrix());
