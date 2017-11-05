from load import *;
from geometry import *;
from property import *;
from model import *;

import csv;

if __name__ == '__main__':
	numpy.set_printoptions(threshold = numpy.nan);
	E = 200e9;
	rho = 7.9e3;
	nu = 0.3;

	model = Model();
	
	m = Material();
	m.setProperty(MaterialProperty.E, E);
	m.setProperty(MaterialProperty.rho, rho);
	m.setProperty(MaterialProperty.nu, nu);

	nodes = dict();
	l = 0;
	with open('demo.msh', 'r') as f:
		reader = csv.reader(f, delimiter = ' ');
		mode = '';
		for row in reader:
			if row[0] == '$Nodes':
				mode = 'Node';
				continue;
			elif row[0] == '$Elements':
				mode = 'Element';
				continue;

			if mode == 'Node':
				if len(row) == 4:
					n = Node(float(row[1]) / 1000, float(row[2]) / 1000, float(row[3]) / 1000);
					if abs(n.x + 5e-2) < 1e-10:
						n.addConstraint(Constraint.X);
						n.addConstraint(Constraint.Y);
						n.addConstraint(Constraint.Z);
					elif abs(n.x - 5e-2) < 1e-10:
						n.addLoad(Load.X, 10);
						l += 10;
					nodes[row[0]] = n;
			elif mode == 'Element':
				if len(row) == 9:
					n1 = nodes[row[5]];
					n2 = nodes[row[6]];
					n3 = nodes[row[7]];
					n4 = nodes[row[8]];
					model.addElement(Tet4Element(n1, n2, n3, n4, m));
	print(l);	
	model.solve();
	#n1 = Node(1, 0, 0);
	#n2 = Node(0, 1, 0);
	#n3 = Node(1, 1, 0);
	#n4 = Node(0, 0, 1);

	#e = Tet4Element(n1, n2, n3, n4, m);
	#n1.addConstraint(Constraint.X);
	#n1.addConstraint(Constraint.Y);
	#n1.addConstraint(Constraint.Z);
	#
	#n3.addConstraint(Constraint.X);
	#n3.addConstraint(Constraint.Y);
	#n3.addConstraint(Constraint.Z);
	#
	#n4.addConstraint(Constraint.X);
	#n4.addConstraint(Constraint.Y);
	#n4.addConstraint(Constraint.Z);

	#n2.addLoad(Load.X, 1000);

	#model.addElement(e);
	#model.solve();
