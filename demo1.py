from geometry import *;
from model import *;
from property import *;
from pyNastran.bdf.bdf import BDF, read_bdf;
import numpy;
import time;
import sys;

if __name__ == '__main__':
	numpy.set_printoptions(threshold = sys.maxsize);
	E = 200e3;
	rho = 7.9e-9;
	nu = 0.3;
	bdf = BDF(debug=False);
	bdf.read_bdf('A1.bdf');
	model = Model();
	m = Material();
	m.setProperty(MaterialProperty.rho, rho);
	m.setProperty(MaterialProperty.E, E);
	m.setProperty(MaterialProperty.nu, nu);
	p = Property3D(m);
	nodes = dict();
	for i, n in bdf.nodes.items():
		node = Node(n.nid, n.xyz[0], n.xyz[1], n.xyz[2]);
		nodes[n.nid] = node;

	for i, e in bdf.elements.items():
		list = [];
		for n in e.nodes:
			list.append(nodes[n]); 
		element = Tet10Element(e.eid, list, p);
		model.addElement(element);


	for i, spc in bdf.spcs.items():
		for x in spc:
			if x.type == 'SPC1':
				for n in x.nodes:
					nodes[n].addConstraint(Constraint.X);
					nodes[n].addConstraint(Constraint.Y);
					nodes[n].addConstraint(Constraint.Z);

	for i, l in bdf.loads.items():
		for x in l:
			if x.type == 'FORCE':
				nodes[x.node].addLoad(Load.X, l[0].mag * x.xyz[0]);
				nodes[x.node].addLoad(Load.Y, l[0].mag * x.xyz[1]);
				nodes[x.node].addLoad(Load.Z, l[0].mag * x.xyz[2]);
	t1 = time.time();
	model.solve();
	t2 = time.time();
	print(t2 - t1);
