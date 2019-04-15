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
	bdf.read_bdf('A4.bdf');
	model = Model();
	m = Material();
	m.setProperty(MaterialProperty.rho, rho);
	m.setProperty(MaterialProperty.E, E);
	m.setProperty(MaterialProperty.nu, nu);
	p = Property3D(m);
	nodes = dict();
	elements = dict();
	for i, n in bdf.nodes.items():
		node = Node(n.nid, n.xyz[0], n.xyz[1], n.xyz[2]);
		nodes[n.nid] = node;

	for i, e in bdf.elements.items():
		list = [];
		for nid in e.nodes:
			list.append(nodes[nid]);
		element = Tet10Element(e.eid, list, p);
		elements[e.eid] = element;
		model.addElement(element);


	for i, spc in bdf.spcs.items():
		for x in spc:
			if x.type == 'SPC1':
				for nid in x.nodes:
					nodes[nid].addConstraint(Constraint.X);
					nodes[nid].addConstraint(Constraint.Y);
					nodes[nid].addConstraint(Constraint.Z);

	for i, l in bdf.loads.items():
		for x in l:
			if x.type == 'FORCE':
				nodes[x.node].addLoad(Load.X, l[0].mag * x.xyz[0]);
				nodes[x.node].addLoad(Load.Y, l[0].mag * x.xyz[1]);
				nodes[x.node].addLoad(Load.Z, l[0].mag * x.xyz[2]);
			elif x.type == 'PLOAD4':
				p = x.pressures[0];
				v1 = p * x.nvector[0];
				v2 = p * x.nvector[1];
				v3 = p * x.nvector[2];

				e = elements[x.eids[0]];
				n1 = nodes[x.g1];
				n2 = nodes[x.g34];
				e.addPload(Load.X, v1, n1, n2);
				e.addPload(Load.Y, v2, n1, n2);
				e.addPload(Load.Z, v3, n1, n2);
	model.solve();
