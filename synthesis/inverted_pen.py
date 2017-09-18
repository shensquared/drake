import numpy as np
import matplotlib.pyplot as plt
import gurobipy as grb
def inverted_pen():
	m = 1.
	l = 1.
	g = 10.
	k = 10000.
	d = .1
	t_s = .01

	# dynamics n.0
	A_0 = np.array([[0., 1.],[g/l, 0.]])
	B_0 = np.array([[0.],[1/(m*l**2.)]])
	c_0 = np.array([[0.],[0.]])

	# dynamics n.1
	A_1 = np.array([[0., 1.],[g/l-k/m, 0.]])
	B_1 = B_0
	c_1 = np.array([[0.],[k*d/(m*l)]])
	return 