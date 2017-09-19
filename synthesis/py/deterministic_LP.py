import numpy as np
import gurobipy as grb
import sympy as sp
from sympy import *
import mpmath
def deterministic_LP():
	rho=.5;
	resolution=4;
	
	row_verts=2*resolution+1;
	num_tris=(2*resolution)^2;
	num_verts=(row_verts)^2;
	
	x1,x2=sp.symbols('x1,x2', real=True)
	xdot=sp.Matrix([[x2], [-x1-x2*(x1*x1-1)]])
	df=xdot.jacobian([x1,x2])

	# [a,b]=np.meshgrid(-rho,rho/resolution,rho,-rho,rho/resolution,rho);
	# [centroid_a,centroid_b]=np.meshgrid((-3+1/resolution)*rho/3:rho/resolution:(3-2/resolution)*rho/3,(-3+1/resolution)*rho/3:rho/resolution:(3-2/resolution)*rho/3);
	# [centroid_c,centroid_d]=np.meshgrid((-3+2/resolution)*rho/3:rho/resolution:(3-1/resolution)*rho/3,(-3+2/resolution)*rho/3:rho/resolution:(3-1/resolution)*rho/3);
	# centroid_a=centroid_a.reshape([num_tris,1]);
	# centroid_b=centroid_b.reshape([num_tris,1]);
	# centroid_c=centroid_c.reshape([num_tris,1]);
	# centroid_d=centroid_d.reshape([num_tris,1]);


	return
deterministic_LP()