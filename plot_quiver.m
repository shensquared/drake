function plot_quiver(x,xdot,sum_rho)
	figure(5)
	% combined=[-x(1)+x(2)-sum_rho;-x(1)-x(2)-sum_rho;x(1)-x(2)-sum_rho;x(1)+x(2)-sum_rho;];
	[a,b]=meshgrid(-sum_rho:sum_rho/100:sum_rho,-sum_rho:sum_rho/100:sum_rho);
	% regional=[min(dmsubs(combined,x,[a(:)';b(:)'])<=0)];
	z=(dmsubs(xdot,x,[a(:)';b(:)']));
	quiver(a,b,z(1),z(2));
	title('$$\dot{x}$$','interpreter','latex','fontsize',20) 
	xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
	ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)


end