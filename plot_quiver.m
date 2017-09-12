function plot_quiver(x,xdot,sum_rho)
	figure(5)
	% combined=[-x(1)+x(2)-sum_rho;-x(1)-x(2)-sum_rho;x(1)-x(2)-sum_rho;x(1)+x(2)-sum_rho;];
	[a,b]=meshgrid(-sum_rho:sum_rho/50:sum_rho,-sum_rho:sum_rho/50:sum_rho);
	% regional=[min(dmsubs(combined,x,[a(:)';b(:)'])<=0)];
	z1=(dmsubs(xdot(1),x,[a(:)';b(:)']));
	z2=(dmsubs(xdot(2),x,[a(:)';b(:)']));
    z1=reshape(z1,size(a));
    z2=reshape(z2,size(a));
	quiver(a,b,z1,z2);
	title('$$\dot{x}$$','interpreter','latex','fontsize',20) 
	xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
	ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)


end