function [V,rho,all_V,sol_OK]=VanDerPol_PWA()
	do_plots=true;
	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);

	% xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
	% xdot=[-x(1);-x(2)];
	xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];

	rho=3e-4;
	delta_rho=3e-4;
	Vertices_values=zeros(4,1);
	sol_OK=true;
	level=0;

	[rho,Vertices_values,w,sol_OK]=cont_diamond_ring(x,xdot,rho,zeros(4,1),delta_rho,do_plots)
	[V,rho]=dis_diamond(x,xdot,do_plots,rho);
	[V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,level,rho,delta_rho,V,V,do_plots)

	% [rho,Vertices_values,w,sol_OK]=cont_diamond_ring(x,xdot,0,zeros(4,1),rho,do_plots);
	% [V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,varargin);
	% [V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,0,level,rho);
	% while(sol_OK)
	% 	disp(level);
	% 	[V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,rho,1,delta_rho,V,all_V);
	% 	level=level+1;
 %    end

end
