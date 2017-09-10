function [V,rho,all_V,sol_OK]=VanDerPol_PWA()
	do_plots=false;
	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);

	% xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
	% xdot=[-x(1);-x(2)];
	xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];

	rho=3e-4;
	delta_rho=3e-4;
	Vertices_values=zeros(4,1);
	w=zeros(2,4);
	sol_OK=true;
	level=0;
	[V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,0,level,rho);
	level=1;
	% while(sol_OK)
	% 	disp(level);
	% 	[V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,rho,1,delta_rho,V,all_V);
	% 	level=level+1;
 %    end

	% diamond_ring(x,xdot,rho,[7.6180;7.5818;7.6180;7.5818],.2);
% 	[rho,Vertices_values,w,sol_OK]=diamond_ring(x,xdot,rho,Vertices_values,1,do_plots);
% 	[rho,Vertices_values,w,sol_OK]=diamond_ring(x,xdot,rho,Vertices_values,.4,do_plots);

end
