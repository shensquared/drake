function [V,rho,all_V,sol_OK]=VanDerPol_PWA()
	do_plots=false;
	verbose=false;
	checkDependency('spotless');
	checkDependency('mosek');
	x=msspoly('x',2);

	% xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
	% xdot=[-x(1);-x(2)];
	xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];

	last_rho=3;
	delta_rho=3e-2;
	last_V=zeros(4,1);
	all_V=[];
	sol_OK=true;
	level=1;
	[last_V,last_rho,all_V,sol_OK]=flat_diamond_rings(x,xdot,level,last_rho,delta_rho,last_V,all_V,do_plots,verbose)
	% [V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,level,last_rho,delta_rho,last_V,all_V,do_plots)
	% level=1;
	% while(sol_OK)
	% 	[last_V,last_rho,all_V,sol_OK]=flat_diamond_rings(x,xdot,level,last_rho,delta_rho,last_V,all_V,do_plots,verbose)
	% 	level=level+1;
	% end

	% [V,rho,all_V,sol_OK]=flat_diamond_rings(x,xdot,level,rho,delta_rho,V,all_V,do_plots)
	% [rho,Vertices_values,w,sol_OK]=cont_diamond_ring(x,xdot,rho,zeros(4,1),delta_rho,do_plots)
	% [V,rho]=dis_diamond(x,xdot,do_plots,rho);
	% [V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,level,rho,delta_rho,V,V,do_plots)

	% [rho,Vertices_values,w,sol_OK]=cont_diamond_ring(x,xdot,0,zeros(4,1),rho,do_plots);
	% [V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,varargin);
	% [V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,0,level,rho);
	% while(sol_OK)
	% 	disp(level);
	% 	[V,rho,all_V,sol_OK]=dis_diamond_ring(x,xdot,do_plots,rho,1,delta_rho,V,all_V);
	% 	level=level+1;
 %    end

end
