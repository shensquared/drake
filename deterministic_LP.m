function deterministic_LP()
	sol_OK=false;

	rho=1.5;
	resolution=8;
	row_verts=2*resolution+1;
	num_tris=(2*resolution)^2;
	num_verts=(row_verts)^2;


	[a,b]=meshgrid(-rho:rho/resolution:rho,-rho:rho/resolution:rho);
	[centroid_a,centroid_b]=meshgrid((-3+1/resolution)*rho/3:rho/resolution:(3-2/resolution)*rho/3,(-3+1/resolution)*rho/3:rho/resolution:(3-2/resolution)*rho/3);
	[centroid_c,centroid_d]=meshgrid((-3+2/resolution)*rho/3:rho/resolution:(3-1/resolution)*rho/3,(-3+2/resolution)*rho/3:rho/resolution:(3-1/resolution)*rho/3);

	centroid_a=reshape(centroid_a,[num_tris,1]);
	centroid_b=reshape(centroid_b,[num_tris,1]);
	centroid_c=reshape(centroid_c,[num_tris,1]);
	centroid_d=reshape(centroid_d,[num_tris,1]);
	% clf
	scatter(reshape(a,[num_verts,1]),reshape(b,[num_verts,1])); hold on
	scatter(centroid_a,centroid_b); hold on
	scatter(centroid_c,centroid_d); hold on

	x=msspoly('x',2);
	% xdot = [-2*x(1); -2*x(2)];
	% xdot = [-2*x(1)+x(1)^3; -2*x(2)+x(2)^3];
	xdot = -[x(2); -x(1)-x(2).*(x(1).^2-1)];

	df=diff(xdot,x);
	A1=dmsubs(df(:,1),x,[centroid_a(:)';centroid_b(:)']);
	A2=dmsubs(df(:,2),x,[centroid_a(:)';centroid_b(:)']);
	A=[];
	for i =1:num_tris
		A=[A,A1(:,i),A2(:,i)];
	end

	% vert_values=find_V(x,xdot,A,rho,resolution,a,b);
	% scatter(centroid_c,centroid_d)
end

function vert_values=find_V(x,xdot,A,rho,resolution,a,b)
	row_verts=2*resolution+1;
	num_tris=(2*resolution)^2;
	num_verts=(row_verts)^2;
	prog=spotsosprog;
	prog=prog.withIndeterminate(x);
	[prog,vert_values]=prog.newPos(num_verts);
	[prog,slacks]=prog.newPos(num_tris);
	slacks=reshape(slacks,[2*resolution,2*resolution]);
	vert_values=reshape(vert_values,[row_verts,row_verts]);
	prog=prog.withEqs(vert_values(resolution+1,resolution+1));
	% the exterior values are the same
	prog=prog.withEqs(vert_values(1,1)-1);

	for i=1:row_verts
		for j=1:row_verts
			if i==1|j==1|i==row_verts|j==row_verts
				prog=prog.withEqs(vert_values(i,j)-1);
			end
		end
	end


	vert_values_no_top=vert_values(2:row_verts,:);
	vert_values_no_right=vert_values(:,1:row_verts-1);

	for i=1:2*resolution
		for j=1:2*resolution
			w1(i,j)=(resolution/(rho)).*(vert_values_no_top(i,j+1)-vert_values_no_top(i,j));
			w2(i,j)=(resolution/(rho)).*(vert_values_no_right(i+1,j)-vert_values_no_right(i,j));
		end
	end
	w1=reshape(w1,num_tris,1);
	w2=reshape(w2,num_tris,1);
	% for i=1:64
	% 	Vdot(i)=[w1(i),w2(i)]*A{i}*x
	% end
	xrep=repmat(x,[num_tris,1]);
	Vdot=[w1,w2]*A*xrep;
	Vdot=reshape(Vdot,[2*resolution,2*resolution]);

	for i=1:2*resolution
		for j=1:2*resolution
			% verts=[(a(i,j);b(i,j));;(a(i+1,j);b(i+1,j))];
			prog=prog.withPos(-slacks(i,j)-(subs(Vdot(i,j),x,[a(i,j);b(i,j)])));
			prog=prog.withPos(-slacks(i,j)-(subs(Vdot(i,j),x,[a(i+1,j);b(i+1,j)])));
			prog=prog.withPos(-slacks(i,j)-(subs(Vdot(i,j),x,[a(i,j);b(i,j)])));

		end
	end
	slacks=reshape(slacks,[num_tris,1]);

	% vdot<0 at the vertices, with f approxiamted at the centroid
	
	options = spot_sdp_default_options();
	options.verbose=1;
	sol=prog.minimize(-sum(slacks),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		values=double(sol.eval(vert_values))
		surf(a,b,values);
		title('piecewise $$\dot{V}$$','interpreter','latex','fontsize',20) 
		xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
		ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
	end

end



