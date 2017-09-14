function deterministic_LP()
	sol_OK=false;

	rho=1;
	resolution=4;
	num_tris=(2*resolution)^2;
	[a,b]=meshgrid(-rho:rho/resolution:rho,-rho:rho/resolution:rho);
	[centroid_a,centroid_b]=meshgrid(-11/12*rho:rho/4:5/6*rho,-11/12*rho:rho/4:5/6*rho);
	% [centroid_c,centroid_d]=meshgrid(-5/6*rho:rho/4:11/12*rho,-5/6*rho:rho/4:11/12*rho);

	% a=reshape(a,[81,1]);
	% b=reshape(b,[81,1]);
	centroid_a=reshape(centroid_a,[num_tris,1]);
	centroid_b=reshape(centroid_b,[num_tris,1]);
	% centroid_c=reshape(centroid_c,[num_tris,1]);
	% centroid_d=reshape(centroid_d,[num_tris,1]);
	% clf
	% scatter(a,b); hold on
	% scatter(centroid_a,centroid_b)
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
	% A=mat2cell(A,2,2*ones(64,1));

	vert_values=find_V(x,xdot,A,rho,resolution,a,b);
	% scatter(centroid_c,centroid_d)
	% A=dmsubs(diff(xdot,x),x,centroids);
end

function vert_values=find_V(x,xdot,A,rho,resolution,a,b)
	row_verts=2*resolution+1;
	num_tris=(2*resolution)^2
	num_verts=(row_verts)^2;
	prog=spotsosprog;
	prog=prog.withIndeterminate(x);
	[prog,vert_values]=prog.newPos(num_verts);
	[prog,slacks]=prog.newPos(num_tris);
	slacks=reshape(slacks,[2*resolution,2*resolution])
	vert_values=reshape(vert_values,[row_verts,row_verts]);
	prog=prog.withEqs(vert_values(resolution+1,resolution+1));
	% the exterior values are the same
	for i=1:row_verts
		for j=1:row_verts
			prog=prog.withEqs(vert_values(1,1)-vert_values(1,2));
			prog=prog.withEqs(vert_values(1,1)-vert_values(1,3));
			prog=prog.withEqs(vert_values(1,1)-vert_values(1,4));
			prog=prog.withEqs(vert_values(1,1)-vert_values(1,5));
			prog=prog.withEqs(vert_values(1,1)-vert_values(1,6));
			prog=prog.withEqs(vert_values(1,1)-vert_values(1,7));
			prog=prog.withEqs(vert_values(1,1)-vert_values(1,8));
			prog=prog.withEqs(vert_values(1,1)-vert_values(1,9));
			prog=prog.withEqs(vert_values(1,1)-vert_values(2,1));
			prog=prog.withEqs(vert_values(1,1)-vert_values(3,1));
			prog=prog.withEqs(vert_values(1,1)-vert_values(4,1));
			prog=prog.withEqs(vert_values(1,1)-vert_values(5,1));
			prog=prog.withEqs(vert_values(1,1)-vert_values(6,1));
			prog=prog.withEqs(vert_values(1,1)-vert_values(7,1));
			prog=prog.withEqs(vert_values(1,1)-vert_values(8,1));
			prog=prog.withEqs(vert_values(1,1)-vert_values(9,1));
			prog=prog.withEqs(vert_values(1,1)-vert_values(9,2));
			prog=prog.withEqs(vert_values(1,1)-vert_values(9,3));
			prog=prog.withEqs(vert_values(1,1)-vert_values(9,4));
			prog=prog.withEqs(vert_values(1,1)-vert_values(9,5));
			prog=prog.withEqs(vert_values(1,1)-vert_values(9,6));
			prog=prog.withEqs(vert_values(1,1)-vert_values(9,7));
			prog=prog.withEqs(vert_values(1,1)-vert_values(9,8));
			prog=prog.withEqs(vert_values(1,1)-vert_values(9,9));
			prog=prog.withEqs(vert_values(1,1)-vert_values(8,9));
			prog=prog.withEqs(vert_values(1,1)-vert_values(7,9));
			prog=prog.withEqs(vert_values(1,1)-vert_values(6,9));
			prog=prog.withEqs(vert_values(1,1)-vert_values(5,9));
			prog=prog.withEqs(vert_values(1,1)-vert_values(4,9));
			prog=prog.withEqs(vert_values(1,1)-vert_values(3,9));
			prog=prog.withEqs(vert_values(1,1)-vert_values(2,9));
			prog=prog.withEqs(vert_values(1,1)-1);
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
	sol=prog.minimize(-sum(0),@spot_mosek,options);
	if sol.status==spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
		sol_OK=true;
		values=double(sol.eval(vert_values))
		surf(a,b,values);
	end
end



