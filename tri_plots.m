function tri_plots(x,xdot,samples,V,flags)
        mins=(min(samples,[],1));
        maxes=(max(samples,[],1));
        x1min=mins(1);x2min=mins(2);
        x1max=maxes(1);x2max=maxes(2);

        [a,b]=meshgrid(x1min:(x1max-x1min)/50:x1max,x2min:(x2max-x2min)/50:x2max);

        figure(1)
        clf
        V_values=(dmsubs(V,x,[a(:)';b(:)']));
        V_1_values=reshape(V_values(1,:),size(a));
        V_2_values=reshape(V_values(2,:),size(a));
        V_3_values=reshape(V_values(3,:),size(a));

        subplot(1,2,1)
        surf(a,b,V_1_values);hold on
        surf(a,b,V_2_values);hold on
        surf(a,b,V_3_values);hold on
        surf(a,b,ones(size(a)));
        title('piecewise affine $$V$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)

        subplot(1,2,2)
        w=diff(V,x);
        constraint1=[V(1)-1;V(2)-V(1);V(3)-V(1)];
        constraint2=[V(2)-1;V(1)-V(2);V(3)-V(2)];
        constraint3=[V(3)-1;V(1)-V(3);V(2)-V(3)];
        regional=[min(dmsubs(constraint1,x,[a(:)';b(:)'])<0);min(dmsubs(constraint2,x,[a(:)';b(:)'])<0);min(dmsubs(constraint3,x,[a(:)';b(:)'])<0);];
        z=(regional).*(dmsubs(w*xdot,x,[a(:)';b(:)']));
        z=sum(z,1);
        z=reshape(z,size(a));
        z(z==0) = NaN;
        surf(a,b,z);
        title('piecewise $$\dot{V}$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)

end