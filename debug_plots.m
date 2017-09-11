function debug_plots(x,f,sum_rho)
        figure(25)
        combined=[-x(1)+x(2)-sum_rho;-x(1)-x(2)-sum_rho;x(1)-x(2)-sum_rho;x(1)+x(2)-sum_rho;];
        [a,b]=meshgrid(-sum_rho:sum_rho/100:sum_rho,-sum_rho:sum_rho/100:sum_rho);
        % regional=[min(dmsubs(combined,x,[a(:)';b(:)'])<=0)];
        z=(dmsubs(f,x,[a(:)';b(:)']));
        z=reshape(z,size(a));
        surf(a,b,z);
        title('$$max(V)$$','interpreter','latex','fontsize',20) 
        xlabel('$$ x_1 $$','interpreter','latex','fontsize',15)
        ylabel('$$ x_2 $$','interpreter','latex','fontsize',15)
end