bound=2.5;
delta=.05;

x=-bound:delta:bound;
y=x;
[X,Y]=meshgrid(x,y);
V=(4.7092e-06)+(0.2837)*X.^2+(0.0087485)*X.^4+(0.18161)*Y.^2+(0.016551)*Y.^4+(-0.34678)*Y.*X+(0.065717)*Y.*(X.^3)+(0.10657)*(Y.^2).*(X.^2)+(-0.061218)*(Y.^3).*X;
% (0.00085577)+(0.28001)*X^2+(0.0098348)*X^4+(0.18945)*Y^2+(0.013964)*Y^4+(-0.35727)*Y*X+(0.068126)*Y*X^3+(0.10163)*Y^2*X^2+(-0.054716)*Y^3*X  ;
contourf(X,Y,100*V);
hold on
x1=X;
x2=Y;
% xdot = [x2; -x1-x2.*(x1.^2-1)];
quiver(X,Y,x2,-x1-x2.*(x1.^2-1));

% cleanVdot=(0.35727)*x1^2+(-0.068126)*x1^4+(0.021633)*x2^2+(0.0011409)*x2^4+(-0.17614)*x2*x1+(0.26147)*x2*x1^3+(-0.068126)*x2*x1^5+(0.19289)*x2^2*x1^2+(-0.20326)*x2^2*x1^4+(-0.016744)*x2^3*x1+(0.16415)*x2^3*x1^3+(-0.055857)*x2^4*x1^2;
