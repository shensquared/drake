classdef d2model < PolynomialSystem
  methods
    function obj = d2model(b)
      % Construct a new PendulumPlant
      obj = obj@PolynomialSystem(4,0,1,2,false,true,false);
    end


    function xdot = dynamicsRHS(~,~,x,u)
      x
      h=x(3:4);
      A_x=[x(1), (x(1)); -x(1),0];
      A_h=[h(1), (h(1)); -h(1),0];
      [K,L]=d2read(x);
      xdot=[A_x*x(1:2);A_h*h]+[K;K;K;K];
      % [zeros(2,1);L]
      % xdot=[(x(2))*x(3);-x(1)*x(3);-p.g*x(1)/p.l-x(3)*p.b/(p.m*p.l.^2);...
      % (x(5))*x(6);-x(4)*x(6);-p.g*x(4)/p.l-x(5)*p.b/(p.m*p.l.^2);]+[0;0;K;0;0;K]+[zeros(3,1);L];
    end
    function y=output(~,~,x,~)
      y=x(1:2);
    end
  end  
end


function [K,L]=d2read(state)
checkDependency('spotless');
if checkDependency('mosek')
    solver = @spot_mosek;
else
    if ~checkDependency('sedumi')
        error('You will need to install Mosek or Sedumi to run this code');
    end
    solver = @spot_sedumi;
end

n=2;
x=msspoly('x',n);
h=msspoly('h',n);
X=[x;h];
load('d2KL.mat');
K=dmsubs(K',X,state);
% testing the controller only
% L=[0;0];
% L2=[0;0];
L1_old=dmsubs(L1,X,state);
L2_old=dmsubs(L2,X,state);
L=L1_old-L2_old;

end


