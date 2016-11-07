classdef MyPendulum < PolynomialSystem
% Defines the dynamics for the Pendulum.
  
  properties
    m = 1;   % kg
    l = .5;  % m
    b = 0.1; % kg m^2 /s
    lc = .5; % m
    I = .25; %m*l^2; % kg*m^2
    g = 9.81; % m/s^2
    
    xG;
    uG;
  end
  
  methods
    function obj = MyPendulum(b)
      % Construct a new PendulumPlant
      obj = obj@PolynomialSystem(6,0,1,3,false,true,false);
      if nargin>0 && ~isempty(b) % accept damping as optional input
        obj.b = b;
      end
      % obj = setStateFrame(obj,PendulumState);
      % obj = setOutputFrame(obj,PendulumState);
      % obj = setStateFrame(obj,CoordinateFrame('MyPendulumState',4,'x',{'theta','thetadot','thetahat','thetahatdot'}));
      obj = setStateFrame(obj,CoordinateFrame('MyPendulumState',6,'x',{'sin','cos','thetadot','sinhat','coshat','thetahatdot'}));
      % obj = setOutputFrame(obj,MyPendulumState);

    end


    function xdot = dynamicsRHS(~,~,x,u)
      % x = [s;c;thetadot];
      x
      p = MyPendulum;
      h=x(4:6);
      A_x=[x(1), 0, (x(2)-1);x(1), 0, -x(1); -p.g/p.l, 0, -p.b/(p.m*p.l.^2)]
      A_h=[h(1), 0, (h(2)-1);h(1), 0, -h(1); -p.g/p.l, 0, -p.b/(p.m*p.l.^2)]
      [K,L]=DOFreadPartition(x);
      xdot=[A_x*x(1:3);A_h*h]+[K;K;K;K;K;K]+[zeros(3,1);L]
      % xdot=[(x(2))*x(3);-x(1)*x(3);-p.g*x(1)/p.l-x(3)*p.b/(p.m*p.l.^2);...
      % (x(5))*x(6);-x(4)*x(6);-p.g*x(4)/p.l-x(5)*p.b/(p.m*p.l.^2);]+[0;0;K;0;0;K]+[zeros(3,1);L];
    end



    function y=output(~,~,x,~)
      y=x(1:3);
    end
  end  
end


function [K,L]=DOFreadPartition(state)

  checkDependency('spotless');
if checkDependency('mosek')
    solver = @spot_mosek;
else
    if ~checkDependency('sedumi')
        error('You will need to install Mosek or Sedumi to run this code');
    end
    solver = @spot_sedumi;
end

n=3;
x=msspoly('x',n);
h=msspoly('h',n);
X=[x;h];
% observed=[x(1);x(2)];
load('partitionDOF_nonuni.mat')

% s1=state(1);
% s2=state(2);
% s3=state(3);
% sh1=state(4);
% sh2=state(5);
% sh3=state(6);

K=dmsubs(K',X,state)
% testing the controller only
% L=[0;0];
% L2=[0;0];

L1_old=dmsubs(L1,X,state);
L2_old=dmsubs(L2,X,state);
% L=L1_old-L2_old;
L=[0;0;0];
end










