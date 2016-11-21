classdef MyPendulum2States < PolynomialSystem
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
    function obj = MyPendulum2States(b)
      % Construct a new PendulumPlant
      obj = obj@PolynomialSystem(2,0,1,2,false,true,false);
      obj = setStateFrame(obj,CoordinateFrame('MyPendulumState',2,'x',{'theta','thetadot'}));
      if nargin>0 && ~isempty(b) % accept damping as optional input
        obj.b = b;
      end
      obj = setOutputFrame(obj,PendulumState);
    end

    function xdot = dynamicsRHS(~,~,x,u)
      p= MyPendulum2States;
      q=x(1);
      qd=x(2);
      qdd = (- p.m*p.g*p.lc*sin(q) - p.b*qd)/p.I;
      xdotopenloop=[qd;qdd];
      B=[0;1];
      K=readK(x);
      xdot=xdotopenloop+B*K;
    end

    function y=output(~,~,x,~)
      y=x(1:2);
    end
  end
end



function K=readK(state)
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
load('k.mat');
q=state(1);
qd=state(2);
newstate=[sin(q);cos(q)+1;qd];
K_gain=dmsubs(K',x,newstate)';
K=K_gain*newstate;
end










