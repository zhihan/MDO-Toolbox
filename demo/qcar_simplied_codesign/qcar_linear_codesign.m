function [system, objective, init] = qcar_linear_codesign()
% Simple linear quarter car model

% The linear model uses a transformation of the state space 
%   x_l' = x - x_s where x_s = A \ d is the steady state solution.
% This is needed because the state feedback is then formulated as 
%   u = k x_l
% and then applied to the original system.

param.ms = 325;               % 1/4 sprung mass (kg)
param.mus = 65;               % 1/4 unsprung mass (kg)
param.kus = 232.5e3;          % tire stiffness (N/m)

param.ct = 0 ;          % Assume
 
rgrade = 25;            % ramp grade for rattle space test (percent)
v1 = 10;                % vehicle velocity for rattle space test (m/s)
param.v = v1;
param.rgrade = rgrade;
xd = [80e+3; 8e+3];
param.xs = steady_state(xd, param);

system = System;
system.parameters =param;
system.deriv  = @(t,x,u,xd) f(t,x,u,xd, param);
system.jacobian = @(t,x,u,xd) jacobian(t,xd, param);
system.jacobian_d = @(t,x,u,xd) jacobian_d(t,x, xd, param);

objective = Objective;
objective.fun = @(x_d,xin, t) rattle_space(xin, t, param);
objective.gradient = @(x_d,xin,t) rattle_space_gradient(xin, t, param);
objective.QR = @Rattle;

init = InitCondition;
init.n_variables = 0;
init.fun = @(xd)X0(xd, param);
init.jacobian = @(xd) X0Jacobian(xd, param);

end

function dx = f(t, x, u, xd, param)
J = jacobian(t, xd, param);
A = J(1:4,1:4);
B = J(1:4, 5);

dx = A* x + B * u ;
end

function xs = steady_state(xd, param)
% Calculate the steady state equation
J = jacobian(0,xd, param);
A = J(1:4, 1:4);
b_ramp = [-1; param.ct/param.mus; 0; 0];
xs = - (A \ b_ramp) * param.rgrade/100 * param.v;
end

function J = jacobian(t, xd, param)
ks = xd(1);
cs = xd(2);

A = [ 0 1 0 0;
    -param.kus/param.mus -cs/param.mus ks/param.mus cs/param.mus;
    0 -1 0 1;
    0 cs/param.ms -ks/param.ms -cs/param.ms];
B = [0 -1/param.mus 0 1/param.ms]'; 
J = [A, B];
end

function J = jacobian_d(t, x, xd, param)
J = [ 0 0; x(3)/param.mus (-x(2) + x(4))/param.mus; ...
       0 0; -x(3)/param.ms (x(2) - x(4))/param.ms];
end

function val = rattle_space(xu, t, param)
xu_frame = xu;
val = 0;
R = Rattle();
for i=1:(size(xu_frame,2)-1)
    val = val+ (t(i+1)-t(i))*xu_frame(:,i)'*R*xu_frame(:,i);
end
end

function J = rattle_space_gradient(xu_frame, t, param)
nt = numel(xu_frame)/5;
R = Rattle();
J = zeros(numel(xu_frame),1);
for i=1:(nt-1)
    J((i-1)*5+1:i*5) = 2*(t(i+1)-t(i))*R*xu_frame(:,i);
end
end

function x0 = X0(xd, param)
x0 = [0;0;0;0] - steady_state(xd, param);
end

function J = X0Jacobian(xd, param)
J = perturb_x(@(xd)X0(xd, param), xd);
end


function R = Rattle()
r1 = 3e+4;
r3 = 3e+4;
q = 1e-6;
R = diag([r1, 0, r3, 0, q]);
end