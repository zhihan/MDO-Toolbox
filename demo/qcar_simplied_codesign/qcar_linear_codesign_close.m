function [system, objective, init] = qcar_linear_codesign_close()
% Simple linear quarter car model

param.ms = 325;               % 1/4 sprung mass (kg)
param.mus = 65;               % 1/4 unsprung mass (kg)
param.kus = 232.5e3;          % tire stiffness (N/m)

param.ct = 0 ;                 % I cannot find ct in James' original code
 

rgrade = 25;            % ramp grade for rattle space test (percent)
v1 = 10;                % vehicle velocity for rattle space test (m/s)
param.v = v1;
param.rgrade = rgrade;
k0 = [0;0;0;0];
xd = [80e+3; 8e+3];
param.xs = steady_state([xd;k0], param);


system = System;
system.parameters =param;
system.deriv  = @(t,x,u,xd) f(t,x,u,xd, param);
system.jacobian = @(t,x,u,xd) jacobian(t,xd, param);
system.jacobian_d = @(t,x,u,xd) jacobian_d(t,x, param);

objective = Objective;
objective.fun = @(x_d, xin, t) rattle_space(x_d, xin, t, param);
objective.QR = @QR;
init = InitCondition;
init.n_variables = 0;
init.fun = @(x_d)X0(x_d, param);
init.jacobian = @X0Jacobian;

end

function xs = steady_state(xd, param)
J = jacobian(0,xd, param);
A = J(1:4, 1:4);
b_ramp = [-1; param.ct/param.mus; 0; 0];
xs = - (A \ b_ramp) * param.rgrade/100 * param.v;
end

function dx = f(t, x, u, xd, param)
A = jacobian(t, xd, param);
dx = A * x;
end


function J = jacobian(t, xd, param)
ks = xd(1);
cs = xd(2);
K = xd(3:end)';
A = [ 0 1 0 0;
    -param.kus/param.mus -cs/param.mus ks/param.mus cs/param.mus;
    0 -1 0 1;
    0 cs/param.ms -ks/param.ms -cs/param.ms];
B = [0 -1/param.mus 0 1/param.ms]'; 
J = A + B*K;
end


function J = jacobian_d(t, x, param)
J1 = [ 0 0; x(3)/param.mus (-x(2) + x(4))/param.mus; ...
       0 0; -x(3)/param.ms (x(2) - x(4))/param.ms];
B = [0 -1/param.mus 0 1/param.ms]'; 
J2 = B * x';

J = [J1,J2];
end

function val = rattle_space(x_d, xu, t, param)
xu_frame = xu;
val = 0;
R = Rattle(x_d, param);
for i=1:(size(xu_frame,2)-1)
    val = val+ (t(i+1)-t(i))*xu_frame(:,i)'*R*xu_frame(:,i);
end
end

function J = rattle_space_gradient(x_d, xu_frame, t, param)
nt = size(xu_frame,2);
R = Rattle(x_d, param);
J = zeros(numel(xu_frame),1);
n = size(xu_frame,1);
for i=1:(nt-1)
    J((i-1)*n+1:i*n) = 2*(t(i+1)-t(i))*R*xu_frame(:,i);
end
end

function x0 = X0(xd, param)
x0 = [0;0;0;0] - steady_state(xd, param);
end

function J = X0Jacobian(xd, param)
J = perturb_x(@(xd)X0(xd,param), xd);
end

function R = QR()
r1 = 3e+4;
r3 = 3e+4;
q = 1e-6;
R = diag([r1, 0, r3, 0, q]);
end

function R = Rattle(xd, param)
c1 = [1 0 0 0];
r1 = 3e+4;
c3 = [0 0 1 0];
r3 = 3e+4;
k = xd(3:end);
q =  1e-6;
R = r1 * (c1' * c1) + r3* (c3'*c3) + q*(k*k');
end
