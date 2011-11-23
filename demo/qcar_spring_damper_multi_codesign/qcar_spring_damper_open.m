function [system, objective, init, constraints] = qcar_spring_damper_open()
% Complex linear quarter car model

% The linear model uses a transformation of the state space 
%   x_l' = x - x_s where x_s = A \ d is the steady state solution.
% This is needed because the state feedback is then formulated as 
%   u = k x_l
% and then applied to the original system.

param.ms = 325;               % 1/4 sprung mass (kg)
param.mus = 65;               % 1/4 unsprung mass (kg)
param.kus = 232.5e3;          % tire stiffness (N/m)

param.ct = 0 ;          % I cannot find ct in James' original code
 
rgrade = 25;            % ramp grade for rattle space test (percent)
v1 = 10;                % vehicle velocity for rattle space test (m/s)
param.v = v1;
param.rgrade = rgrade;
% Design variables
d = 0.01;   % d: spring wire diameter 
D = 0.129;  % D: spring helix diameter
p = 0.106;  % p: spring pitch
Na = 3.57;   % Na: number of active spring coils 

Do = 0.006;          % Do: orifice diameter (m)
Dp = 0.035;          % Dp: piston diameter (m)
Ds = 0.17;           % Ds: damper stroke (m)

spring_vars = [d;D;p;Na];
damper_vars = [Do; Dp; Ds];
xd = [spring_vars; damper_vars];

system = System;
system.parameters =param;
system.deriv  = @(t,x,u,xd, input) f(t,x,u,xd, param, input);
system.jacobian = @(t,x,u,xd) jacobian(t,xd, param);
system.jacobian_d = @(t,x,u,xd) jacobian_d(t,x, xd, param);

objective = Objective;
objective.fun = @(x_d,xin, t) rattle_space(xin, t, param);
objective.gradient = @(x_d,xin,t) rattle_space_gradient(xin, t, param);
objective.QR = @Rattle;
objective.hessian =  @(x_d,xin,t) rattle_space_hessian(xin, t, param);

init = InitCondition;
init.n_variables = 0;
init.fun = @(xd)X0(xd, param);
init.jacobian = @(xd) X0Jacobian(xd, param);

constraints = Constraint;
constraints.fun = @(xd) spring_damper_constraints(xd, param);
constraints.jacobian = @(xd)spring_damper_constraints_jacobian(xd,param);
constraints.hessian = @(xd, lambda)constraint_hessian(xd, param, lambda);
end

function dx = f(t, x, u, xd, param, input)

J = jacobian(t, xd, param);
A = J(1:4,1:4);
B = J(1:4, 5);
b_ramp = [-1; param.ct/param.mus; 0; 0];
if ~isempty(input)
    dx= A* x + B * u + b_ramp * input(t);
else
    dx = A*x + B*u ;
end
end


function J = jacobian(t, xd, param)
ks = spring_constant(xd);
cs = damper_constant(xd);

A = [ 0 1 0 0;
    -param.kus/param.mus -cs/param.mus ks/param.mus cs/param.mus;
    0 -1 0 1;
    0 cs/param.ms -ks/param.ms -cs/param.ms];
B = [0 -1/param.mus 0 1/param.ms]'; 
J = [A, B];
end

function ks = spring_constant(xd)
% Calculate the spring constant ks from the design variables
d = xd(1);
D = xd(2);
p = xd(3);
Na = xd(4);
G = 77.2e+9; 
C = D/d;
ks = d^4 * G/ (8* D^3*Na*(1+1/2/C^2));
end

function cs = damper_constant(xd)
% Calculate the damper constant cs from the design variables
Do = xd(5);     % orifice diameter
Dp = xd(6);     % piston diameter
Ds = xd(7);     % damper stroke 

rho1 = 850;             % density of oil (kg/m^3)
kv = 7500;              % damper valve spring rate (N/m)
Cd = 0.7;               % valve discharge coefficient
Afa = 0.1;              % damper valve area factor 
Pallow = 4.75e6;        % maximum allowed damper pressure (Pa)
xi = 0.9;               % damper valve circumference factor

Ao = pi*Do^2/4;         % valve pressure surface area (M^2)
xm = Ao*Pallow/kv;      % lift @ Pallow (m)
C2 = xi*Afa*sqrt(xm);   % damper valve parameter (m^0.5)

cs = Dp^4*sqrt(kv* pi*rho1/2)/(8*Cd*Do^2*C2); %

end

function J = jacobian_d(t, x, xd, param)
J = perturb_x(@(xd)f(t,x,0, xd, param, []), xd);
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

function H = rattle_space_hessian(xin, t, param)
nt =  numel(xin)/5;
R = Rattle();
H = sparse([],[], numel(xin),numel(xin), nt*5);
for i=1:nt
    for j=1:5
        H((i-1)*5+j, (i-1)*5+j) = R(j,j);
    end
end
end

function H = constraint_hessian(xd, param, lambda)
% Hessian for algebraic constraints
F = hessian_x(@(xd) spring_damper_constraints(xd, param), xd);
H = zeros(size(F{1}));
for k=1:length(F)
    H = H + F{k} * lambda(k);
end
end

function x0 = X0(xd, param)
x0 = [0;0;0;0];
end

function J = X0Jacobian(xd, param)
x0 = X0(xd, param);
J = zeros(size(x0,1), size(xd,1));
end


function R = Rattle()
r1 = 3e+4;
r3 = 3e+4;
q = 1e-6;
R = diag([r1, 0, r3, 0, q]);
end

function g = spring_damper_constraints(xd, param)
d = xd(1);
D = xd(2);
p = xd(3);
Na = xd(4);
Do = xd(5);     % orifice diameter
Dp = xd(6);     % piston diameter
Ds = xd(7);     % damper stroke 

G = 77.2e+9; 
L0max = 0.40;
DoMax = 0.25;   
nd = 1.2;               % design factor (safety factor)
A = 1974;
m = 0.108;
Sut = A*1e6/(d^m);
Ssy = 0.65*Sut; 
dsc = 0.009;            % damper-spring clearance (m)
dwt = 0.002;            % damper wall thickness (m)
LB = 0.02;              % bump stop length (m)
ld1 = 0.02;             % axial length required for upper damper mounting/valving
ld2 = 0.04;             % axial length required for lower damper mounting/valving
ld3 = 0.02;             % axial length required for lower damper casing extension

C = D/d;
L0 = p*Na + 2*d;
Ls = d*(Na + 1.75 - 1);
ks = d^4 * G/ (8* D^3*Na*(1+1/2/C^2));
Fs = ks*(L0 - Ls);
KB = (4*C+2)/(4*C-3);  % Bergstrausser factor
tau_s = KB*8*Fs*D/(pi*d^3); 
delta_g = param.ms*9.81/ks; 

g = zeros(9,1);
g(1) = 4 - C;           % spring manufacturing constraint
g(2) = C - 12;          % spring manufacturing constraint
g(3) = L0 - 5.26*D;     % spring stability/buckling constraint for squared ground ends
g(4) = L0 - L0max;      % spring packaging constraint
g(5) = d + D - DoMax;   % spring packaging constraint (shouldn't be active: constrained by g3)
g(6) = d - D + Dp + 2*(dsc+dwt);        % spring-damper interference constraint
% g(7) = deltamax - L0 + Ls + LB + delta_g;
% suspension stop constraint (ensures don't hit stops during test bump).
g(7) = (tau_s*nd - Ssy)/Ssy;            % shear stress constraint (note scaling)
% g(9) = 0.15+1-(L0-Ls)/(delta_g+deltamaxf*1.1); % MULT OF STROKE: ensures lineari

g(8) = L0 - Ls - Ds;                   % ensures adequate damper stroke
g(9) = 2*Ds + ld1 + ld2 - L0max;       % ensures enough space for damper


end

function J = spring_damper_constraints_jacobian(xd, param)
J = perturb_x(@(xd) spring_damper_constraints(xd, param), xd, param);
end

