%%
% Prepare the road data
load IRI_737b
zdot = calc_zdot(road_x, road_z, 10); % Vehicle speed assumed to be 10

%%
% Simple demo of the quarter car model
%clear all
[system, objective, init, constraint] = qcar_spring_damper_open();
tf = 2;

% open loop
xd = zeros(7,1);
xd(1) = 0.010;          % d: spring wire diameter (m)
xd(2) = 0.129;          % D: spring helix diameter (m)
xd(3) = 0.106;          % p: spring pitch (m)
xd(4) = 3.57;           % Na: number of active spring coils (3 <= Na <= 15 for
xd(5) = 0.006;          % Do: orifice diameter (m)
xd(6) = 0.035;          % Dp: piston diameter (m)
xd(7) = 0.17;           % Ds: damper stroke (m)

opt = odeset('RelTol', 1e-2, 'AbsTol', 1e-3);
ramp_in = 25/100 * 10;
[t01,x01] = ode23(@(t,x)system.deriv(t,x,0, xd, @(t)ramp_in), [0 3], init.fun(xd),opt);
[t02,x02] = ode23(@(t,x)system.deriv(t,x,0, xd, @(t)lookup_u(zdot,t)), [0 2], init.fun(xd),opt);

%%
% Map the state back
%
xu01 = [x01'; zeros(1,length(t01))];
f0 = cell(1);
f0{1} = objective.fun([], xu01,t01);
xu02 = [x02'; zeros(1,length(t02))];
f0{2} = objective.fun([], xu02,t02);

%%
%{

s = DTSolverWithInputs();
s.system = system;
s.objective = objective;
s.initCondition = init;
s.designConstraint = constraint;


%%
% Optimization
s.n_control = 1;
s.t ={t01, t02};
s.input = {@(t)ramp_in, @(t)lookup_u(zdot, t)};
s.weight = {1e-2, 1};
s.pinit =xd; % initial parameterc
s.x0 = [];
s.lb = [0.005 0.05 0.05 3 0.005 0.03 0.1];
s.ub = [0.02 0.4 0.5 12 0.008 0.06 0.3];
fprintf('\nDirect Transcription\n');
starttime = cputime;
out1 = s.f_solve()
out1.cputime = cputime - starttime;
xu1 = out1.xu;

save('road_result');
f1 = cell(1);
f1{1} = objective.fun(out1.param, out1.xu{1},t01);
f1{2} = objective.fun(out1.param, out1.xu{2},t02);
%%


%%
% Generate initial control using lqr
[system2, objective2, init2] =  qcar_spring_damper_close();
% Get from original system
xd = out1.param;
J = system.jacobian(0,[],[],xd);
a = J(:,1:4);
b = J(:,5);
nx = 4;
QR = objective.QR();
Q = QR(1:nx, 1:nx);
R = QR(nx+1:end, nx+1:end);
k2 = lqr(a,b,Q,R);


[t22,x22] = ode23s(@(t,x)system2.deriv(t,x,[], [xd;-k2'],@(t)lookup_u(zdot,t)), t02, init2.fun(xd));
xu2 = cell(1);
xu2{2} = x22';

%%
[system2l, objective2l, init2l] = qcar_spring_damper_linear_close();
xd = out1.param;
J = system.jacobian(0,[],[],xd);
a = J(:,1:4);
b = J(:,5);
[f21, t21, xu21, k21] = lqr_control_design(xd, system2l, objective2l, init2l, b, t01);
xu2{1} = xu21;
f2 = cell(2);
f2{1} = f21;
f2{2} = objective2.fun([xd;-k2'], xu2{2},t22);


%}

%%
ss = SequentialSolverWithInputs();
ss.system = system;
ss.objective = objective;
ss.initCondition = init;
 
ss.n_control = 1;
ss.t = {t01,t02};
ss.input = {@(t)ramp_in, @(t)lookup_u(zdot, t)};
ss.weight = {1e-2, 1};
ss.xd0 =xd; % initial parameterc
 
ss.lb = [0.005 0.05 0.05 3 0.005 0.03 0.1];
ss.ub = [0.02 0.4 0.5 12 0.008 0.06 0.3];
fprintf('\nSequential Solver\n');
 
starttime = cputime;
out3 = ss.f_solve()
out3.cputime = cputime - starttime;
xu3 = out3.xu_frames;
 

f3{1} = objective.fun(out3.xd, out3.xu_frames{1},t01);
f3{2} = objective.fun(out3.xd, out3.xu_frames{2},t02);

%%

%plot_result(t02, x02, t02, out1.xu{2}, t02, xu2{2});
%plot_result(t01, x01, t01, out1.xu{1}, t01, xu2{1});

