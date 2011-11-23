% Simple demo of the quarter car model
%clear all
[system, objective, init] = qcar_linear_codesign();
tf = 2;

% open loop
xd = [80e+3; 4e+3];
[t0,x0] = ode23s(@(t,x)system.deriv(t,x,0, xd), [0 tf], init.fun(xd));
%%
% Map the state back
%
xu0 = [x0'; zeros(1,length(t0))];
f0 = objective.fun([], xu0,t0);
xu0 = xu0 + repmat([system.parameters.xs; 0], 1, size(xu0,2));


%%
s = DTSolver();
s.system = system;
s.objective = objective;
s.initCondition = init;


s.n_control = 1;
s.t =t0;
s.xd0 =xd; % initial parameterc
s.lb = [40e+3 2e+3];
s.ub = [120e+3 8e+3];
fprintf('\nDirect Transcription\n');
tic,
out1 = s.f_solve()
toc

%%
% Map the state back to the original system
%
% 
t1 = t0;
xd = out1.xd;
f1 = objective.fun(out1.xd, out1.xu,t1);
xu1 = out1.xu + repmat([init.fun(out1.xd); 0], 1, size(out1.xu,2));


%%
% Nested optimization
J = system.jacobian(0,[],[],xd);
a = J(:,1:4);
b = J(:,5);
options = optimset('display','iter','algorithm','interior-point',...
                       'MaxIter', 5000, 'MaxFunEvals', 100, ...
                       'UseParallel','never', 'TolX', 1e-6,...
                       'GradConstr','off', 'GradObj','off');
lb = [40e+3 2e+3];
ub = [120e+3 8e+3];
xd = [80e+3; 4e+3];
[system2, objective2, init2] = qcar_linear_codesign_close();



[xopt,fopt,flag,output]=fmincon(@(xd)lqr_control_design(xd, system2, objective2, init2, b, tf), ...
    xd,[],[],[],[], ...
      lb,ub,[],options);
 out2.xopt = xopt;
 out2.fopt = fopt;
 out2.flag = flag;
 out2.output = output;
 
%%
%
% Compute the result
% lqr result
[f2, t2, xu2] = lqr_control_design(out2.xopt, system2, objective2, init2, b, tf);


%%
%
ss = SequentialSolver();
ss.system = system;
ss.objective = objective;
ss.initCondition = init;

ss.n_control = 1;
ss.t = t0;
ss.xd0 =[80e+3; 4e+3]; % initial parameterc
ss.lb = [40e+3 2e+3];
ss.ub = [120e+3 8e+3];

fprintf('\nSequential Solver\n');
tic,
out3 = ss.f_solve();
toc
%
t3 = t0;
f3 = objective.fun(xd, out3.xu_frames,t3);
xu3 = out3.xu_frames + repmat([init.fun(out3.xd); 0], 1, size(out3.xu_frames,2));


%%
h1 = figure;
plot(t0,x0(:,1), 'r-.', 'LineWidth',2);
hold on;
plot(t1, xu1(1,:), 'k-', 'LineWidth',2);
plot(t2, xu2(1,:), 'b--', 'LineWidth',2 );
plot(t3, xu3(1,:), 'g--', 'LineWidth',2 );
ylabel('z_{us} - z_g', 'FontSize',12);
xlabel('t', 'FontSize', 12)
legend('Open-loop (no control)', 'DT', 'LQR', 'Sequential');
set(get(h1,'CurrentAxes'), 'FontSize', 12);
grid on;

h2 = figure;
plot(t0,x0(:,3), 'r-.', 'LineWidth', 2);
hold on;
plot(t1, xu1(3,:), 'k-', 'LineWidth',2);
plot(t2, xu2(3,:), 'b--', 'LineWidth', 2)
plot(t3, xu3(3,:), 'g--', 'LineWidth',2 );
ylabel('z_s - z_{us}', 'FontSize', 12);
legend('Open-loop (no control)', 'DT', 'LQR',  'Sequential');
xlabel('t');
set(get(h2, 'CurrentAxes'), 'FontSize', 12);
grid on;

h3 = figure;
plot(t1, xu1(5,:), 'k-', 'LineWidth',2);
hold on;
plot(t2, xu2(5,:), 'b--','LineWidth',2);
plot(t3, xu3(5,:), 'g--','LineWidth',2);
legend('DT', 'LQR', 'Sequential');
xlabel('t', 'FontSize',12);
ylabel('u', 'FontSize', 12);
set(get(h3,'CurrentAxes'), 'FontSize', 12);
grid on;

