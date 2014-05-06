function [obj, t, xu, k] = lqr_control_design(xd, system, objective, init, b, tf)
% Design an LQR 
K0 = [0;0;0;0];
Ja = system.jacobian(0, [],[],[xd;K0]);

% nd = length(xd);
x0 = init.fun([xd; K0]);   % Open-loop initial state
nx = length(x0);
% nu = size(J,2) - nx;
a = Ja(:,1:nx);

QR = objective.QR();
Q = QR(1:nx, 1:nx);
R = QR(nx+1:end, nx+1:end);
k = lqr(a,b,Q,R);

if numel(tf)==1
    tspan = [0 tf];
else
    tspan = tf;
end
[t,x] = ode23s(@(t,x)system.deriv(t,x,[], [xd;-k']), tspan, x0);
xu1 = x';
obj = objective.fun([xd;-k'], xu1,t);
xu2 = -k * xu1;

xu1 = xu1 + repmat(x0, 1, size(xu1,2));

xu = [xu1; xu2];


end