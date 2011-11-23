function J  = perturb_x(fun, x0, varargin)
% Finite difference

if nargin < 5
    pert = 1e-5;
else
    pert = varargin{1};
end
f0 = fun(x0);
m = length(f0);
n = length(x0);
% Simple single-sided perturbation
J = zeros(m,n);
for i=1: n
    if length(pert) >1
        delta = pert(i);
    else
        delta = pert;
    end
    x = x0;
    x(i) = x(i) + delta;
    J(:,i) = (fun(x) - f0)/delta;
end
end