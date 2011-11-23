classdef DTSolver <handle
    % DTSOLVER Direct transcription solver for autonomous systems (no inputs)
    %
    
    % Created by Zhi Han
    
    properties
        system      %An ODE system
        objective   %Control objective
        initCondition
        
        designConstraint % Design constraints
    end
    
    properties
        t   % Vector of time
        
        xd0% % Initial guess of design variables
        
        x0 % initial states
        n_control % number of control variables
        
        lb  % lower bound
        ub  % upper bound
        
        n_design  % number of parameters
        n_state  % number of states
    end
    
    methods
        function [obj, Gradient] = f_obj(obj, x)
            [obj, Gradient] = f_obj_impl(obj,x);
        end
        function [G, H, JG, JH] = f_con(obj, x)
            [G, H, JG, JH] = f_con_impl(obj, x);
        end
        function [out] = f_solve(obj)
            out = f_solve_impl(obj);
        end
        function h = f_hessian(obj, x, lambda)
            h = f_hessian_impl(obj, x, lambda);
        end
    end
end

function [out] = f_solve_impl(obj)
obj.n_design = length(obj.xd0);
xinit = obj.initCondition.fun(obj.xd0);
obj.n_state = length(xinit);
t = obj.t;
% Initial control is zero
nt = length(t);
if isempty(obj.x0)
    % Step 1,
    t = obj.t;
    [~,x_state] = ode23(...
        @(t,x)obj.system.deriv(t,x,zeros(obj.n_control,1),obj.xd0),...
        t, xinit);
    nt = length(t);
    u0 = ones(obj.n_control, length(t));
    xu = reshape([x_state'; u0], (obj.n_state + obj.n_control) * nt, 1);
    x0 = [obj.xd0; xu];
else
    x0 = obj.x0;
end

options = optimset('display','iter','algorithm','interior-point',...
    'MaxIter', 1e+6, 'MaxFunEvals', 1e+6, ...
    'LargeScale', 'on', ...
    'TolFun',1e-6, 'UseParallel','always', ...
    'GradConstr','on', 'GradObj','on', 'DerivativeCheck','off');
if ~isempty(obj.objective.hessian)
    options = optimset(options, 'Hessian', 'on', ...
        'HessFcn', @(x, lambda) obj.f_hessian(x, lambda));
end
[xopt,fopt,flag,output]=fmincon(@(x)obj.f_obj(x),x0,[],[],[],[], ...
    obj.lb,obj.ub,@(x)obj.f_con(x),options);
out.xopt = xopt;
out.xd = xopt(1:obj.n_design);
out.xu    = reshape(xopt(obj.n_design+1:end), (obj.n_state+obj.n_control), nt);
out.fopt = fopt;
out.flag = flag;
out.output = output;

end

function H = f_hessian_impl(obj, x, lambda)
x_state = x(obj.n_design+1:end);
x_d = x(1:obj.n_design);
nt = length(x_state)/(obj.n_state + obj.n_control);
x_state_frames = reshape(x_state, obj.n_state + obj.n_control, nt);

G = obj.designConstraint.fun(x_d);
m = size(G,1);
F = obj.designConstraint.hessian(x_d, lambda.ineqnonlin(1:m));
H =  [F, sparse(obj.n_design,length(x_state));
    sparse(length(x_state), obj.n_design),  ...
    obj.objective.hessian(x_d, x_state_frames, obj.t)];
end

function [obj_f, Gradient, H] = f_obj_impl(obj ,x)
x_state = x(obj.n_design+1:end);
x_d = x(1:obj.n_design);
nt = length(x_state)/(obj.n_state + obj.n_control);
x_state_frames = reshape(x_state, obj.n_state + obj.n_control, nt);
obj_control = obj.objective.fun(x_d, x_state_frames, obj.t);
obj_control_Gradient = obj.objective.gradient(x_d, x_state_frames, obj.t);

obj_f = obj_control;
if isempty(obj.objective.fun_d)
    gradient_d = zeros(obj.n_design,1);
else
    gradient_d = obj.objective.gradient_d(x_d, x_state_frames, obj.t);
end
Gradient = [gradient_d; obj_control_Gradient];
end

function [G, H, JG, JH] = f_con_impl(obj,x)
x_state = x(obj.n_design+1:end);
nt = length(x_state)/(obj.n_state + obj.n_control);
x_state_frames = reshape(x_state, obj.n_state + obj.n_control, nt);

% Calculate h and Jh
params = x(1:obj.n_design);
xinit = obj.initCondition.fun(params);
H0 = xinit - x_state_frames(1:obj.n_state,1);

H_state = zeros(obj.n_state, nt);
JH_state = sparse(obj.n_state*nt, ...
    (obj.n_state + obj.n_control)*nt);
for i=2:nt
    state_i = x_state_frames(1:obj.n_state,i);
    control_i = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control,i);
    
    xd_i = obj.system.deriv(obj.t(i), state_i, control_i, params);
    state_i_1 = x_state_frames(1:obj.n_state,i-1);
    control_i_1 = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control,i-1);
    
    xd_i_1 = obj.system.deriv(obj.t(i-1) ,state_i_1, control_i_1, params);
    dti = obj.t(i) - obj.t(i-1);
    H_state(:,i) = state_i - state_i_1 ...
        - dti/2 * (xd_i + xd_i_1);
    
    dfdx_i = obj.system.jacobian(obj.t(i), state_i, control_i, params);
    dfdx_i_1 = obj.system.jacobian(obj.t(i-1), state_i_1, control_i_1,params);
    JH_state((i-1)*obj.n_state+1:i*obj.n_state, ...
        (i-1)*(obj.n_state+obj.n_control)+1:i*(obj.n_state +obj.n_control)) ...
        = [speye(obj.n_state), zeros(obj.n_state, obj.n_control)] - dti/2*dfdx_i;
    JH_state((i-1)*obj.n_state+1:i*obj.n_state, ...
        (i-2)*(obj.n_state + obj.n_control) +1:(i-1)*(obj.n_state + obj.n_control)) ...
        = [-speye(obj.n_state), zeros(obj.n_state, obj.n_control)] - dti/2*dfdx_i_1;
end

JH_param = sparse(obj.n_state*nt, obj.n_design);
for i=2:nt
    state_i = x_state_frames(1:obj.n_state,i);
    control_i = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control,i);
    state_i_1 = x_state_frames(1:obj.n_state,i-1);
    control_i_1 = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control,i-1);
    dti = obj.t(i) - obj.t(i-1);
    if ~isempty(obj.system.jacobian_d)
        % Has non-empty design variable Jacobian, i.e.,
        % Derivative explicitly depends on the design variable
        dfdpi = obj.system.jacobian_d(obj.t(i), state_i, control_i, params);
        dfdpi_1 = obj.system.jacobian_d(obj.t(i-1), state_i_1, control_i_1,params);
        JH_param((i-1)*obj.n_state +1: i*obj.n_state, :) = -dti/2* (dfdpi + dfdpi_1);
    else
        % Derivative does not depend on design variables
        % These entries are left to be zeros.
    end
end

JH = [JH_param, JH_state; ...
    obj.initCondition.jacobian(params), ...
    -speye(obj.n_state), zeros(obj.n_state,obj.n_control), ...
    zeros(obj.n_state, ...
    (nt-1)*(obj.n_state+obj.n_control))]';
H = [reshape(H_state, numel(H_state), 1); H0];


% Add design constraints
if ~isempty(obj.designConstraint)
    G = obj.designConstraint.fun(params);
    JG = [obj.designConstraint.jacobian(params), zeros(size(G,1), length(x_state))]';
else
    G = [];
    JG = [];
end
end



