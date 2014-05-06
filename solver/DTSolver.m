classdef DTSolver <handle
    % DTSOLVER Direct transcription solver for autonomous systems (no inputs)
    %
    
    % Created by Zhi Han
    
    properties
        system      %An ODE system
        objective   %Control objective
        initCondition
        
        options  % optimization options
        
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
checkSolver(obj);

obj.n_design = length(obj.xd0);
xinit = obj.initCondition.fun(obj.xd0);
obj.n_state = length(xinit);
t = obj.t;
% Initial control is zero
nt = length(t);
if isempty(obj.x0)
    % Solve the ODe to get the time points
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

options = obj.options;

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

function [obj_f, Gradient] = f_obj_impl(obj ,x)
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
dimensions = struct('n_state', obj.n_state, 'n_control', obj.n_control, ...
    'n_input', 0);
[H, JH] = defectEquation(obj.system, obj.initCondition, [], ...
    obj.t, x_state_frames, params, dimensions);

% Add design constraints
if ~isempty(obj.designConstraint)
    G = obj.designConstraint.fun(params);
    JG = [obj.designConstraint.jacobian(params), zeros(size(G,1), length(x_state))]';
else
    G = [];
    JG = [];
end
end



