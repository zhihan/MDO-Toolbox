classdef DTSolverWithInputs <handle
% DTSOLVERWITHINPUTS Direct transcription solver for systems with open input channel
% The optimization can be formulated for multiple simulation runs. 
% 

% See qcar_spring_damper_multi_codesign for sample code of using this solver.
%
    
% Created by Zhi Han
    
    properties
        system      %An ODE system
        objective   %Control objective
        initCondition
        
        designConstraint % Design constraints
    end
    
    properties
        t  % Cell array
        input % cell array of function handles
        weight % cell array of weights
        
        pinit % Initial guess of parameter
        
        x0 % initial states
        n_control % number of control variables
        
        lb  % lower bound
        ub  % upper bound
        
        n_param  % number of parameters
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
    obj.n_param = length(obj.pinit);
    xinit = obj.initCondition.fun(obj.pinit);
    obj.n_state = length(xinit);
    
    t = obj.t;
    % Initial control is zero
    if isempty(obj.x0)
        if iscell(obj.t)
            xu = [];
            for i=1:length(obj.t)
                t = obj.t{i};
                input = obj.input{i};
                xu = [xu; initial_sim(obj, t, input)];
            end
        else
            xu = initial_sim(obj, t);
        end
        x0 = [obj.pinit; xu];
    else
        x0 = obj.x0;
    end
    
    options = optimset('display','iter','algorithm','interior-point',...
                       'MaxIter', 5e+3, 'MaxFunEvals', 1e+6, ...
                       'LargeScale', 'on', ...
                       'TolFun',1e-4, 'UseParallel','always', ...
                       'GradConstr','on', 'GradObj','on', 'DerivativeCheck','off');
   if ~isempty(obj.objective.hessian)
       options = optimset(options, 'Hessian', 'on', ...
           'HessFcn', @(x, lambda) obj.f_hessian(x, lambda));
   end
    [xopt,fopt,flag,output]=fmincon(@(x)obj.f_obj(x),x0,[],[],[],[], ...
        obj.lb,obj.ub,@(x)obj.f_con(x),options);
    out.xopt = xopt;
    out.param = xopt(1:obj.n_param);
    offset = 0;
    for i=1: length(obj.t)
        nt = length(obj.t{i});
        total_length = (obj.n_state+obj.n_control )* nt;
        out.xu{i} = reshape(xopt(obj.n_param + offset +1: obj.n_param+offset+total_length), ...
            (obj.n_state+obj.n_control), nt);
        offset = offset + total_length;
    end
    out.fopt = fopt;
    out.flag = flag;
    out.output = output;
    
end
function xu = initial_sim(obj,t, input)
xinit = obj.initCondition.fun(obj.pinit);
[~,x_state] = ode23(...
    @(t,x)obj.system.deriv(t,x,zeros(obj.n_control,1),obj.pinit, input),...
    t, xinit);
nt = length(t);
u0 = zeros(obj.n_control, length(t));
xu = reshape([x_state'; u0], (obj.n_state + obj.n_control) * nt, 1);
end
function H = f_hessian_impl(obj, x, lambda)
x_state = x(obj.n_param+1:end);
x_d = x(1:obj.n_param);
G = obj.designConstraint.fun(x_d);
m = size(G,1);
F = obj.designConstraint.hessian(x_d, lambda.ineqnonlin(1:m));

state_Hess = sparse(length(x_state), length(x_state));
offset = 0;
for i=1:length(obj.t)
    total_state = length(obj.t{i}) * (obj.n_state + obj.n_control);
    nt = length(obj.t{i});
    w = obj.weight{i};
        
    x_state_frames = reshape(x_state(offset+1:offset+total_state), ...
        obj.n_state + obj.n_control, nt);
    state_Hess(offset+1:offset+total_state, offset+1:offset+total_state) = ...
        w*obj.objective.hessian(x_d, x_state_frames, obj.t{i});
    offset = offset + total_state;
end
H =  [F, sparse(obj.n_param,length(x_state));
    sparse(length(x_state), obj.n_param),  ...
    state_Hess];
end

function [total_f, Gradient] = f_obj_impl(obj ,x)
    x_d = x(1:obj.n_param);
    offset = 0;
    total_control_Gradient = [];
    total_f = 0;
    total_gradient_d = zeros(obj.n_param,1);
    for i=1:length(obj.t)
        t = obj.t{i};
        nt = length(t);
        w = obj.weight{i};
        total_length = (obj.n_state+obj.n_control)*nt;
        
        x_state = x(obj.n_param+ offset +1: offset + obj.n_param+total_length);
        offset = offset + total_length;
        
        x_state_frames = reshape(x_state, obj.n_state + obj.n_control, nt);
        obj_control = w*obj.objective.fun(x_d, x_state_frames, t);
        obj_control_Gradient = w*obj.objective.gradient(x_d, x_state_frames, t);
        obj_f = obj_control;
        
        total_control_Gradient = [total_control_Gradient; obj_control_Gradient];
        total_f = total_f + obj_f;
    
        
        if isempty(obj.objective.fun_d)
            gradient_d = zeros(obj.n_param,1);
        else
            gradient_d = w*obj.objective.gradient_d(x_d, x_state_frames, t);
        end
        total_gradient_d = total_gradient_d + gradient_d;
    end
    

    Gradient = [gradient_d; total_control_Gradient];
end

function [G, H, JG, JH] = f_con_impl(obj,x)
    x_d = x(1:obj.n_param);
    T = [];
    JH_state = [];
    H_state = [];
    
    % calculate total size
    NT = 0;
    M = 0;
    for i=1:length(obj.t)
        t = obj.t{i};
        nt = length(t);
        NT = NT + nt;
        M = M + nt +1; 
    end
    
    H_state = sparse(obj.n_state, M);
    JH_state = sparse(obj.n_state*M, ...
            (obj.n_state + obj.n_control)*NT);
    JH_param = sparse(obj.n_state*M, obj.n_param);
    offset = 0;
    m_offset = 0;
    n_offset = 0;
    for i=1:length(obj.t)
        t = obj.t{i};
        input = obj.input{i};
        nt = length(t);
        total_length = nt * (obj.n_state + obj.n_control);
        
        x_state = x(obj.n_param+offset+1:obj.n_param+offset + total_length);
        offset = offset + total_length;
        x_state_frames = reshape(x_state, obj.n_state + obj.n_control, nt);
        
        % Calculate h and Jh
        params = x(1:obj.n_param);
        xinit = obj.initCondition.fun(params);
        H0 = xinit - x_state_frames(1:obj.n_state,1);
        
        for i=2:nt
            state_i = x_state_frames(1:obj.n_state,i);
            control_i = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control,i);
            
            xdi = obj.system.deriv(t(i), state_i, control_i, params, input);
            state_i_1 = x_state_frames(1:obj.n_state,i-1);
            control_i_1 = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control,i-1);
            
            xdi_1 = obj.system.deriv(t(i-1) ,state_i_1, control_i_1, params, input);
            dti = t(i) - t(i-1);
            H_state(:,m_offset/obj.n_state +i) = state_i - state_i_1 ...
                - dti/2 * (xdi + xdi_1);
            
            dfdxi = obj.system.jacobian(t(i), state_i, control_i, params);
            dfdxi_1 = obj.system.jacobian(t(i-1), state_i_1, control_i_1,params);
            JH_state(m_offset + (i-1)*obj.n_state+1:m_offset+i*obj.n_state, ...
                n_offset+(i-1)*(obj.n_state+obj.n_control)+1:n_offset+i*(obj.n_state +obj.n_control)) ...
                = [speye(obj.n_state), zeros(obj.n_state, obj.n_control)] - dti/2*dfdxi;
            JH_state(m_offset + (i-1)*obj.n_state+1:m_offset + i*obj.n_state, ...
                n_offset+(i-2)*(obj.n_state + obj.n_control) +1:n_offset+(i-1)*(obj.n_state + obj.n_control)) ...
                = [-speye(obj.n_state), zeros(obj.n_state, obj.n_control)] - dti/2*dfdxi_1;
        end

        for i=2:nt
            state_i = x_state_frames(1:obj.n_state,i);
            control_i = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control,i);
            state_i_1 = x_state_frames(1:obj.n_state,i-1);
            control_i_1 = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control,i-1);
            dti = t(i) - t(i-1);
            if ~isempty(obj.system.jacobian_d)
                % Has non-empty design variable Jacobian, i.e.,
                % Derivative explicitly depends on the design variable
                dfdpi = obj.system.jacobian_d(t(i), state_i, control_i, params);
                dfdpi_1 = obj.system.jacobian_d(t(i-1), state_i_1, control_i_1,params);
                JH_param(m_offset + (i-1)*obj.n_state +1: m_offset+i*obj.n_state, :) = -dti/2* (dfdpi + dfdpi_1);
            else
                % Derivative does not depend on design variables
                % These entries are left to be zeros.
            end
        end
        JH_param(m_offset + i*obj.n_state+1: m_offset+(i+1)*obj.n_state, :) = obj.initCondition.jacobian(params);
        JH_state(m_offset + i*obj.n_state+1:m_offset+ (i+1)*obj.n_state, ...
            n_offset +1: n_offset + obj.n_state) = -speye(obj.n_state);
        
        H_state(:, m_offset/obj.n_state + nt +1) = H0;
        m_offset = m_offset + (nt + 1)*obj.n_state;
        n_offset = n_offset + nt*(obj.n_state + obj.n_control);
    end
    
    H = reshape(H_state, numel(H_state),1);
    JH = [JH_param, JH_state]';
    % Add design constraints
    if ~isempty(obj.designConstraint)
        G = obj.designConstraint.fun(params);
        JG = [obj.designConstraint.jacobian(params), zeros(size(G,1), length(x)-obj.n_param)]';
    else
        G = [];
        JG = [];
    end
end



