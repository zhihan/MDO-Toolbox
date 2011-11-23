classdef DTSolver2 < handle
%  DTSolver2 (experimental) direct transcription solver with tFinal constraints
%
    
% DTSolver2 was a temporary solution for formulating the final time 
% as a constraint rather than a fixed value. See demo golfball/demo_golfball.m
% for further demonstration.
%

% Created by Zhi Han
    
    properties
        system      %An ODE system
        objective   %Control objective
        initCondition
        stateConstraint
    end
    
    properties
        nt_1    % number of steps
        tf    % hit time
        pinit % Initial guess of parameter (theta)
        
        x0 % initial states
        
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
    end
end

function [out] = f_solve_impl(obj)
    obj.n_param = length(obj.pinit);
    xinit = obj.initCondition.fun(obj.pinit);
    obj.n_state = length(xinit);
    
    h = obj.tf/obj.nt_1;
    t = 0:h:obj.tf;
    nt = obj.nt_1 + 1;
    if isempty(obj.x0) 
% Step 1,
        [~,x_state] = ode23(obj.system.deriv, t, xinit);
        nt = length(t);
        x0 = [obj.pinit; obj.tf; reshape(x_state', obj.n_state * nt, 1)];
    else
        x0 = obj.x0;
    end
    
    options = optimset('display','iter','algorithm','interior-point',...
                       'maxfunevals',20000,'DiffMinChange',1e-7,'UseParallel','never', ...
                       'GradConstr','on', 'GradObj','on');%, 'DerivativeCheck','on');
    
    [xopt,fopt,flag,output]=fmincon(@(x)obj.f_obj(x),x0,[],[],[],[], ...
        obj.lb,obj.ub,@(x)obj.f_con(x),options);
    out.xopt = xopt;
    out.param = xopt(1:obj.n_param);
    out.x    = reshape(xopt(obj.n_param+2:end), obj.n_state, nt);
    out.fopt = fopt;
    out.flag = flag;
    out.output = output;
    
end

function [obj_f, Gradient] = f_obj_impl(obj ,x)
    x_state = x(obj.n_param+2:end);
    nt = length(x_state)/obj.n_state;
    x_state_frames = reshape(x_state, obj.n_state, nt);
    
    obj_f = obj.objective.fun(x_state_frames);
    obj_control_Gradient = obj.objective.gradient(x_state_frames);
    Gradient = [zeros(obj.n_param,1); ...
        0;
        obj_control_Gradient];     
end

function [G, H, JG, JH] = f_con_impl(obj,x)
    x_state = x(obj.n_param+2:end);
    nt = length(x_state)/obj.n_state;
    tf = x(obj.n_param+1);
    h = tf/obj.nt_1;
    x_state_frames = reshape(x_state, obj.n_state, nt);

    % Calculate h and Jh
    params = x(1:obj.n_param);
    xinit = obj.initCondition.fun(params);
    H0 = xinit - x_state_frames(:,1);
    
    H_state = zeros(obj.n_state, nt);
    JH_state = sparse(obj.n_state*nt, obj.n_state*nt);
    JH_tf_frame = sparse(obj.n_state, nt);
    for i=2:nt
        xdi = obj.system.deriv(i*h, x_state_frames(:,i));
        xdi_1 = obj.system.deriv((i-1)*h ,x_state_frames(:,i-1));
        H_state(:,i) = x_state_frames(:,i) - x_state_frames(:,i-1) ...
            - h/2 * (xdi + xdi_1);
        JH_tf_frame(:,i) = - 1/2/obj.nt_1*(xdi + xdi_1);
        dfdxi = obj.system.jacobian(i*h, x_state_frames(:,i));
        dfdxi_1 = obj.system.jacobian((i-1)*h, x_state_frames(:,i-1));
        JH_state((i-1)*obj.n_state+1:i*obj.n_state, (i-1)*obj.n_state+1:i*obj.n_state) ...
            = speye(obj.n_state) - h/2*dfdxi;
        JH_state((i-1)*obj.n_state+1:i*obj.n_state, (i-2)*obj.n_state+1:(i-1)*obj.n_state) ...
            = -speye(obj.n_state) - h/2*dfdxi_1;
    end
    
    Hfinal_state = obj.stateConstraint.eqFun(x_state_frames);
    JHfinal_state = obj.stateConstraint.eqJacobian(x_state_frames);
    JH_tf = reshape(JH_tf_frame,obj.n_state*nt,1);
    
    JH = [zeros(obj.n_state*nt,1), JH_tf, JH_state; ...
         obj.initCondition.jacobian(params), zeros(obj.n_state,1)...
         -speye(obj.n_state), zeros(obj.n_state, ...
         (nt-1)*obj.n_state);
         zeros(1,1), zeros(1,1), JHfinal_state ]';
    H = [reshape(H_state, numel(H_state), 1); H0; Hfinal_state];
    
    
    G = []; 
    JG = [];
end



