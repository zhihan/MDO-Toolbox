classdef NestedSolver
    properties
        system      % An ODE system
        objective   % Control objective
        initCondition
    end
    
    properties
        h    % step size
        tf   % final time

        pinit % Initial guess of parameter
        
        lb  % lower bound
        ub  % upper bound
        
        n_param  % number of parameters
        n_control % number of control signals
        
    end
    methods
        function [obj] = f_obj(obj, x)
            [obj] = f_obj_impl(obj,x);
        end
        function [G, H] = f_con(obj, x)
            [G, H] = f_con_impl(obj, x);
        end
        function [out] = f_solve(obj)
            out = f_solve_impl(obj);
        end
    end    
end

function [G, H] = f_con_impl(obj, x)
G = [];
H = [];
end

function [obj_f] = f_obj_impl(obj, x)
    params = x;
    xinit = obj.initCondition.fun(params);
    t = 0:obj.h:obj.tf;
    n_state = length(xinit);
    nt = length(t);
    [~,x_state] = ode23(obj.system.deriv, t, xinit);
    x_state_frames = x_state';
    
    obj_f = obj.objective.fun(x_state_frames);
     
end

function out = f_solve_impl(obj)
    obj.n_param = length(obj.pinit);
    options = optimset('display','iter','algorithm','interior-point',...
                       'maxfunevals',20000,'DiffMinChange',1e-7,...
                       'UseParallel','never', ...
                       'GradConstr','off', 'GradObj','off'); %, 'DerivativeCheck','on');

    [xopt,fopt,flag,output]=fmincon(@(x)obj.f_obj(x), obj.pinit,[],[],[],[], ...
        obj.lb,obj.ub,@(x)obj.f_con(x),options);
    out.xopt = xopt;
    out.param = xopt;
    out.fopt = fopt;
    out.flag = flag;
    out.output = output;

end