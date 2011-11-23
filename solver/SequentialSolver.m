classdef SequentialSolver < handle
    % SEQUENTIALSOLVER A Sequential solver for MDO autonomous systems
    
    % Created by Zhi Han
    properties
        system
        objective
        initCondition
        
        designConstraint
    end
    
    properties
        t
        pinit
        
        xd0 % initial guess of design variables
        
        x0
        n_control
        
        lb
        ub
        
        n_design
        n_state
    end
    
    methods
        function out = f_solve(obj)
            % Sequential solver
            fprintf('\n Step 1: Solve the design problem');
            out1 = obj.f_solve_xd();
            xd = out1.xd;
            fprintf('\n Step 2: Solve the control problem');
            out2 = obj.f_solve_xu(xd);
            
            out.xd = out1.xd;
            out.xu_frames = out2.xu;
            out.fopt = out2.fopt; 
            out.step1 = out1;
            out.step2 = out2;
        end
    end
    %% Prototypes
    % These are function prototypes. They act like the declaration in C/C++.
    % The actual implemenation of the functions are found in this file with
    % suffix _impl.
    methods
        function [obj, Gradient] = f_obj_design(obj, xd, u, xinit)
            % Objective for the plant design proble
            obj = f_obj_design_impl(obj, xd, u, xinit);
            Gradient = []; % (XXX) not implemented
        end
        
        function [G, H, JG, JH] = f_con_design(obj, x)
            % Objective for the plant design proble
            error('DT:SequentialSolver:NoImplementation',...
                'Not implemented yet');
            [G, H, JG, JH] = f_con_design_impl(obj, x);
        end
        function out = f_solve_xd(obj)
            % Solve the plant design problem
            out = f_solve_xd_impl(obj);
        end
        
        function [obj, Gradient] = f_obj_control(obj, xd, xu)
            % Objective for the control design problem
            [obj, Gradient] = f_obj_control_impl(obj,xd, xu);
        end
        function [G, H, JG, JH] = f_con_control(obj, xd, xu)
            % Constraints for the control design problem
            [G, H, JG, JH] = f_con_control_impl(obj,xd, xu);
        end
        function out = f_solve_xu(obj, xd)
            out = f_solve_xu_impl(obj, xd);
        end
        
    end
    
end



function out = f_solve_xd_impl(obj)
% Solve optimization for design variable xd
obj.n_design = length(obj.xd0);
xinit = obj.initCondition.fun(obj.xd0);
obj.n_state = length(xinit);
t = obj.t;

options = optimset('display','iter','algorithm','sqp',...
    'MaxIter', 1e+6, 'MaxFunEvals', 1e+6, ...
    'LargeScale', 'on', ...
    'TolFun',1e-4, 'UseParallel','always', ...
    'GradConstr','off', 'GradObj','off','DerivativeCheck','off');

u = zeros(obj.n_control,1);  % Assume no control input
xd0 = obj.xd0; %initial guess

[xdopt,fopt,flag,output]=fmincon(@(xd)obj.f_obj_design(xd, u, xinit), ...
    xd0,[],[],[],[], ...
    obj.lb,obj.ub,[],options);

out.xd = xdopt;
out.fopt = fopt;
out.flag = flag;
out.output = output;
end



function [obj_f] = f_obj_design_impl(obj, xd, u, xinit)
t = obj.t;
nt = length(t);
[~,x_state] = ode23(@(t, x)obj.system.deriv(t, x, u, xd) , ...
    t, xinit);
x_state_frames = reshape(x_state', obj.n_state, nt);
xu_frames = [x_state_frames; repmat(u, 1, nt)];
obj_f = obj.objective.fun(xd, xu_frames, obj.t);
end

function out = f_solve_xu_impl(obj, xd)
% Solve optimization for design variable xd
xinit = obj.initCondition.fun(obj.xd0);
t = obj.t;
nt = length(t);
% Simulate for initial state trajectory
[~,x_state] = ode23(...
    @(t,x)obj.system.deriv(t,x,zeros(obj.n_control,1), xd),...
    t, xinit);
xu0 = reshape([x_state'; ones(obj.n_control, nt)], (obj.n_state+obj.n_control)*nt,1);
options = optimset('display','iter','algorithm','interior-point',...
    'MaxIter', 1e+6, 'MaxFunEvals', 1e+6, ...
    'LargeScale', 'on', ...
    'TolFun',1e-6, 'UseParallel','never', ...
    'GradConstr','on', 'GradObj','on','DerivativeCheck','off');
ulb = -2000; % Some arbitrary lower bound and upper bound
uub = 2000;
[xuopt,fopt,flag,output]=fmincon(@(xu)obj.f_obj_control(xd, xu), ...
    xu0,[],[],[],[], ...
    ulb,uub,@(xu)obj.f_con_control(xd, xu),options);

out.xu = reshape(xuopt, obj.n_control+obj.n_state, nt);
out.fopt = fopt;
out.flag = flag;
out.output = output;
end


function [obj_f, Gradient] = f_obj_control_impl(obj, xd, xu)
% Direct-transcription evaluation of control objective
nt = length(obj.t);
x_state_frames = reshape(xu, obj.n_state + obj.n_control, nt);

obj_control = obj.objective.fun(xd, x_state_frames, obj.t);
obj_control_Gradient = obj.objective.gradient(xd, x_state_frames, obj.t);

obj_f = obj_control;

Gradient = obj_control_Gradient;
end

function [G, H, JG, JH] = f_con_control_impl(obj, xd, xu)
% Direct-transcription evaluation of constraints
nt = length(obj.t);
xu_frames = reshape(xu, obj.n_state + obj.n_control, nt);
xinit = obj.initCondition.fun(xd);
H0 = xinit - xu_frames(1:obj.n_state, 1);

H_state = zeros(obj.n_state, nt);
JH_state = sparse(obj.n_state*nt, ...
    (obj.n_state + obj.n_control)*nt); % pre-allocate memory

for i=2:nt
    % Calculate the defect constraints
    x_i = xu_frames(1:obj.n_state,i);
    u_i = xu_frames(obj.n_state+1:obj.n_state+obj.n_control, i);
    xd_i = obj.system.deriv(obj.t(i), x_i, u_i, xd);
    
    x_i_1 = xu_frames(1:obj.n_state,i-1);
    u_i_1 = xu_frames(obj.n_state+1:obj.n_state+obj.n_control, i-1);
    xd_i_1 = obj.system.deriv(obj.t(i-1), x_i_1, u_i_1, xd);
    
    dti = obj.t(i) - obj.t(i-1);
    
    H_state(:,i) = x_i - x_i_1 - dti/2 * (xd_i + xd_i_1);

    dfdx_i = obj.system.jacobian(obj.t(i), x_i, u_i, xd);
    dfdx_i_1 = obj.system.jacobian(obj.t(i-1), x_i_1, u_i_1, xd);
    
    JH_state((i-1)*obj.n_state+1:i*obj.n_state, ...
        (i-1)*(obj.n_state+obj.n_control)+1:i*(obj.n_state +obj.n_control)) ...
        = [speye(obj.n_state), zeros(obj.n_state, obj.n_control)] - dti/2*dfdx_i;
    JH_state((i-1)*obj.n_state+1:i*obj.n_state, ...
        (i-2)*(obj.n_state + obj.n_control) +1:(i-1)*(obj.n_state + obj.n_control)) ...
        = [-speye(obj.n_state), zeros(obj.n_state, obj.n_control)] - dti/2*dfdx_i_1;
end

H = [reshape(H_state, numel(H_state), 1); H0];

JH = [JH_state; ...
      -speye(obj.n_state), zeros(obj.n_state,obj.n_control), ...
    zeros(obj.n_state, ...
    (nt-1)*(obj.n_state+obj.n_control))]'; % (XXX) to implement
G = []; %(XXX) to implement
JG = [];
end
