classdef SequentialSolver < handle
    % SEQUENTIALSOLVER A Sequential solver for MDO autonomous systems
    
    % Created by Zhi Han
    % This solver contains design constraints 
    % Revised by Tinghao 
    properties
        system
        objective
        initCondition
        
        % Design constraint is a constraint in terms of the design
        % variables
        designConstraint
        
        % State constraint is a constraint in terms of the state and
        % control variables
        stateConstraint
        options
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
            out.flag = out2.flag;
            out.step1 = out1;
            out.step2 = out2;
        end
    end
    %% Prototypes
    % These are function prototypes. They act like the declaration in C/C++.
    % The actual implemenation of the functions are found in this file with
    % suffix _impl.
    methods
        function [f, Gradient] = f_obj_design(obj, xd, u, xinit)
            % Objective for the plant design proble
            [f, Gradient] = f_obj_design_impl(obj, xd, u, xinit);
            
        end
        
        function [G, H, JG, JH] = f_con_design(obj, x)
            % Objective for the plant design proble
            
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

function yesno = myCheckSolver(obj)
if strcmp(optimget(obj.options, 'GradObj'), 'on') && ...
    isempty(obj.objective.gradient_d)
    error('SequentialSolver:NoDesignGradient', ...
        'Design gradient not provided in the user code');
end
end

function out = f_solve_xd_impl(obj)
checkSolver(obj);
myCheckSolver(obj);

% Solve optimization for design variable xd
obj.n_design = length(obj.xd0);
xinit = obj.initCondition.fun(obj.xd0);
obj.n_state = length(xinit);
t = obj.t;

options = obj.options(1);

u = zeros(obj.n_control,1);  % Assume no control input
xd0 = obj.xd0; %initial guess

[xdopt,fopt,flag,output]=fmincon(@(xd)obj.f_obj_design(xd, u, xinit), ...
    xd0,[],[],[],[], ...
    obj.lb,obj.ub,@(xd)f_con_design(obj, xd),options);

out.xd = xdopt;
out.fopt = fopt;
out.flag = flag;
out.output = output;
end



function [obj_f, Gradient] = f_obj_design_impl(obj, xd, u, xinit)
t = obj.t;
nt = length(t);
[~,x_state] = ode23(@(t, x)obj.system.deriv(t, x, u, xd) , ...
    t, xinit);
x_state_frames = reshape(x_state', obj.n_state, nt);
xu_frames = [x_state_frames; repmat(u, 1, nt)];
obj_f = obj.objective.fun(xd, xu_frames, obj.t);

if ~isempty(obj.objective.gradient_d)
    Gradient = obj.objective.gradient_d(xd,xu_frames, obj.t)';
else
    % No design gradient function
    Gradient = [];
end
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
xu0 = reshape([x_state'; zeros(obj.n_control, nt)], ...
    (obj.n_state+obj.n_control)*nt,1);

options = obj.options(2);
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

function [G, H, JG, JH] = f_con_design_impl(obj,xd)
% Evaluate the design constraints, and their gradient
H = []; % (XXX)
JH = [];

if ~isempty(obj.designConstraint)
    G = obj. designConstraint.fun(xd,[]);
    JG = obj. designConstraint.jacobian(xd,[])';
else
    G = [];
    JG = [];
end
end


function [G, H, JG, JH] = f_con_control_impl(obj, xd, xu)
% Direct-transcription evaluation of constraints
nt = length(obj.t);
xu_frames = reshape(xu, obj.n_state + obj.n_control, nt);
xinit = obj.initCondition.fun(xd);

dimensions = struct('n_state', obj.n_state, 'n_control', obj.n_control, ...
    'n_input', 0);

[H, JH] = defectEquation(obj.system, obj.initCondition, [], ...
                         obj.t, xu_frames, xd, dimensions);
J = JH';

JH = J(:, length(xd)+1:end)'; % extract the parts

if ~isempty(obj.designConstraint)
    
    g = obj.designConstraint.fun(xd, xu_frames(3,:));
    G = g(10:end);
    JG = obj.designConstraint.state_jacobian(xd, xu_frames(3,:))';
else
    G = []; %(XXX) to implement
    JG = [];
    
end
end
