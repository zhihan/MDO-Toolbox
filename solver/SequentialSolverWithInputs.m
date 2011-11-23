classdef SequentialSolverWithInputs < handle
    % SEQUENTIALSOLVER A Sequential solver for MDO autonomous systems
    
    % Created by Zhi Han
    properties
        system
        objective
        initCondition
        
        designConstraint
    end
    
    properties
        t  %Cell array
        
        input %Cell array of function handles
        weight %Cell array of weights
        pinit
        
        xd0 % initial guess of design variables
        
        x0 %initial states
        n_control % number of control variables
        
        lb %lower bound
        ub %upper bound
        
        n_design % number of parameters
        n_state  % number of states
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
            
        end
        
        function [G, H, JG, JH] = f_con_design(obj, xd)
            % Objective for the plant design problem
            
            [G, H, JG, JH] = f_con_design_impl(obj, xd);
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
    obj.lb,obj.ub,@(xd)obj.f_con_design(xd),options);

out.xd = xdopt;
out.fopt = fopt;
out.flag = flag;
out.output = output;
end



function [total_f] = f_obj_design_impl(obj, xd, u, xinit)

obj.n_state=length(xinit);
x_state = [];
total_f = 0;
total_length = 0;
offset = 0;


for i=1:length(obj.t)
    t = obj.t{i};
    nt = length(t);
    w = obj.weight{i};
    input = obj.input{i};
    [~,x_state_sim] = ode23(@(t, x)obj.system.deriv(t, x, u, xd, input) , ...
    t, xinit);
    x_state = [x_state; x_state_sim];
    total_length = obj.n_state*nt;
    
    x_state_frames = reshape(x_state(offset/obj.n_state+1:offset/obj.n_state...
        +total_length/obj.n_state,:)', obj.n_state, nt);
    
    offset = offset+total_length;
    xu_frames = [x_state_frames; repmat(u, 1 ,nt)];
    obj_f = w*obj.objective.fun(xd, xu_frames, t);
    total_f = total_f + obj_f;
end  

end


function out = f_solve_xu_impl(obj, xd)
% Solve optimization for design variable xd
t = obj.t;
obj.n_design = length(obj.xd0);
xinit = obj.initCondition.fun(obj.xd0);
obj.n_state = length(xinit);
xu0 = [];

for i = 1:length(obj.t)
    t = obj.t{i};
    input = obj.input{i};
    xu0 = [xu0; initial_sim(obj, t, xd, input)];
end
options = optimset('display','iter','algorithm','interior-point',...
    'MaxIter', 1e+6, 'MaxFunEvals', 1e+6, ...
    'LargeScale', 'on', ...
    'TolFun',1e-6, 'UseParallel','never', ...
    'GradConstr','on', 'GradObj','on','DerivativeCheck','off');
ulb = -2000; % Some arbitrary lower bound and upper bound
uub = 2000;
% Simulate for initial state trajectory
[xuopt,fopt,flag,output]=fmincon(@(xu)obj.f_obj_control(xd, xu), ...
    xu0,[],[],[],[], ...
    ulb,uub,@(xu)obj.f_con_control(xd, xu),options);

offset = 0;
total_length = 0;

for i = 1:length(obj.t)
    nt = length(obj.t{i});
    total_length = (obj.n_state+obj.n_control)*nt;
    out.xu{i} = reshape(xuopt(offset+1: offset+total_length),...
        (obj.n_state+obj.n_control),nt);
    offset = offset + total_length;
end
out.fopt = fopt;
out.flag = flag;
out.output = output;
end

function xu0 = initial_sim(obj, t, xd, input)
xinit = obj.initCondition.fun(obj.xd0);
obj.n_state = length(xinit);
[~,x_state] = ode23(...
    @(t,x)obj.system.deriv(t,x,zeros(obj.n_control,1), xd, input),...
    t, xinit);
nt = length(t);
xu0 = reshape([x_state'; ones(obj.n_control, nt)], (obj.n_state+obj.n_control)*nt,1);
end



function [total_f, Gradient] = f_obj_control_impl(obj, xd, xu)
% Direct-transcription evaluation of control objective
t=obj.t;
offset = 0;
total_length = 0;
total_control_Gradient = [];
total_f = 0;


for i=1:length(obj.t)
    t = obj.t{i};
    nt = length(t);
    w = obj.weight{i};
    total_length = (obj.n_state+obj.n_control) * nt;
    x_state = xu(offset+1: offset+total_length);
    offset = offset + total_length;
    
    x_state_frames = reshape(x_state, obj.n_state + obj.n_control, nt);
    
    obj_control = w *obj.objective.fun(xd,x_state_frames,t); 
    obj_control_Gradient = w * obj.objective.gradient(xd, x_state_frames, t);
    obj_f = obj_control;
    total_control_Gradient = [total_control_Gradient; obj_control_Gradient];
    total_f = total_f + obj_f;   
    Gradient = total_control_Gradient;
end

end

function [G, H, JG, JH] = f_con_design_impl(obj,xd)
% Add design constraints
H = []; % (XXX)
JH = [];

if ~isempty(obj.designConstraint)
    G = obj. designConstraint.fun(xd);
    JG = [obj. designConstraint.jacobian(xd)];
else
    G = [];
    JG = [];
end
end


function [G, H, JG, JH] = f_con_control_impl(obj, xd, xu)
% Direct-transcription evaluation of constraints

T = [];
JH_state = [];
H_state = [];


% calculate total size

NT = 0;
M = 0;
for i = length (obj.t);
    t = obj.t{i};
    nt = length (t);
    NT = NT + nt;
    M = M + nt +1; 
end


H_state = sparse(obj.n_state, M);
JH_state = sparse(obj.n_state*M,...
    (obj.n_state + obj.n_control)*NT); % pre-allocate memory
offset = 0;
m_offset = 0;
n_offset = 0;
total_length = 0;

for j = 1: length(obj.t)
    t = obj.t{j};
    input = obj.input{j};
    nt = length(t);
    total_length = nt * (obj.n_state+ obj.n_control);
    
    x_state = xu(offset+1:offset+total_length);
    offset = offset + total_length;
    
    x_state_frames = reshape(x_state, obj.n_state + obj.n_control, nt);
    
    % calculate h and Jh
    xinit = obj.initCondition.fun(xd);
    H0 = xinit - x_state_frames(1:obj.n_state, 1);
    for i = 2:nt
        %calculate the defect contraints
        
        x_i = x_state_frames(1:obj.n_state,i);
        u_i = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control, i);
        xd_i = obj.system.deriv(t(i), x_i, u_i, xd, input);
    
        x_i_1 = x_state_frames(1:obj.n_state,i-1);
        u_i_1 = x_state_frames(obj.n_state+1:obj.n_state+obj.n_control, i-1);
        xd_i_1 = obj.system.deriv(t(i-1), x_i_1, u_i_1, xd, input);
        
        dti = t(i) - t(i-1);
    
        H_state(:,m_offset/obj.n_state+i) = x_i - x_i_1 ...
            - dti/2 * (xd_i + xd_i_1);
        
        dfdx_i = obj.system.jacobian(t(i), x_i, u_i, xd);
        dfdx_i_1 = obj.system.jacobian(t(i-1), x_i_1, u_i_1, xd);
        
        JH_state(m_offset + (i-1)*obj.n_state+1:m_offset + i*obj.n_state, ...
        n_offset+(i-1)*(obj.n_state+obj.n_control)+1:n_offset+i*(obj.n_state +obj.n_control)) ...
        = [speye(obj.n_state), zeros(obj.n_state, obj.n_control)] - dti/2*dfdx_i;
        JH_state(m_offset + (i-1)*obj.n_state+1:m_offset + i*obj.n_state, ...
        n_offset + (i-2)*(obj.n_state + obj.n_control) +1:n_offset +(i-1)*(obj.n_state + obj.n_control)) ...
        = [-speye(obj.n_state), zeros(obj.n_state, obj.n_control)] - dti/2*dfdx_i_1;
        
        
    end
    JH_state(m_offset + i * obj.n_state+1: m_offset + (i+1)*obj.n_state,...
        n_offset + 1: n_offset + obj.n_state) = -speye(obj.n_state);
    H_state(:, m_offset/obj.n_state + nt +1) = H0;
    m_offset = m_offset + (nt +1)*obj.n_state;
    n_offset = n_offset + nt*(obj.n_state + obj.n_control);
    
end

    H = reshape(H_state, numel(H_state), 1);

    JH = JH_state'; % (XXX) to implement  
    
    G = [];
    JG = [];

end
    
   
    
    
    
    

