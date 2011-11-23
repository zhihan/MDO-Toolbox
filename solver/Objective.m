classdef Objective
% An objective function formulated for state trajectories

% The trajectory must be in the format of [x(1),x(2),x(3)...] where
% each x could in turns be vectors. For systems with inputs the state
% vector has the form [x;u];
    methods
    end
    
    properties
        fun      % xout = fun(x)
        gradient % 
        hessian
        
        fun_d
        gradient_d
        
        QR % may be empty
        
    end
end
