classdef DesignConstraint
    % Design constraints
    % Design constraints are constraints on the design variables.
    
    properties
        % nonlinear functions
        fun
        
        % The Jacobian of the function fun
        jacobian
        
        % The Hessian of the Lagrangian function
        hessian  
        
    end
end