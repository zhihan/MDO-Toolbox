classdef StateConstraint
    % A state Jacobian 
    properties
        
        % An must implements nonlinear constraint g(xd, x) <= 0
        fun
        
        % An must implements the Jacobian of the nonlinear
        % constraint
        % [ J_d, J_s] = jacobian(xd, x)
        jacobian
        
        % To-do: Hessian of the Lagrangian 
        hessian
    end
    
end

