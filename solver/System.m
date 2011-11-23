classdef System
% A System is essentially an ODE system
    methods
        
    end
    properties
        %%
        % Functions
        %
        % These properties are handles of the functions that defines the 
        % system
        
        deriv      % dx = f(t,x, x_d, parameters)
        
        % Until automatic differentiation is available, we require that
        % the user provide analytic Jacobians.
        
        % state Jacobian
        jacobian   % dx = d f/d x(t,x, x_d, parameters) 
        
        % design variable Jacobian
        jacobian_d   % dx_d = d f/d x_d(t,x,x_d, parameters)
        
        %%
        % Parameters are the variables excluding design variable and state
        % variables. These usually include constants used in the equations.
        parameters % constant parameters needed
        
        input      % u = g(t)
    end
end


