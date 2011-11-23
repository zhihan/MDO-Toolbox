classdef Constraint
    % Design constraints
    
    properties
        ub
        lb
        A
        b
        fun
        jacobian
        hessian
    end
end