function checkSolver(obj)
% Check solver option 

% A solver must have setup its options
if isempty(obj.options)
    error('Solver:NoOption', 'Solver has not specified option');
end

end