[system, objective, init, designC, stateC] = simple_sys();

s = SequentialSolverWithInputs();
s.system = system;
s.objective = objective;
s.initCondition = init;
s.designConstraint = designC;
s.stateConstraint = stateC;

t = 0:0.1:1;
s.n_control = 1;
s.t ={t};
s.input = {@(t)[t;0]};
s.pinit = 0.6;
s.x0 = [];
s.lb = [];
s.ub = [];
s.weight = {1};

s.options = OptionFactory.getInstance.makeOption('');
s.options(2) = OptionFactory.getInstance.makeOption('');
out = s.f_solve();