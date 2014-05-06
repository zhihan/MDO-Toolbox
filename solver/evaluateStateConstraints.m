function y = evaluateStateConstraints(t, xd, xu, stateCon)

y1 = stateCon.fun(xd, xu);
m = length(y1);
y = zeros(m, length(t));
for i=1:length(t)
    y_i = stateCon.fun(xd, xu(:,i));
    y(:,i) = reshape(y_i, m,1);
end
end