function out = calc_zdot(x, z, v)

n = length(x);
t = zeros(length(x),1);
zdot = zeros(length(x)-1, 1);
for i=2:n
    dist = sqrt((z(i) - z(i-1))^2 + (x(i)-x(i-1))^2);
    t(i) = t(i-1) + dist/v;
    zdot(i-1) = (z(i) - z(i-1))/(t(i) - t(i-1));
end

out = struct('t', t, 'u', zdot);

end