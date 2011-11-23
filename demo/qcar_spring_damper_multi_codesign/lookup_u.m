function u = lookup_u(tu, t)
i = find(tu.t <= t, 1, 'last');
u = tu.u(i);
end