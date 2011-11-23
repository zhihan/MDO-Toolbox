function plot_result(t0, x0, t1, xu1, t2, xu2)

h2 = figure;
plot(t0,x0(:,3) + x0(:,1), 'r-.', 'LineWidth', 2);
hold on;
plot(t1, xu1(3,:) + xu1(1,:), 'k-', 'LineWidth',2);
plot(t2, xu2(3,:) + xu2(1,:), 'b--', 'LineWidth',2);

ylabel('z_s - z_0', 'FontSize', 14);
legend('Open-loop (no control)', 'DT', 'LQR');
xlabel('t', 'FontSize', 14);
set(get(h2, 'CurrentAxes'), 'FontSize', 14);
grid on;

end