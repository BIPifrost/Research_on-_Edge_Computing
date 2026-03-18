clc; clear; close all;

data1 = load('data_original.mat');
data2 = load('data_optimized.mat');

T = length(data1.time_history); 

figure('Name', 'Efficiency Comparison');
time_axis = 1:T;

cumsum_old = cumsum(data1.time_history);
cumsum_new = cumsum(data2.time_history);

plot(time_axis, cumsum_old, 'r-', 'LineWidth', 1.5); hold on;
plot(time_axis, cumsum_new, 'b-', 'LineWidth', 1.5);

legend('Original (fsolve)', 'Optimized (Bisection)', 'Location', 'NorthWest');
title('Cumulative Computation Time Comparison');
xlabel('Time Slot');
ylabel('Total Time (seconds)');
grid on;

speedup = cumsum_old(end) / cumsum_new(end);
text(T*0.1, cumsum_old(end)*0.8, ...
    ['Speedup: ', num2str(speedup, '%.1f'), 'x'], ...
    'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', 'w');

figure('Name', 'Battery Dynamics Comparison');
plot(time_axis, data1.B(1:T), 'r-', 'LineWidth', 0.5); hold on;
plot(time_axis, data2.B(1:T), 'b-', 'LineWidth', 0.5);

legend('Original (Conservative)', 'Optimized (Active/Aggressive)');
title('Battery Energy Level Evolution');
xlabel('Time Slot');
ylabel('Battery Level (J)');
grid on;
xlim([0, 2000]);

figure('Name', 'Cost Convergence Comparison');
plot(time_axis, data1.average_cost, 'r--', 'LineWidth', 1.5); hold on;
plot(time_axis, data2.average_cost, 'b-', 'LineWidth', 1.5);

legend('Original Cost', 'Optimized Cost');
title('Average Execution Cost Convergence');
xlabel('Time Slot');
ylabel('Average Cost (Delay + Weight*Energy)');
grid on;

disp('Comparison plots generated successfully.');
