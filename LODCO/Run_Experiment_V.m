clc; clear; close all;
V_list = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4];
result_cost = [];
result_battery = [];

disp('=== Starting V-parameter trade-off analysis ===');

for i = 1:length(V_list)
    v_val = V_list(i);
    fprintf('Testing V = %e ... \n', v_val);
    [c, b] = LODCO_Real_world_Modeling(v_val);
    
    result_cost(end+1) = c;
    result_battery(end+1) = b;
end

disp('Testing complete, generating plot...');
figure('Name', 'Energy-Delay Tradeoff');
plot(result_battery, result_cost, '-ro', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
xlabel('Average Battery Energy (J) [Energy Constraint]');
ylabel('Average Execution Cost (Delay) [Performance]');
title('Trade-off between Energy and Delay (Parameter V)');
for i = 1:length(V_list)
    text(result_battery(i), result_cost(i), ['  V=', num2str(V_list(i))], ...
        'VerticalAlignment', 'bottom', 'FontSize', 10);
end
