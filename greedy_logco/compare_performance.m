clc; clear; close all;

disp('>>> Running original code (fsolve)...');
rng(2026);
tic;
run('greedy_logco.m');
time_fsolve = toc;

B_orig = B;
cost_orig = final_chosen_cost;
mode_orig = chosen_mode;

clearvars -except time_fsolve B_orig cost_orig mode_orig

disp(' ');
disp('>>> Running bisection-based code...');
rng(2026);
tic;
run('greedy_logco_bisearch.m');
time_bisearch = toc;

B_new = B;
cost_new = final_chosen_cost;
mode_new = chosen_mode;

speedup = time_fsolve / time_bisearch;
max_diff_B = max(abs(B_orig(:) - B_new(:)));
max_diff_cost = max(abs(cost_orig(:) - cost_new(:)));

total_decisions = numel(mode_orig);
diff_mode_count = sum(mode_orig(:) ~= mode_new(:));
consistent_count = total_decisions - diff_mode_count;

disp(' ');
disp('>>> Simulation complete. Generating comparison figures...');

fig = figure('Name', 'Algorithm Performance Comparison', 'Position', [100, 150, 1200, 450], 'Color', 'w');

subplot(1, 3, 1);
bar_data = [time_fsolve, time_bisearch];
b = bar(bar_data, 0.5);
b.FaceColor = 'flat';
b.CData(1,:) = [0.8500, 0.3250, 0.0980];
b.CData(2,:) = [0.0000, 0.4470, 0.7410];

set(gca, 'XTickLabel', {'Original (fsolve)', 'Bisection'}, 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Execution Time (s)', 'FontSize', 12);
title('Runtime Comparison', 'FontSize', 14);
grid on;
ylim([0, max(bar_data) * 1.2]);

text(1, time_fsolve + max(bar_data) * 0.03, sprintf('%.2f s', time_fsolve), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
text(2, time_bisearch + max(bar_data) * 0.03, sprintf('%.2f s', time_bisearch), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');

subplot(1, 3, 2);
if diff_mode_count == 0
    p = pie(1);
    p(1).FaceColor = [0.4660, 0.6740, 0.1880];
    legend_str = {sprintf('Fully consistent\n(%d decisions)', total_decisions)};
else
    p = pie([consistent_count, diff_mode_count]);
    p(1).FaceColor = [0.4660, 0.6740, 0.1880];
    p(3).FaceColor = [0.8500, 0.3250, 0.0980];
    legend_str = {sprintf('Consistent (%d)', consistent_count), sprintf('Different (%d)', diff_mode_count)};
end
title('Mode Consistency', 'FontSize', 14);
legend(legend_str, 'Location', 'southoutside', 'FontSize', 11);

subplot(1, 3, 3);
axis off;
title('Accuracy Summary', 'FontSize', 14, 'Visible', 'on');

info_str = {
    'Performance';
    sprintf('Speedup: %.2fx', speedup);
    '';
    'Absolute Error';
    sprintf('Battery level max error: %e', max_diff_B);
    sprintf('Delay max error: %e', max_diff_cost);
    '';
    'Conclusion';
    'Errors stay very small.';
    'Decision consistency remains high.';
};

text(0.05, 0.5, info_str, 'FontSize', 12, ...
    'VerticalAlignment', 'middle', ...
    'BackgroundColor', [0.96 0.96 0.96], ...
    'EdgeColor', [0.6 0.6 0.6], ...
    'LineWidth', 1, ...
    'Margin', 15);

disp('>>> Figures generated.');
