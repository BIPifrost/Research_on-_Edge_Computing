clc;
close all;

set(0, 'DefaultFigureVisible', 'off');

% Build a shared physical environment so both algorithms see exactly the same
% channels, task arrivals, and harvested energy in every slot/device.
seed_env = 20260307;
seed_origin_policy = 20260308;
seed_greedy_policy = 20260309;
env_cfg = struct( ...
    'T', 2000, ...
    'N', 10, ...
    'M', 8, ...
    'min_distance', 10, ...
    'max_distance', 50, ...
    'E_H_max', 48e-6, ...
    'rho', 0.6 ...
);
fixed_env = generate_fixed_env(seed_env, env_cfg);

% Run each script in an isolated function workspace to avoid variable
% cross-contamination between runs.
origin_raw = run_algorithm('eps_greedy_logco_origin.m', 'origin', fixed_env, seed_origin_policy);
greedy_raw = run_algorithm('greedy_logco_T2000.m', 'greedy', fixed_env, seed_greedy_policy);

% Align shape for fair comparison.
T_cmp = min(origin_raw.T, greedy_raw.T);
N_cmp = min(origin_raw.N, greedy_raw.N);
origin = trim_raw(origin_raw, T_cmp, N_cmp);
greedy = trim_raw(greedy_raw, T_cmp, N_cmp);

% Build multidimensional metrics.
metrics_origin = build_metrics(origin);
metrics_greedy = build_metrics(greedy);

summary = build_summary(metrics_origin, metrics_greedy, T_cmp, N_cmp, seed_env);

% Figure 1: Main dashboard.
fig_main = figure('Name', 'Comparison Main Dashboard', 'Position', [80 80 1450 860]);

subplot(3, 2, 1);
plot(1:T_cmp, metrics_origin.slot_cost_mean, 'b-', 'LineWidth', 1.0); hold on;
plot(1:T_cmp, metrics_greedy.slot_cost_mean, 'r-', 'LineWidth', 1.0);
plot(1:T_cmp, movmean(metrics_origin.slot_cost_mean, 50, 'omitnan'), 'b--', 'LineWidth', 1.6);
plot(1:T_cmp, movmean(metrics_greedy.slot_cost_mean, 50, 'omitnan'), 'r--', 'LineWidth', 1.6);
title('Slot Mean Cost + Moving Average');
xlabel('Time Slot');
ylabel('Cost (s)');
legend('eps_greedy', 'greedy', 'eps_greedy MA(50)', 'greedy MA(50)', 'Location', 'best');
grid on;

subplot(3, 2, 2);
plot(1:T_cmp, metrics_origin.cum_mean_cost, 'b-', 'LineWidth', 1.6); hold on;
plot(1:T_cmp, metrics_greedy.cum_mean_cost, 'r-', 'LineWidth', 1.6);
title('Cumulative Mean Cost');
xlabel('Time Slot');
ylabel('Cost (s)');
legend('eps_greedy', 'greedy', 'Location', 'best');
grid on;

subplot(3, 2, 3);
plot(1:T_cmp, metrics_origin.slot_energy_mean, 'b-', 'LineWidth', 1.0); hold on;
plot(1:T_cmp, metrics_greedy.slot_energy_mean, 'r-', 'LineWidth', 1.0);
plot(1:T_cmp, movmean(metrics_origin.slot_energy_mean, 50, 'omitnan'), 'b--', 'LineWidth', 1.6);
plot(1:T_cmp, movmean(metrics_greedy.slot_energy_mean, 50, 'omitnan'), 'r--', 'LineWidth', 1.6);
title('Slot Mean Energy + Moving Average');
xlabel('Time Slot');
ylabel('Energy (J)');
legend('eps_greedy', 'greedy', 'eps_greedy MA(50)', 'greedy MA(50)', 'Location', 'best');
grid on;

subplot(3, 2, 4);
plot(1:T_cmp, metrics_origin.cum_mean_energy, 'b-', 'LineWidth', 1.6); hold on;
plot(1:T_cmp, metrics_greedy.cum_mean_energy, 'r-', 'LineWidth', 1.6);
title('Cumulative Mean Energy');
xlabel('Time Slot');
ylabel('Energy (J)');
legend('eps_greedy', 'greedy', 'Location', 'best');
grid on;

subplot(3, 2, 5);
plot_battery_band(metrics_origin, [0.2 0.4 0.95]); hold on;
plot_battery_band(metrics_greedy, [0.95 0.3 0.3]);
title('Battery Mean and 10%-90% Band');
xlabel('Time Slot');
ylabel('Battery (J)');
legend('eps_greedy p10-p90', 'eps_greedy mean', 'greedy p10-p90', 'greedy mean', 'Location', 'best');
grid on;

subplot(3, 2, 6);
mode_bar = [metrics_origin.mode_ratio_overall; metrics_greedy.mode_ratio_overall];
h = bar(mode_bar, 'grouped');
h(1).FaceColor = [0.2 0.5 0.9];
h(2).FaceColor = [0.9 0.4 0.2];
h(3).FaceColor = [0.4 0.75 0.4];
set(gca, 'XTickLabel', {'eps_greedy', 'greedy'});
ylim([0 1]);
title('Overall Mode Ratio on Requested Tasks');
ylabel('Ratio');
legend({'Local', 'Remote', 'Drop'}, 'Location', 'best');
grid on;

sgtitle(sprintf('Main Dashboard | seed_env=%d | T=%d | N=%d', seed_env, T_cmp, N_cmp));
saveas(fig_main, 'comparison_dashboard_main.png');
savefig(fig_main, 'comparison_dashboard_main.fig');

% Figure 2: Detail dashboard.
fig_detail = figure('Name', 'Comparison Detail Dashboard', 'Position', [90 90 1450 860]);

subplot(3, 2, 1);
[xo, yo] = empirical_cdf(metrics_origin.request_cost_samples);
[xg, yg] = empirical_cdf(metrics_greedy.request_cost_samples);
plot(xo, yo, 'b-', 'LineWidth', 1.5); hold on;
plot(xg, yg, 'r-', 'LineWidth', 1.5);
title('Request-Level Cost CDF');
xlabel('Cost (s)');
ylabel('CDF');
legend('eps_greedy', 'greedy', 'Location', 'best');
grid on;

subplot(3, 2, 2);
[xo, yo] = empirical_cdf(metrics_origin.request_energy_samples);
[xg, yg] = empirical_cdf(metrics_greedy.request_energy_samples);
plot(xo, yo, 'b-', 'LineWidth', 1.5); hold on;
plot(xg, yg, 'r-', 'LineWidth', 1.5);
title('Request-Level Energy CDF');
xlabel('Energy (J)');
ylabel('CDF');
legend('eps_greedy', 'greedy', 'Location', 'best');
grid on;

subplot(3, 2, 3);
bar([metrics_origin.per_device_avg_cost, metrics_greedy.per_device_avg_cost], 'grouped');
title('Per-Device Average Cost');
xlabel('Device Index');
ylabel('Cost (s)');
legend('eps_greedy', 'greedy', 'Location', 'best');
grid on;

subplot(3, 2, 4);
bar([metrics_origin.per_device_avg_energy, metrics_greedy.per_device_avg_energy], 'grouped');
title('Per-Device Average Energy');
xlabel('Device Index');
ylabel('Energy (J)');
legend('eps_greedy', 'greedy', 'Location', 'best');
grid on;

subplot(3, 2, 5);
plot(1:T_cmp, metrics_origin.request_count_slot, 'b-', 'LineWidth', 1.2); hold on;
plot(1:T_cmp, metrics_greedy.request_count_slot, 'r-', 'LineWidth', 1.2);
plot(1:T_cmp, movmean(metrics_origin.request_count_slot, 50, 'omitnan'), 'b--', 'LineWidth', 1.6);
plot(1:T_cmp, movmean(metrics_greedy.request_count_slot, 50, 'omitnan'), 'r--', 'LineWidth', 1.6);
title('Request Load Per Slot');
xlabel('Time Slot');
ylabel('# Requested Devices');
legend('eps_greedy', 'greedy', 'eps_greedy MA(50)', 'greedy MA(50)', 'Location', 'best');
grid on;

subplot(3, 2, 6);
plot(1:T_cmp, metrics_origin.remote_ratio_cum, 'b-', 'LineWidth', 1.4); hold on;
plot(1:T_cmp, metrics_greedy.remote_ratio_cum, 'r-', 'LineWidth', 1.4);
plot(1:T_cmp, metrics_origin.drop_ratio_cum, 'b--', 'LineWidth', 1.4);
plot(1:T_cmp, metrics_greedy.drop_ratio_cum, 'r--', 'LineWidth', 1.4);
title('Cumulative Mode Ratio Evolution');
xlabel('Time Slot');
ylabel('Ratio');
legend('eps_greedy remote', 'greedy remote', 'eps_greedy drop', 'greedy drop', 'Location', 'best');
grid on;

sgtitle(sprintf('Detail Dashboard | seed_env=%d | T=%d | N=%d', seed_env, T_cmp, N_cmp));
saveas(fig_detail, 'comparison_dashboard_detail.png');
savefig(fig_detail, 'comparison_dashboard_detail.fig');

% Save result data.
result = struct();
result.timestamp = char(datetime('now'));
result.seed_env = seed_env;
result.seed_origin_policy = seed_origin_policy;
result.seed_greedy_policy = seed_greedy_policy;
result.T = T_cmp;
result.N = N_cmp;
result.origin = metrics_origin;
result.greedy = metrics_greedy;
result.summary = summary;
result.figures = {
    'comparison_dashboard_main.png', ...
    'comparison_dashboard_detail.png'
};
save('comparison_result.mat', 'result');

fprintf('Comparison completed.\n');
fprintf('Mean slot cost: origin=%.8f, greedy=%.8f\n', summary.mean_slot_cost_origin, summary.mean_slot_cost_greedy);
fprintf('Mean slot energy: origin=%.8e, greedy=%.8e\n', summary.mean_slot_energy_origin, summary.mean_slot_energy_greedy);
fprintf('Drop ratio: origin=%.6f, greedy=%.6f\n', summary.drop_ratio_origin, summary.drop_ratio_greedy);
fprintf('Seeds: env=%d, origin_policy=%d, greedy_policy=%d\n', seed_env, seed_origin_policy, seed_greedy_policy);

set(0, 'DefaultFigureVisible', 'on');

function raw = run_algorithm(script_name, tag, fixed_env, policy_seed)
    preserve_workspace = true;
    use_fixed_env = true;
    fixed_env_data = fixed_env;
    policy_stream = RandStream('mt19937ar', 'Seed', policy_seed);

    evalc(sprintf('run(''%s'');', script_name));

    required = {'T', 'N', 'B', 'final_chosen_cost', 'final_chosen_E', 'chosen_mode'};
    for k = 1:numel(required)
        if ~exist(required{k}, 'var')
            error('Script %s did not produce required variable %s.', script_name, required{k});
        end
    end
    raw = collect_raw(tag, T, N, B, final_chosen_cost, final_chosen_E, chosen_mode);
end

function fixed_env = generate_fixed_env(seed_env, cfg)
    rng(seed_env, 'twister');

    fixed_env = struct();
    fixed_env.seed_env = seed_env;
    fixed_env.T = cfg.T;
    fixed_env.N = cfg.N;
    fixed_env.M = cfg.M;
    fixed_env.distances = unifrnd(cfg.min_distance, cfg.max_distance, cfg.N, cfg.M, cfg.T);
    fixed_env.gamma = exprnd(1, cfg.N, cfg.M, cfg.T);
    fixed_env.E_H = unifrnd(0, cfg.E_H_max, cfg.T, cfg.N);
    fixed_env.zeta = binornd(1, cfg.rho, cfg.T, cfg.N);
end

function raw = collect_raw(tag, T, N, B, final_cost, final_energy, chosen_mode)
    raw = struct();
    raw.tag = tag;
    raw.T = T;
    raw.N = N;
    raw.B = B;
    raw.final_cost = final_cost;
    raw.final_energy = final_energy;
    raw.chosen_mode = chosen_mode;
end

function out = trim_raw(raw, T_cmp, N_cmp)
    out = struct();
    out.tag = raw.tag;
    out.T = T_cmp;
    out.N = N_cmp;
    out.B = raw.B(1:T_cmp, 1:N_cmp);
    out.final_cost = raw.final_cost(1:T_cmp, 1:N_cmp);
    out.final_energy = raw.final_energy(1:T_cmp, 1:N_cmp);
    out.chosen_mode = raw.chosen_mode(1:T_cmp, 1:N_cmp);
end

function m = build_metrics(data)
    T = data.T;
    N = data.N;
    requested = (data.chosen_mode ~= 4);

    m = struct();
    m.T = T;
    m.N = N;
    m.request_count_slot = sum(requested, 2);
    m.slot_cost_mean = nan(T, 1);
    m.slot_energy_mean = nan(T, 1);
    m.slot_cost_p50 = nan(T, 1);
    m.slot_cost_p90 = nan(T, 1);
    m.slot_energy_p50 = nan(T, 1);
    m.slot_energy_p90 = nan(T, 1);

    for t = 1:T
        idx = requested(t, :);
        if any(idx)
            costs = data.final_cost(t, idx);
            energies = data.final_energy(t, idx);
            m.slot_cost_mean(t) = mean(costs, 'omitnan');
            m.slot_energy_mean(t) = mean(energies, 'omitnan');
            m.slot_cost_p50(t) = simple_percentile(costs, 50);
            m.slot_cost_p90(t) = simple_percentile(costs, 90);
            m.slot_energy_p50(t) = simple_percentile(energies, 50);
            m.slot_energy_p90(t) = simple_percentile(energies, 90);
        end
    end

    req_count_cum = cumsum(m.request_count_slot);
    req_cost_sum_cum = cumsum(replace_nan_with_zero(m.slot_cost_mean .* m.request_count_slot));
    req_energy_sum_cum = cumsum(replace_nan_with_zero(m.slot_energy_mean .* m.request_count_slot));
    m.cum_mean_cost = req_cost_sum_cum ./ max(req_count_cum, 1);
    m.cum_mean_energy = req_energy_sum_cum ./ max(req_count_cum, 1);

    m.battery_mean = mean(data.B, 2, 'omitnan');
    m.battery_p10 = nan(T, 1);
    m.battery_p90 = nan(T, 1);
    for t = 1:T
        m.battery_p10(t) = simple_percentile(data.B(t, :), 10);
        m.battery_p90(t) = simple_percentile(data.B(t, :), 90);
    end

    req_modes = data.chosen_mode(requested);
    req_count_total = numel(req_modes);
    if req_count_total == 0
        m.mode_ratio_overall = [0, 0, 0];
    else
        m.mode_ratio_overall = [
            sum(req_modes == 1) / req_count_total, ...
            sum(req_modes == 2) / req_count_total, ...
            sum(req_modes == 3) / req_count_total
        ];
    end

    local_slot = sum((data.chosen_mode == 1), 2);
    remote_slot = sum((data.chosen_mode == 2), 2);
    drop_slot = sum((data.chosen_mode == 3), 2);
    req_cum = cumsum(m.request_count_slot);
    m.local_ratio_cum = cumsum(local_slot) ./ max(req_cum, 1);
    m.remote_ratio_cum = cumsum(remote_slot) ./ max(req_cum, 1);
    m.drop_ratio_cum = cumsum(drop_slot) ./ max(req_cum, 1);

    m.per_device_avg_cost = nan(N, 1);
    m.per_device_avg_energy = nan(N, 1);
    for i = 1:N
        idx = requested(:, i);
        if any(idx)
            m.per_device_avg_cost(i) = mean(data.final_cost(idx, i), 'omitnan');
            m.per_device_avg_energy(i) = mean(data.final_energy(idx, i), 'omitnan');
        else
            m.per_device_avg_cost(i) = 0;
            m.per_device_avg_energy(i) = 0;
        end
    end

    m.request_cost_samples = data.final_cost(requested);
    m.request_energy_samples = data.final_energy(requested);
end

function s = build_summary(o, g, T_cmp, N_cmp, seed)
    s = struct();
    s.seed = seed;
    s.T_compared = T_cmp;
    s.N_compared = N_cmp;
    % Use request-weighted averages to avoid bias from unequal per-slot load.
    s.mean_slot_cost_origin = o.cum_mean_cost(end);
    s.mean_slot_cost_greedy = g.cum_mean_cost(end);
    s.mean_slot_energy_origin = o.cum_mean_energy(end);
    s.mean_slot_energy_greedy = g.cum_mean_energy(end);
    s.mean_battery_origin = mean(o.battery_mean, 'omitnan');
    s.mean_battery_greedy = mean(g.battery_mean, 'omitnan');
    s.local_ratio_origin = o.mode_ratio_overall(1);
    s.local_ratio_greedy = g.mode_ratio_overall(1);
    s.remote_ratio_origin = o.mode_ratio_overall(2);
    s.remote_ratio_greedy = g.mode_ratio_overall(2);
    s.drop_ratio_origin = o.mode_ratio_overall(3);
    s.drop_ratio_greedy = g.mode_ratio_overall(3);
end

function plot_battery_band(m, color_base)
    x = (1:m.T)';
    c_fill = min(color_base + 0.35, 1);
    fill([x; flipud(x)], [m.battery_p10; flipud(m.battery_p90)], c_fill, ...
        'FaceAlpha', 0.30, 'EdgeColor', 'none'); hold on;
    plot(x, m.battery_mean, '-', 'Color', color_base, 'LineWidth', 1.6);
end

function y = replace_nan_with_zero(x)
    y = x;
    y(isnan(y)) = 0;
end

function p = simple_percentile(vec, pct)
    v = vec(~isnan(vec));
    if isempty(v)
        p = nan;
        return;
    end
    v = sort(v(:));
    n = numel(v);
    idx = max(1, min(n, round((pct / 100) * n)));
    p = v(idx);
end

function [x, y] = empirical_cdf(samples)
    v = samples(~isnan(samples));
    if isempty(v)
        x = [0; 1];
        y = [0; 1];
        return;
    end
    x = sort(v(:));
    y = (1:numel(x))' / numel(x);
end
