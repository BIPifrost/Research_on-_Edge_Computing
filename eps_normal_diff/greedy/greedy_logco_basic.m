function result = greedy_logco_basic(cfg)
% GREEDY_LOGCO_BASIC
% Configurable greedy LODCO simulation that preserves the math/decision
% logic from greedy_logco.m while exposing parameters for sweeps.

if nargin < 1
    cfg = struct();
end

def = struct();
def.seed = 20260307;
def.k = 1e-28;
def.tau = 0.002;
def.phi = 0.002;
def.omega = 1e6;
def.sigma = 1e-13;
def.p_tx_max = 1;
def.f_max = 1.5e9;
def.E_max = 0.002;
def.L = 1000;
def.X = 737.5;
def.E_H_max = 48e-6;
def.g0 = 1e-4;
def.d0 = 1;
def.N = 10;
def.M = 8;
def.T = 1000;
def.tau_d = 0.002;
def.d = 50;
def.E_min = 0.02e-3;
def.V = 1e-5;
def.rho = 0.6;
def.max_connects = 4;
def.min_distance = 10;
def.max_distance = 50;
def.show_progress = false;
def.make_plot = false;
def.plot_title = 'greedy mode selection';
def.save_plot = false;
def.output_prefix = 'greedy_basic';

p = apply_overrides(def, cfg);
validate_params(p);
rng(p.seed, 'twister');

W = p.L * p.X;
E_max_hat = min(max(p.k * W * (p.f_max)^2, p.p_tx_max * p.tau), p.E_max);
theta = E_max_hat + p.V * p.phi / p.E_min;

T = p.T;
N = p.N;
M = p.M;

B = zeros(T + 1, N);
B_hat = zeros(T, N);
e = zeros(T, N);
chosen_mode = zeros(T, N);
chosen_server = zeros(T, N);
f = zeros(T, N);
ptx = zeros(T, N);

mobile_exe_cost = zeros(T, N);
server_exe_cost = zeros(T, N);
final_chosen_cost = zeros(T, N);
mobile_exe_E = zeros(T, N);
server_exe_E = zeros(T, N);
final_chosen_E = zeros(T, N);

opt = optimset('Display', 'none');

for t = 1:T
    if p.show_progress && mod(t, 50) == 0
        fprintf('slot %d / %d\n', t, T);
    end

    device_server_pairs = [];
    remained_connects = p.max_connects * ones(M, 1);
    J_m = zeros(N, 1);
    J_s = zeros(N, M);
    p_mat = zeros(N, M);
    server_cost_mat = zeros(N, M);
    server_E_mat = zeros(N, M);

    B_hat(t, :) = B(t, :) - theta;
    distances = unifrnd(p.min_distance, p.max_distance, N, M);
    gamma = exprnd(1, N, M);
    h_mat = p.g0 .* gamma .* power(p.d0 ./ distances, 4);

    for i = 1:N
        E_H_t = unifrnd(0, p.E_H_max);
        if B_hat(t, i) <= 0
            e(t, i) = E_H_t;
        end

        zeta = binornd(1, p.rho);
        if zeta == 0
            chosen_mode(t, i) = 4;
            continue;
        end

        f_L = max(sqrt(p.E_min / (p.k * W)), W / p.tau_d);
        f_U = min(sqrt(p.E_max / (p.k * W)), p.f_max);
        if f_L <= f_U
            if B_hat(t, i) < 0
                f_0 = (p.V / (-2 * B_hat(t, i) * p.k))^(1/3);
            else
                f_0 = -(p.V / (2 * B_hat(t, i) * p.k))^(1/3);
            end

            if (f_0 > f_U && B_hat(t, i) < 0) || (B_hat(t, i) >= 0)
                f(t, i) = f_U;
            elseif f_0 >= f_L && f_0 <= f_U && B_hat(t, i) < 0
                f(t, i) = f_0;
            else
                f(t, i) = f_L;
            end

            mobile_exe_cost(t, i) = W / f(t, i);
            mobile_exe_E(t, i) = p.k * W * (f(t, i)^2);
            J_m(i) = -B_hat(t, i) * p.k * W * (f(t, i)^2) + p.V * W / f(t, i);
        else
            J_m(i) = inf;
        end

        for j = 1:M
            h = h_mat(i, j);
            E_tmp = p.sigma * p.L * log(2) / (p.omega * h);
            p_L_taud = (power(2, p.L / (p.omega * p.tau_d)) - 1) * p.sigma / h;

            if E_tmp >= p.E_min
                p_L = p_L_taud;
            else
                y = @(x) x * p.L - p.omega * log2(1 + h * x / p.sigma) * p.E_min;
                p_E_min = fsolve(y, 0.2, opt);
                p_L = max(p_L_taud, p_E_min);
            end

            if E_tmp >= p.E_max
                p_U = 0;
            else
                y = @(x) x * p.L - p.omega * log2(1 + h * x / p.sigma) * p.E_max;
                p_E_max = fsolve(y, 100, opt);
                p_U = min(p.p_tx_max, p_E_max);
            end

            if p_L <= p_U
                virtual_battery = B_hat(t, i);
                y = @(x) virtual_battery * log2(1 + h * x / p.sigma) + ...
                    h * (p.V - virtual_battery * x) / log(2) / (p.sigma + h * x);
                p_0 = fsolve(y, 0.5, opt);

                if (p_U < p_0 && B_hat(t, i) < 0) || B_hat(t, i) >= 0
                    p_mat(i, j) = p_U;
                elseif p_0 < p_L && B_hat(t, i) < 0
                    p_mat(i, j) = p_L;
                else
                    p_mat(i, j) = p_0;
                end

                server_cost_mat(i, j) = p.L / (p.omega * log2(1 + h * p_mat(i, j) / p.sigma));
                server_E_mat(i, j) = p_mat(i, j) * server_cost_mat(i, j);
                J_s(i, j) = (-B_hat(t, i) * p_mat(i, j) + p.V) * server_cost_mat(i, j);
            else
                J_s(i, j) = inf;
            end
        end

        J_d = p.V * p.phi;
        [~, mode] = min([J_m(i), J_s(i, :), J_d]);
        if mode == 1
            chosen_mode(t, i) = 1;
            final_chosen_cost(t, i) = mobile_exe_cost(t, i);
            final_chosen_E(t, i) = mobile_exe_E(t, i);
        elseif mode == (M + 2)
            chosen_mode(t, i) = 3;
            final_chosen_cost(t, i) = p.phi;
            final_chosen_E(t, i) = 0;
        else
            chosen_mode(t, i) = 2;
            device_server_pairs = [device_server_pairs; [i, mode - 1, J_s(i, mode - 1)]]; %#ok<AGROW>
        end
    end

    while ~isempty(device_server_pairs)
        for j = 1:M
            device_j_pairs = device_server_pairs(device_server_pairs(:, 2) == j, :);
            is = device_j_pairs(:, 1);
            if isempty(is)
                continue;
            end

            if remained_connects(j) >= numel(is)
                chosen_mode(t, is) = 2;
                ptx(t, is) = reshape(p_mat(is, j), 1, []);
                server_exe_cost(t, is) = reshape(server_cost_mat(is, j), 1, []);
                server_exe_E(t, is) = reshape(server_E_mat(is, j), 1, []);
                chosen_server(t, is) = j;
                final_chosen_cost(t, is) = server_exe_cost(t, is);
                final_chosen_E(t, is) = server_exe_E(t, is);
                remained_connects(j) = remained_connects(j) - numel(is);
                device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
            elseif remained_connects(j) == 0
                device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
                J_s(is, j) = inf;
                J_d = p.V * p.phi;
                for idx = 1:numel(is)
                    ii = is(idx);
                    [~, mode2] = min([J_m(ii), J_s(ii, :), J_d]);
                    if mode2 == 1
                        chosen_mode(t, ii) = 1;
                        final_chosen_cost(t, ii) = mobile_exe_cost(t, ii);
                        final_chosen_E(t, ii) = mobile_exe_E(t, ii);
                    elseif mode2 == (M + 2)
                        chosen_mode(t, ii) = 3;
                        final_chosen_cost(t, ii) = p.phi;
                        final_chosen_E(t, ii) = 0;
                    else
                        chosen_mode(t, ii) = 2;
                        device_server_pairs = [device_server_pairs; [ii, mode2 - 1, J_s(ii, mode2 - 1)]]; %#ok<AGROW>
                    end
                end
            else
                [~, idxs] = sort(device_j_pairs(:, 3));
                keep_k = remained_connects(j);
                for idx = 1:keep_k
                    ii = device_j_pairs(idxs(idx), 1);
                    chosen_mode(t, ii) = 2;
                    ptx(t, ii) = p_mat(ii, j);
                    server_exe_cost(t, ii) = server_cost_mat(ii, j);
                    server_exe_E(t, ii) = server_E_mat(ii, j);
                    chosen_server(t, ii) = j;
                    final_chosen_cost(t, ii) = server_exe_cost(t, ii);
                    final_chosen_E(t, ii) = server_exe_E(t, ii);
                    device_server_pairs(device_server_pairs(:, 1) == ii, :) = [];
                end

                remained_connects(j) = 0;
                residual_is = device_server_pairs(device_server_pairs(:, 2) == j, 1);
                J_s(residual_is, j) = inf;
                J_d = p.V * p.phi;
                for idx = 1:numel(residual_is)
                    ii = residual_is(idx);
                    % Remove the stale queue entry on server j before re-queuing ii.
                    device_server_pairs(device_server_pairs(:, 1) == ii & device_server_pairs(:, 2) == j, :) = [];
                    [~, mode3] = min([J_m(ii), J_s(ii, :), J_d]);
                    if mode3 == 1
                        chosen_mode(t, ii) = 1;
                        final_chosen_cost(t, ii) = mobile_exe_cost(t, ii);
                        final_chosen_E(t, ii) = mobile_exe_E(t, ii);
                    elseif mode3 == (M + 2)
                        chosen_mode(t, ii) = 3;
                        final_chosen_cost(t, ii) = p.phi;
                        final_chosen_E(t, ii) = 0;
                    else
                        chosen_mode(t, ii) = 2;
                        device_server_pairs = [device_server_pairs; [ii, mode3 - 1, J_s(ii, mode3 - 1)]]; %#ok<AGROW>
                    end
                end
            end
        end
    end

    B(t + 1, :) = B(t, :) - final_chosen_E(t, :) + e(t, :);
end

requested = (chosen_mode ~= 4);
req_per_slot = sum(requested, 2);
local_slot = sum(chosen_mode == 1, 2);
remote_slot = sum(chosen_mode == 2, 2);
drop_slot = sum(chosen_mode == 3, 2);

safe_req_slot = max(req_per_slot, 1);
mode_ratio_slot = [local_slot ./ safe_req_slot, remote_slot ./ safe_req_slot, drop_slot ./ safe_req_slot];

req_cum = cumsum(req_per_slot);
safe_req_cum = max(req_cum, 1);
mode_ratio_cum = [cumsum(local_slot) ./ safe_req_cum, cumsum(remote_slot) ./ safe_req_cum, cumsum(drop_slot) ./ safe_req_cum];

req_modes = chosen_mode(requested);
if isempty(req_modes)
    mode_ratio_overall = [0, 0, 0];
else
    mode_ratio_overall = [sum(req_modes == 1), sum(req_modes == 2), sum(req_modes == 3)] ./ numel(req_modes);
end

result = struct();
result.params = p;
result.theta = theta;
result.B = B(1:T, :);
result.chosen_mode = chosen_mode;
result.chosen_server = chosen_server;
result.final_chosen_cost = final_chosen_cost;
result.final_chosen_E = final_chosen_E;
result.mode_ratio_slot = mode_ratio_slot;
result.mode_ratio_cum = mode_ratio_cum;
result.mode_ratio_overall = mode_ratio_overall;
result.req_per_slot = req_per_slot;
result.mean_cost = mean(final_chosen_cost(requested), 'omitnan');
result.mean_energy = mean(final_chosen_E(requested), 'omitnan');

if p.make_plot
    plot_mode_profile(result, p.plot_title);
    if p.save_plot
        saveas(gcf, [p.output_prefix, '_mode_profile.png']);
        savefig(gcf, [p.output_prefix, '_mode_profile.fig']);
    end
end
end

function out = apply_overrides(base, cfg)
out = base;
if isempty(cfg)
    return;
end

fields = fieldnames(cfg);
for k = 1:numel(fields)
    out.(fields{k}) = cfg.(fields{k});
end
end

function validate_params(p)
if p.min_distance <= 0 || p.max_distance <= p.min_distance
    error('Invalid distance bounds: require 0 < min_distance < max_distance.');
end
if p.N < 1 || p.M < 1 || p.T < 1
    error('N, M, T must be positive integers.');
end
if p.max_connects < 1
    error('max_connects must be >= 1.');
end
if p.rho < 0 || p.rho > 1
    error('rho must be in [0, 1].');
end
if p.tau_d <= 0
    error('tau_d must be > 0.');
end
end

function plot_mode_profile(result, fig_title)
T = size(result.mode_ratio_slot, 1);
figure('Name', fig_title, 'Position', [120 120 1300 500]);

subplot(1, 2, 1);
plot(1:T, result.mode_ratio_slot(:, 1), 'LineWidth', 1.0); hold on;
plot(1:T, result.mode_ratio_slot(:, 2), 'LineWidth', 1.0);
plot(1:T, result.mode_ratio_slot(:, 3), 'LineWidth', 1.0);
plot(1:T, movmean(result.mode_ratio_slot(:, 1), 50, 'omitnan'), '--', 'LineWidth', 1.4);
plot(1:T, movmean(result.mode_ratio_slot(:, 2), 50, 'omitnan'), '--', 'LineWidth', 1.4);
plot(1:T, movmean(result.mode_ratio_slot(:, 3), 50, 'omitnan'), '--', 'LineWidth', 1.4);
legend('local', 'remote', 'drop', 'local MA50', 'remote MA50', 'drop MA50', 'Location', 'best');
title('Per-Slot Mode Ratio');
xlabel('Time Slot');
ylabel('Ratio');
grid on;

subplot(1, 2, 2);
bar(result.mode_ratio_overall, 'FaceColor', [0.2 0.5 0.85]);
xticks([1, 2, 3]);
xticklabels({'Local', 'Remote', 'Drop'});
ylim([0, 1]);
title('Overall Mode Ratio');
ylabel('Ratio');
grid on;

sgtitle(fig_title);
end
