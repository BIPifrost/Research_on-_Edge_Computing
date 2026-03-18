clc; clear; close all;
opt = optimset('Display', 'none');

%% Basic parameter settings (aligned with eps_greedy_logco_bi)
params.k = 1e-28;
params.tau = 0.002;
params.phi = 0.002;
params.omega = 1e6;
params.sigma = 1e-13;
params.p_tx_max = 1;
params.f_max = 1.5e9;
params.E_max = 0.002;
params.L = 1000;
params.X = 737.5;
params.W = params.L * params.X;
params.E_H_max = 48e-6;
params.p_H = params.E_H_max / (2 * params.tau);
params.g0 = 10^-4;
params.d0 = 1;

%% Control parameters
params.N = 10;                 % number of mobile devices
params.T = 1000;               % time slots (reduce for multi-group simulation speed)
params.tau_d = 0.002;
params.d = 50;
params.E_min = 0.02e-3;
params.V = 1e-5;
params.rho = 0.6;
params.max_connects = 4;
params.min_distance = 10;
params.max_distance = 50;
params.eps = 0.8;

M_list = [8, 6, 4, 2];         % multi-group MEC server counts (decreasing)

num_cases = numel(M_list);
avg_cost_all = zeros(params.T, num_cases);
final_cost_all = zeros(num_cases, 1);
final_ratio_all = zeros(num_cases, 3);  % [mobile, server, drop]

for c = 1:num_cases
    M = M_list(c);
    rng(20260307, 'twister');  % fixed seed for fair comparison
    result = simulate_one_case(params, M, opt);
    avg_cost_all(:, c) = result.avg_network_cost;
    final_cost_all(c) = result.final_avg_cost;
    final_ratio_all(c, :) = result.final_ratio;
    fprintf('M=%d | Final Avg Cost=%.6f | Ratio[m,s,d]=[%.3f, %.3f, %.3f]\n', ...
        M, result.final_avg_cost, result.final_ratio(1), result.final_ratio(2), result.final_ratio(3));
end

%% Figure 1: average execution cost vs time slot
figure('Name', 'Average Cost Comparison under Different MEC Counts');
plot(1:params.T, avg_cost_all, 'LineWidth', 1.2);
grid on;
xlabel('time slot');
ylabel('network average execution cost (s)');
title('Average Execution Cost under Different MEC Server Counts');
legend(compose('M=%d', M_list), 'Location', 'best');

%% Figure 2: final mode ratio comparison
figure('Name', 'Mode Ratio Comparison under Different MEC Counts');
bar(categorical(compose('M=%d', M_list)), final_ratio_all, 'grouped');
grid on;
xlabel('MEC server count');
ylabel('final average ratio');
title('Final Mode Ratio Comparison');
legend({'mobile execution', 'MEC server execution', 'drop'}, 'Location', 'best');

%% Figure 3: final average cost vs MEC count
figure('Name', 'Final Average Cost vs MEC Count');
plot(M_list, final_cost_all, '-o', 'LineWidth', 1.4, 'MarkerSize', 6);
grid on;
xlabel('MEC server count');
ylabel('final network average execution cost (s)');
title('Final Average Cost vs MEC Server Count');
set(gca, 'XDir', 'reverse');

%% --------------------------- Local functions ----------------------------
function result = simulate_one_case(p, M, opt)
N = p.N; T = p.T;

E_max_hat = min(max(p.k * p.W * (p.f_max)^2, p.p_tx_max * p.tau), p.E_max);
theta = E_max_hat + p.V * p.phi / p.E_min;

B = zeros(T + 1, N);
B_hat = zeros(T, N);
e = zeros(T, N);
chosen_mode = zeros(T, N);      % 1:local, 2:server, 3:drop, 4:no request
f = zeros(T, N);
ptx = zeros(T, N);
final_cost = zeros(T, N);
final_E = zeros(T, N);

mobile_exe_cost = zeros(T, N);
mobile_exe_E = zeros(T, N);
server_cost_mat = zeros(N, M);
server_E_mat = zeros(N, M);

for t = 1:T
    B_hat(t, :) = B(t, :) - theta;
    device_server_pairs = [];
    remained_connects = p.max_connects * ones(M, 1);
    J_m = zeros(N, 1);
    J_s = inf(N, M);
    p_mat = zeros(N, M);
    server_cost_mat(:, :) = 0;
    server_E_mat(:, :) = 0;
    J_d = p.V * p.phi;

    distances = unifrnd(p.min_distance, p.max_distance, N, M);
    gamma = exprnd(1, N, M);
    h_mat = p.g0 * gamma .* (p.d0 ./ distances).^4;

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

        % Solve P_ME
        f_L = max(sqrt(p.E_min / (p.k * p.W)), p.W / p.tau_d);
        f_U = min(sqrt(p.E_max / (p.k * p.W)), p.f_max);
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
            mobile_exe_cost(t, i) = p.W / f(t, i);
            mobile_exe_E(t, i) = p.k * p.W * f(t, i)^2;
            J_m(i) = -B_hat(t, i) * p.k * p.W * f(t, i)^2 + p.V * p.W / f(t, i);
        else
            J_m(i) = inf;
        end

        % Solve P_SE for all MEC servers
        for j = 1:M
            h = h_mat(i, j);
            E_tmp = p.sigma * p.L * log(2) / (p.omega * h);
            p_L_taud = (2^(p.L / (p.omega * p.tau_d)) - 1) * p.sigma / h;

            if E_tmp >= p.E_min
                p_L = p_L_taud;
            else
                y_E_min = @(x) x * p.L - p.omega * log2(1 + h * x / p.sigma) * p.E_min;
                a_min = max(0, (p.omega * p.E_min) / (p.L * log(2)) - p.sigma / h);
                p_E_min = my_bisection(y_E_min, a_min, 100, 1e-7);
                if isinf(p_E_min) || isnan(p_E_min)
                    [p_E_min, ~, flag1] = fsolve(y_E_min, 0.2, opt);
                    if flag1 <= 0 || p_E_min < 0 || ~isreal(p_E_min)
                        p_E_min = inf;
                    end
                end
                p_L = max(p_L_taud, p_E_min);
            end

            if E_tmp >= p.E_max
                p_U = 0;
            else
                y_E_max = @(x) x * p.L - p.omega * log2(1 + h * x / p.sigma) * p.E_max;
                a_max = max(0, (p.omega * p.E_max) / (p.L * log(2)) - p.sigma / h);
                if y_E_max(p.p_tx_max) <= 0
                    p_E_max = inf;
                else
                    p_E_max = my_bisection(y_E_max, a_max, p.p_tx_max, 1e-7);
                end
                if isinf(p_E_max) || isnan(p_E_max)
                    [p_E_max, ~, flag2] = fsolve(y_E_max, 100, opt);
                    if flag2 <= 0 || p_E_max < 0 || ~isreal(p_E_max)
                        p_E_max = inf;
                    end
                end
                p_U = min(p.p_tx_max, p_E_max);
            end

            p_0 = inf;
            if p_L <= p_U
                vb = B_hat(t, i);
                y_p0 = @(x) vb * log2(1 + h * x / p.sigma) + ...
                    h * (p.V - vb * x) / log(2) / (p.sigma + h * x);
                p_0 = my_bisection(y_p0, 0, 100, 1e-7);
                if isinf(p_0) || isnan(p_0)
                    [p_0, ~, flag3] = fsolve(y_p0, 0.5, opt);
                    if flag3 <= 0 || p_0 < 0 || ~isreal(p_0)
                        p_0 = inf;
                    end
                end

                if (p_U < p_0 && B_hat(t, i) < 0) || B_hat(t, i) >= 0
                    p_mat(i, j) = p_U;
                elseif p_0 < p_L && B_hat(t, i) < 0
                    p_mat(i, j) = p_L;
                else
                    p_mat(i, j) = p_0;
                end

                if p_mat(i, j) > 0 && isfinite(p_mat(i, j))
                    server_cost_mat(i, j) = p.L / (p.omega * log2(1 + h * p_mat(i, j) / p.sigma));
                    server_E_mat(i, j) = p_mat(i, j) * server_cost_mat(i, j);
                    J_s(i, j) = (-B_hat(t, i) * p_mat(i, j) + p.V) * server_cost_mat(i, j);
                end
            end
        end

        [~, mode] = min([J_m(i), J_s(i, :), J_d]);
        if mode == 1
            chosen_mode(t, i) = 1;
            final_cost(t, i) = mobile_exe_cost(t, i);
            final_E(t, i) = mobile_exe_E(t, i);
        elseif mode == (M + 2)
            chosen_mode(t, i) = 3;
            final_cost(t, i) = p.phi;
            final_E(t, i) = 0;
        else
            chosen_mode(t, i) = 2;
            device_server_pairs = [device_server_pairs; [i, mode - 1, J_s(i, mode - 1)]];
        end
    end

    % Step 2: allocate MEC connections
    while ~isempty(device_server_pairs)
        if rand() <= p.eps
            [~, idx] = min(device_server_pairs(:, 3));
            i = device_server_pairs(idx, 1);
            j = device_server_pairs(idx, 2);
            if remained_connects(j) >= 1
                chosen_mode(t, i) = 2;
                ptx(t, i) = p_mat(i, j);
                final_cost(t, i) = server_cost_mat(i, j);
                final_E(t, i) = server_E_mat(i, j);
                remained_connects(j) = remained_connects(j) - 1;
                device_server_pairs(device_server_pairs(:, 1) == i, :) = [];
            else
                J_s(i, j) = inf;
                device_server_pairs(device_server_pairs(:, 1) == i, :) = [];
                if min(J_s(i, :)) ~= inf
                    [second_min_Js, second_min_j] = min(J_s(i, :));
                    device_server_pairs = [device_server_pairs; [i, second_min_j, second_min_Js]];
                else
                    [~, mode] = min([J_m(i), inf, J_d]);
                    if mode == 1
                        chosen_mode(t, i) = 1;
                        final_cost(t, i) = mobile_exe_cost(t, i);
                        final_E(t, i) = mobile_exe_E(t, i);
                    else
                        chosen_mode(t, i) = 3;
                        final_cost(t, i) = p.phi;
                        final_E(t, i) = 0;
                    end
                end
            end
        else
            for j = 1:M
                device_j_pairs = device_server_pairs(device_server_pairs(:, 2) == j, :);
                is = device_j_pairs(:, 1);
                if isempty(is)
                    continue;
                end

                if remained_connects(j) >= length(is)
                    chosen_mode(t, is) = 2;
                    ptx(t, is) = p_mat(is, j).';
                    final_cost(t, is) = server_cost_mat(is, j).';
                    final_E(t, is) = server_E_mat(is, j).';
                    remained_connects(j) = remained_connects(j) - length(is);
                    device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
                else
                    if remained_connects(j) == 0
                        device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
                        J_s(is, j) = inf;
                        for idx = 1:numel(is)
                            ii = is(idx);
                            [~, mode] = min([J_m(ii), J_s(ii, :), J_d]);
                            if mode == 1
                                chosen_mode(t, ii) = 1;
                                final_cost(t, ii) = mobile_exe_cost(t, ii);
                                final_E(t, ii) = mobile_exe_E(t, ii);
                            elseif mode == (M + 2)
                                chosen_mode(t, ii) = 3;
                                final_cost(t, ii) = p.phi;
                                final_E(t, ii) = 0;
                            else
                                chosen_mode(t, ii) = 2;
                                device_server_pairs = [device_server_pairs; [ii, mode - 1, J_s(ii, mode - 1)]];
                            end
                        end
                    else
                        [~, idxs] = sort(device_j_pairs(:, 3));
                        for idx = 1:remained_connects(j)
                            ii = device_j_pairs(idxs(idx), 1);
                            chosen_mode(t, ii) = 2;
                            ptx(t, ii) = p_mat(ii, j);
                            final_cost(t, ii) = server_cost_mat(ii, j);
                            final_E(t, ii) = server_E_mat(ii, j);
                            device_server_pairs(device_server_pairs(:, 1) == ii, :) = [];
                        end
                        remained_connects(j) = 0;

                        residual_is = device_server_pairs(device_server_pairs(:, 2) == j, 1);
                        J_s(residual_is, j) = inf;
                        for idx = 1:numel(residual_is)
                            ii = residual_is(idx);
                            [~, mode] = min([J_m(ii), J_s(ii, :), J_d]);
                            if mode == 1
                                chosen_mode(t, ii) = 1;
                                final_cost(t, ii) = mobile_exe_cost(t, ii);
                                final_E(t, ii) = mobile_exe_E(t, ii);
                            elseif mode == (M + 2)
                                chosen_mode(t, ii) = 3;
                                final_cost(t, ii) = p.phi;
                                final_E(t, ii) = 0;
                            else
                                chosen_mode(t, ii) = 2;
                                device_server_pairs = [device_server_pairs; [ii, mode - 1, J_s(ii, mode - 1)]];
                            end
                        end
                    end
                end
            end
        end
    end

    B(t + 1, :) = B(t, :) - final_E(t, :) + e(t, :);
end

% Build network-level metrics
avg_network_cost = zeros(T, 1);
sum_cost = 0;
sum_req = 0;
mode_count = [0, 0, 0]; % mobile, server, drop
ratio_ts = zeros(T, 3);

for t = 1:T
    requested = (chosen_mode(t, :) ~= 4);
    sum_cost = sum_cost + sum(final_cost(t, requested));
    sum_req = sum_req + sum(requested);
    if sum_req > 0
        avg_network_cost(t) = sum_cost / sum_req;
    end

    mode_count(1) = mode_count(1) + sum(chosen_mode(t, requested) == 1);
    mode_count(2) = mode_count(2) + sum(chosen_mode(t, requested) == 2);
    mode_count(3) = mode_count(3) + sum(chosen_mode(t, requested) == 3);
    if sum_req > 0
        ratio_ts(t, :) = mode_count / sum_req;
    end
end

result.avg_network_cost = avg_network_cost;
result.final_avg_cost = avg_network_cost(end);
result.final_ratio = ratio_ts(end, :);
result.ratio_ts = ratio_ts;
end

function root = my_bisection(func, lower_bound, upper_bound, tol)
if lower_bound >= upper_bound
    root = inf;
    return;
end

a = lower_bound;
b = upper_bound;
fa = func(a);
fb = func(b);

if abs(fa) < tol
    root = a;
    return;
end
if abs(fb) < tol
    root = b;
    return;
end
if fa * fb > 0
    root = inf;
    return;
end

while (b - a) > tol
    c = (a + b) / 2;
    fc = func(c);
    if abs(fc) < tol
        root = c;
        return;
    elseif fa * fc < 0
        b = c;
    else
        a = c;
        fa = fc;
    end
end

root = (a + b) / 2;
end
