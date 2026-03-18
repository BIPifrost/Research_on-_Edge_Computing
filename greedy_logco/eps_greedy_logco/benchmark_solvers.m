clc; clear; close all;
disp('==================================================');
disp('开始严格基准测�? 原版 fsolve vs 优化�?bisection');
disp('==================================================');

T_test = 2000;
global_seed = 42;

disp(['[1/2] 正在运行原版 fsolve 方法 (T=', num2str(T_test), '), 请耐心等待...']);
[mode_fsolve, time_fsolve] = run_simulation(false, T_test, global_seed);
disp(['原版运行完毕，耗时: ', num2str(time_fsolve), ' �?]);

disp('--------------------------------------------------');
disp(['[2/2] 正在运行优化�?my_bisection 方法 (T=', num2str(T_test), ')...']);
[mode_bisect, time_bisect] = run_simulation(true, T_test, global_seed);
disp(['优化版运行完毕，耗时: ', num2str(time_bisect), ' �?]);

disp('==================================================');
disp('最终对比报�?');
disp(['1. 速度提升倍数: ', num2str(time_fsolve / time_bisect), ' �?]);

task_generated_idx = (mode_fsolve ~= 4); 
total_decisions = sum(task_generated_idx, 'all');
diff_count = sum(mode_fsolve(task_generated_idx) ~= mode_bisect(task_generated_idx), 'all');
consistency_rate = 100 * (1 - diff_count / total_decisions);

disp(['2. 总决策次�?(剔除无任�?: ', num2str(total_decisions)]);
disp(['3. 实际决策不同的数�? ', num2str(diff_count)]);
disp(['   决策一致率: ', num2str(consistency_rate), ' %']);
disp('==================================================');

%% =========================================================================
function [chosen_mode, time_cost] = run_simulation(use_bisection, T, global_seed)
    opt = optimset('Display', 'none');
    
    k = 1e-28; tau = 0.002; phi = 0.002; omega = 1e6; sigma = 1e-13;
    p_tx_max = 1; f_max = 1.5e9; E_max = 0.002; L = 1000; X = 737.5;
    W = L * X; E_H_max = 48e-6; g0 = power(10, -4); d0 = 1;
    N = 10; M = 8; tau_d = 0.002; E_min = 0.02e-3; V = 1e-5;
    rho = 0.6; max_connects = 4; min_distance = 10; max_distance = 50; eps = 0.8;
    
    E_max_hat = min(max(k * W * (f_max)^2, p_tx_max * tau), E_max);
    theta = E_max_hat + V * phi / E_min;

    B = zeros(T, N); B_hat = zeros(T, N); e = zeros(T, N);
    chosen_mode = zeros(T, N); chosen_server = zeros(T, N);
    f = zeros(T, N); p = zeros(T, N);
    mobile_exe_cost = zeros(T, N); server_exe_cost = zeros(T, N); final_chosen_cost = zeros(T, N);
    mobile_exe_E = zeros(T, N); server_exe_E = zeros(T, N); final_chosen_E = zeros(T, N);

    tic;
    
    for t = 1: T
        rng(global_seed + t); 
        
        device_server_pairs = []; 
        remained_connects = max_connects * ones(M, 1);
        J_m = zeros(N, 1); J_s = zeros(N, M);
        p_mat = zeros(N, M); server_cost_mat = zeros(N, M); server_E_mat = zeros(N, M);
        
        B_hat(t, :) = B(t, :) - theta;
        distances = unifrnd(min_distance, max_distance, N, M);
        gamma = exprnd(1, N, M);
        h_mat = g0 * gamma .* power(d0 ./ distances, 4);
        
        for i = 1: N
            E_H_t = unifrnd(0, E_H_max);
            if B_hat(t, i) <= 0
                e(t, i) = E_H_t;
            end
            
            zeta = binornd(1, rho);
            if zeta == 0
                chosen_mode(t, i) = 4;
                continue;
            end
            
            f_L = max(sqrt(E_min / (k * W)), W / tau_d);
            f_U = min(sqrt(E_max / (k * W)), f_max);
            if f_L <= f_U
                if B_hat(t, i) < 0
                    f_0 = (V / (-2 * B_hat(t, i) * k))^(1/3);
                else
                    f_0 = -(V / (2 * B_hat(t, i) * k))^(1/3);
                end
                if (f_0 > f_U && B_hat(t, i) < 0) || (B_hat(t, i) >= 0)
                    f(t, i) = f_U;
                elseif f_0 >= f_L && f_0 <= f_U && B_hat(t, i) < 0
                    f(t, i) = f_0;
                elseif f_0 < f_L && B_hat(t, i) < 0
                    f(t, i) = f_L;
                end
                mobile_exe_cost(t, i) = W / f(t, i);
                mobile_exe_E(t, i) = k * W * (f(t, i)^2);
                J_m(i) = -B_hat(t, i) * k * W * f(t, i)^2 + V * W / f(t, i); 
            else
                f(t, i) = 0; mobile_exe_cost(t, i) = 0; mobile_exe_E(t, i) = 0;
                J_m(i) = inf;
            end
            
            for j = 1: M
                h = h_mat(i, j);
                E_tmp = sigma * L * log(2) / (omega * h);
                p_L_taud = (power(2, L / (omega * tau_d)) - 1) * sigma / h;
                
                if E_tmp >= E_min
                    p_L = p_L_taud;
                else
                    y_E_min = @(x) x * L - omega * log2(1 + h*x/sigma) * E_min;
                    if use_bisection
                        a_min = max(0, (omega * E_min) / (L * log(2)) - sigma / h);
                        if y_E_min(p_tx_max) <= 0
                            p_E_min = inf;
                        else
                            p_E_min = my_bisection(y_E_min, a_min, p_tx_max, 1e-7);
                        end
                        
                        if isinf(p_E_min) || isnan(p_E_min)
                            [p_E_min, ~, flag1] = fsolve(y_E_min, 0.2, opt);
                            if flag1 <= 0 || p_E_min < 0 || ~isreal(p_E_min), p_E_min = inf; end
                        end
                    else
                        [p_E_min, ~, flag1] = fsolve(y_E_min, 0.2, opt);
                        if flag1 <= 0 || p_E_min < 0 || ~isreal(p_E_min), p_E_min = inf; end
                    end
                    p_L = max(p_L_taud, p_E_min);
                end
                
                if E_tmp >= E_max
                    p_U = 0;
                else
                    y_E_max = @(x) x * L - omega * log2(1 + h*x/sigma) * E_max;
                    if use_bisection
                        a_max = max(0, (omega * E_max) / (L * log(2)) - sigma / h);
                        if y_E_max(p_tx_max) <= 0
                            p_E_max = inf;
                        else
                            p_E_max = my_bisection(y_E_max, a_max, p_tx_max, 1e-7);
                        end
                        
                        if isinf(p_E_max) || isnan(p_E_max)
                            [p_E_max, ~, flag2] = fsolve(y_E_max, 100, opt);
                            if flag2 <= 0 || p_E_max < 0 || ~isreal(p_E_max), p_E_max = inf; end
                        end
                    else
                        [p_E_max, ~, flag2] = fsolve(y_E_max, 100, opt);
                        if flag2 <= 0 || p_E_max < 0 || ~isreal(p_E_max), p_E_max = inf; end
                    end
                    p_U = min(p_tx_max, p_E_max);
                end
                
                p_0 = inf;
                if p_L <= p_U
                    virtual_battery = B_hat(t, i);
                    y_p0 = @(x) virtual_battery * log2(1 + h*x/sigma) + ...
                        h * (V - virtual_battery*x) / log(2) / (sigma + h*x);
                    
                    % ==========================================================
                    [p_0_fsolve_test, fval_test, exitflag_test] = fsolve(y_p0, 0.5, opt);
                    p_0_bi_test = my_bisection(y_p0, 0, 100, 1e-7);
                    
                    if abs(y_p0(p_0_fsolve_test)) > 1e-4 && abs(y_p0(p_0_bi_test)) < 1e-4
                        disp(['[抓到 Bug!] 时隙 t=', num2str(t), ', 设备 i=', num2str(i)]);
                        disp([' -> [原版 fsolve] �? ', num2str(p_0_fsolve_test), ' | 误差: ', num2str(y_p0(p_0_fsolve_test))]);
                        disp([' -> [优化 二分法] �? ', num2str(p_0_bi_test),     ' | 误差: ', num2str(y_p0(p_0_bi_test))]);
                    end
                    % ==========================================================

                    if use_bisection
                        p_0 = my_bisection(y_p0, 0, 100, 1e-7);
                        if isinf(p_0) || isnan(p_0)
                            [p_0, ~, flag3] = fsolve(y_p0, 0.5, opt);
                            if flag3 <= 0 || p_0 < 0 || ~isreal(p_0), p_0 = inf; end
                        end
                    else
                        [p_0, ~, flag3] = fsolve(y_p0, 0.5, opt);
                        if flag3 <= 0 || p_0 < 0 || ~isreal(p_0), p_0 = inf; end
                    end
                    
                    if (p_U < p_0 && B_hat(t, i) < 0) || B_hat(t, i) >= 0
                        p_mat(i, j) = p_U;
                    elseif p_0 < p_L && B_hat(t, i) < 0
                        p_mat(i, j) = p_L;
                    elseif p_0 >= p_L && p_0 <= p_U && B_hat(t, i) < 0
                        p_mat(i, j) = p_0;
                    end
                    
                    server_cost_mat(i, j) = L / (omega * log2(1 + h*p_mat(i, j)/sigma));
                    server_E_mat(i, j) = p_mat(i, j) * server_cost_mat(i, j);
                    J_s(i, j) = (-B_hat(t, i) * p_mat(i, j) + V) * server_cost_mat(i, j);
                else
                    p_mat(i, j) = 0; server_cost_mat(i, j) = 0; server_E_mat(i, j) = 0;
                    J_s(i, j) = inf;
                end
            end
            
            J_d = V * phi;
            [~, mode] = min([J_m(i), J_s(i, :), J_d]);
            if mode == 1
                chosen_mode(t, i) = 1;
                final_chosen_cost(t, i) = mobile_exe_cost(t, i); final_chosen_E(t, i) = mobile_exe_E(t, i);
            elseif mode == (M+2)
                chosen_mode(t, i) = 3;
                final_chosen_cost(t, i) = phi; final_chosen_E(t, i) = 0;
            else
                chosen_mode(t, i) = 2;
                device_server_pairs = [device_server_pairs; [i, mode-1, J_s(i, mode-1)]];
            end
        end
        
        while ~isempty(device_server_pairs)
            if rand() <= eps
                [min_Js, idx] = min(device_server_pairs(:, 3));
                i = device_server_pairs(idx, 1); j = device_server_pairs(idx, 2);
                if remained_connects(j) >= 1
                    chosen_mode(t, i) = 2; p(t, i) = p_mat(i, j);
                    server_exe_cost(t, i) = server_cost_mat(i, j); server_exe_E(t, i) = server_E_mat(i, j);
                    chosen_server(t, i) = j; final_chosen_cost(t, i) = server_exe_cost(t, i); final_chosen_E(t, i) = server_exe_E(t, i);
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
                        chosen_mode(t, i) = mode;
                        if mode == 1
                            final_chosen_cost(t, i) = mobile_exe_cost(t, i); final_chosen_E(t, i) = mobile_exe_E(t, i);
                        else
                            final_chosen_cost(t, i) = phi; final_chosen_E(t, i) = 0;
                        end
                    end
                end
            else
                for j = 1: M
                    device_j_pairs = device_server_pairs(device_server_pairs(:, 2) == j, :);
                    is = device_j_pairs(:, 1);
                    if isempty(is), continue; end
                    if remained_connects(j) >= length(is)
                        chosen_mode(t, is) = 2; p(t, is) = transp(p_mat(is, j)); 
                        server_exe_cost(t, is) = transp(server_cost_mat(is, j)); server_exe_E(t, is) = transp(server_E_mat(is, j));
                        chosen_server(t, is) = repmat(j, 1, length(is)); 
                        final_chosen_cost(t, is) = server_exe_cost(t, is); final_chosen_E(t, is) = server_exe_E(t, is);
                        remained_connects(j) = remained_connects(j) - length(is);
                        device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
                    else
                        if remained_connects(j) == 0
                            device_server_pairs(device_server_pairs(:, 2) == j, :) = [];
                            J_s(is, j) = inf;
                            for idx = 1: numel(is)
                                i = is(idx);
                                [~, mode] = min([J_m(i), J_s(i, :), J_d]);
                                if mode == 1
                                    chosen_mode(t, i) = 1; final_chosen_cost(t, i) = mobile_exe_cost(t, i); final_chosen_E(t, i) = mobile_exe_E(t, i);
                                elseif mode == (M+2)
                                    chosen_mode(t, i) = 3; final_chosen_cost(t, i) = phi; final_chosen_E(t, i) = 0;
                                else
                                    chosen_mode(t, i) = 2; device_server_pairs = [device_server_pairs; [i, mode-1, J_s(i, mode-1)]];
                                end
                            end
                        else
                            [~, idxs] = sort(device_j_pairs(:, 3));
                            sorted_pairs = device_j_pairs(idxs, :); 
                            
                            for conn_idx = 1: remained_connects(j)
                                i = sorted_pairs(conn_idx, 1);
                                chosen_mode(t, i) = 2; p(t, i) = p_mat(i, j); 
                                server_exe_cost(t, i) = server_cost_mat(i, j); server_exe_E(t, i) = server_E_mat(i, j);
                                chosen_server(t, i) = j; final_chosen_cost(t, i) = server_exe_cost(t, i); final_chosen_E(t, i) = server_exe_E(t, i);
                                device_server_pairs(device_server_pairs(:, 1) == i, :) = [];
                            end
                            remained_connects(j) = 0;
                            residual_is = device_server_pairs(device_server_pairs(:, 2) == j, 1);
                            J_s(residual_is, j) = inf;
                            for idx = 1: numel(residual_is)
                                residual_i = residual_is(idx);
                                [~, mode] = min([J_m(residual_i), J_s(residual_i, :), J_d]);
                                if mode == 1
                                    chosen_mode(t, residual_i) = 1; final_chosen_cost(t, residual_i) = mobile_exe_cost(t, residual_i); final_chosen_E(t, residual_i) = mobile_exe_E(t, residual_i);
                                elseif mode == (M+2)
                                    chosen_mode(t, residual_i) = 3; final_chosen_cost(t, residual_i) = phi; final_chosen_E(t, residual_i) = 0;
                                else
                                    chosen_mode(t, residual_i) = 2; device_server_pairs = [device_server_pairs; [residual_i, mode-1, J_s(residual_i, mode-1)]];
                                end
                            end
                        end
                    end
                end
            end
        end
        B(t + 1, :) = B(t, :) - final_chosen_E(t, :) + e(t, :);
    end
    
    time_cost = toc;
end

%% =========================================================================
function root = my_bisection(func, lower_bound, upper_bound, tol)
    if lower_bound >= upper_bound
        root = inf; 
        return; 
    end

    a = lower_bound; 
    b = upper_bound; 
    
    fa = func(a); 
    fb = func(b);
    
    if abs(fa) < 1e-6, root = a; return; end
    if abs(fb) < 1e-6, root = b; return; end
    
    if fa * fb > 0
        root = inf; 
        return; 
    end
    
    while (b - a) > 1e-12
        c = (a + b) / 2; 
        fc = func(c);
        if abs(fc) < 1e-12
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
