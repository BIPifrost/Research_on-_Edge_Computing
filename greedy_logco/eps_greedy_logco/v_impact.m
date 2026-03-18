clc, clear
opt = optimset('Display', 'none');

err_fsolve_pEmin = [];
err_fsolve_pEmax = [];
err_fsolve_p0 = [];

%% basic parameter settings (had better not change those paras)
k = 1e-28;                % effective switched capacitance (a constant decided by the chip architecture)
tau = 0.002;              % the length of time slot (in second)
phi = 0.002;              % the cost of task dropping (in second)
omega = 1e6;              % the bandwidth of MEC server (in Hz)
sigma = 1e-13;            % the noise power of the receiver (in W)
p_tx_max = 1;             % the maximum transmit power of mobile device (in W)
f_max = 1.5e9;            % the maximum CPU-cycle frequency of mobile device (in Hz)
E_max = 0.002;            % the maximum amout of battery output energy (in J)
L = 1000;                 % the input size of the computation task (in bit)
X = 737.5;                % the number of CPU cycles needed on processing one bit of task
W = L * X;                % the number of CPU cycles needed on processing one task
E_H_max = 48e-6;          % the upper bound of the energy arrive at the mobile device (in J)
p_H = E_H_max / (2*tau);  % the average Energy Harvesting (EH) power (in W)
g0 = power(10, -4);       % the path-loss constant
d0 = 1;                   % the relative distance between each mobile device and each MEC server

%% parameter control
N = 10;                   % the number of mobile devices
M = 8;                    % the number of MEC servers
T = 2000;                 % the number of time slot (a.k.a. the size of the time horizon)
tau_d = 0.002;            % execution deadline (in second)
d = 50;                   % the distance between the mobile device and the MEC server (in meter)
E_min = 0.02e-3;          % the minimum amout of battery output energy (in J)
rho = 0.6;                % the probability that the computation task is requested
max_connects = 4;         % the maximum number of processible mobile devices for each MEC server ($ \frac{f_s^{max} \tau}{L X} $)
min_distance = 10;        % the minimum distance from mobile device to MEC server
max_distance = 50;        % the maximum distance from mobile device to MEC server
eps = 0.8;                % if rand() <= eps, for the best i-j pair, we always choose MEC server execution mode (except no MEC server to choose)

V_list = [1e-6, 1e-5, 5e-5, 1e-4]; 
num_V = length(V_list);

avg_latency_V = zeros(num_V, 1);
avg_energy_V = zeros(num_V, 1);
battery_history_V = zeros(T, num_V); 
simulation_time_V = zeros(num_V, 1);

for v_idx = 1:num_V
    V = V_list(v_idx);
    fprintf('\n======================================================\n');
    fprintf('正在进行仿真测试，当前控制参�?V = %e ...\n', V);
    
    tic;

    E_max_hat = min(max(k * W * (f_max)^2, p_tx_max * tau), E_max);
    theta = E_max_hat + V * phi / E_min;

    B = zeros(T, N);                   
    B_hat = zeros(T, N);               
    e = zeros(T, N);                   
    chosen_mode = zeros(T, N);         
    chosen_server = zeros(T, N);       
    f = zeros(T, N);                   
    p = zeros(T, N);                   

    mobile_exe_cost = zeros(T, N);     
    server_exe_cost = zeros(T, N);     
    final_chosen_cost = zeros(T, N);   

    mobile_exe_E = zeros(T, N);        
    server_exe_E = zeros(T, N);        
    final_chosen_E = zeros(T, N);      

    t = 1;
    while t <= T
        
        device_server_pairs = [];      
        remained_connects = max_connects * ones(M, 1);    
        J_m = zeros(N, 1); J_s = zeros(N, M);             
        p_mat = zeros(N, M);                              
        server_cost_mat = zeros(N, M);                    
        server_E_mat = zeros(N, M);                       
        
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
            else
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
                    f(t, i) = 0;
                    mobile_exe_cost(t, i) = 0;
                    mobile_exe_E(t, i) = 0;
                    J_m(i) = inf;
                end
                
                for j = 1: M
                    h = h_mat(i, j);
                    E_tmp = sigma * L * log(2) / (omega * h);
                    p_L_taud = (power(2, L / (omega * tau_d)) - 1) * sigma / h;
                    
                    if E_tmp >= E_min
                        p_L = p_L_taud;
                    else
                        y = @(x) x * L - omega * log2(1 + h*x/sigma) * E_min;
                        [p_E_min, ~, exitflag] = fsolve(y, 0.2, opt);
                        err_fsolve_pEmin = [err_fsolve_pEmin, abs(y(p_E_min))];
                        p_L = max(p_L_taud, p_E_min);
                    end
                    
                    if E_tmp >= E_max
                        p_U = 0;
                    else
                        y = @(x) x * L - omega * log2(1 + h*x/sigma) * E_max;
                        [p_E_max, ~, ~] = fsolve(y, 100, opt);
                        err_fsolve_pEmax = [err_fsolve_pEmax, abs(y(p_E_max))];
                        p_U = min(p_tx_max, p_E_max);
                    end
                    
                    if p_L <= p_U
                        virtual_battery = B_hat(t, i);
                        y = @(x) virtual_battery * log2(1 + h*x/sigma) + ...
                            h * (V - virtual_battery*x) / log(2) / (sigma + h*x);
                        [p_0, ~, exitflag] = fsolve(y, 0.5, opt);
                        err_fsolve_p0 = [err_fsolve_p0, abs(y(p_0))];

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
                        p_mat(i, j) = 0;
                        server_cost_mat(i, j) = 0;
                        server_E_mat(i, j) = 0;
                        J_s(i, j) = inf;
                    end
                end
                
                J_d = V * phi;
                [~, mode] = min([J_m(i), J_s(i, :), J_d]);
                if mode == 1
                    chosen_mode(t, i) = 1;
                    final_chosen_cost(t, i) = mobile_exe_cost(t, i);
                    final_chosen_E(t, i) = mobile_exe_E(t, i);
                elseif mode == (M+2)
                    chosen_mode(t, i) = 3;
                    final_chosen_cost(t, i) = phi;
                    final_chosen_E(t, i) = 0;
                else
                    chosen_mode(t, i) = 2;
                    device_server_pairs = [device_server_pairs; [i, mode-1, J_s(i, mode-1)]];
                end
            end
        end
        
        while ~isempty(device_server_pairs)
            if rand() <= eps
                [min_Js, idx] = min(device_server_pairs(:, 3));
                i = device_server_pairs(idx, 1); j = device_server_pairs(idx, 2);    
                if remained_connects(j) >= 1
                    chosen_mode(t, i) = 2;          
                    p(t, i) = p_mat(i, j);
                    server_exe_cost(t, i) = server_cost_mat(i, j);
                    server_exe_E(t, i) = server_E_mat(i, j);
                    chosen_server(t, i) = j;
                    final_chosen_cost(t, i) = server_exe_cost(t, i);
                    final_chosen_E(t, i) = server_exe_E(t, i);
                    
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
                            final_chosen_cost(t, i) = mobile_exe_cost(t, i);
                            final_chosen_E(t, i) = mobile_exe_E(t, i);
                        else
                            final_chosen_cost(t, i) = phi;
                            final_chosen_E(t, i) = 0;
                        end
                    end
                end
            else
                for j = 1: M
                    device_j_pairs = device_server_pairs(device_server_pairs(:, 2) == j, :);
                    is = device_j_pairs(:, 1);
                    if isempty(is)
                        continue;
                    end
                    if remained_connects(j) >= length(is)
                        chosen_mode(t, is) = 2;     
                        p(t, is) = transp(p_mat(is, j));
                        server_exe_cost(t, is) = transp(server_cost_mat(is, j));
                        server_exe_E(t, is) = transp(server_E_mat(is, j));
                        chosen_server(t, is) = repmat(j, 1, length(is));
                        final_chosen_cost(t, is) = server_exe_cost(t, is);
                        final_chosen_E(t, is) = server_exe_E(t, is);

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
                                    chosen_mode(t, i) = 1;
                                    final_chosen_cost(t, i) = mobile_exe_cost(t, i);
                                    final_chosen_E(t, i) = mobile_exe_E(t, i);
                                elseif mode == (M+2)
                                    chosen_mode(t, i) = 3;
                                    final_chosen_cost(t, i) = phi;
                                    final_chosen_E(t, i) = 0;
                                else
                                    chosen_mode(t, i) = 2;
                                    device_server_pairs = [device_server_pairs; [i, mode-1, J_s(i, mode-1)]];
                                end
                            end
                        else
                            [~, idxs] = sort(device_j_pairs(:, 3));
                            for idx = 1: remained_connects(j)
                                i = device_j_pairs(idxs(idx), 1);
                                chosen_mode(t, i) = 2;
                                p(t, i) = p_mat(i, j);
                                server_cost_mat_val = server_cost_mat(i, j);
                                server_exe_cost(t, i) = server_cost_mat_val;
                                server_exe_E(t, i) = server_E_mat(i, j);
                                chosen_server(t, i) = j;
                                final_chosen_cost(t, i) = server_exe_cost(t, i);
                                final_chosen_E(t, i) = server_exe_E(t, i);

                                device_server_pairs(device_server_pairs(:, 1) == i, :) = [];
                            end

                            remained_connects(j) = 0;
                            residual_is = device_server_pairs(device_server_pairs(:, 2) == j, 1);
                            J_s(residual_is, j) = inf;
                            for idx = 1: numel(residual_is)
                                residual_i = residual_is(idx);
                                [~, mode] = min([J_m(residual_i), J_s(residual_i, :), J_d]);
                                if mode == 1
                                    chosen_mode(t, residual_i) = 1;
                                    final_chosen_cost(t, residual_i) = mobile_exe_cost(t, residual_i);
                                    final_chosen_E(t, residual_i) = mobile_exe_E(t, residual_i);
                                elseif mode == (M+2)
                                    chosen_mode(t, residual_i) = 3;
                                    final_chosen_cost(t, residual_i) = phi;
                                    final_chosen_E(t, residual_i) = 0;
                                else
                                    chosen_mode(t, residual_i) = 2;
                                    device_server_pairs = [device_server_pairs; [residual_i, mode-1, J_s(residual_i, mode-1)]];
                                end
                            end
                        end
                    end
                end
            end
        end
        
        if t < T
            B(t + 1, :) = B(t, :) - final_chosen_E(t, :) + e(t, :);
        end
        t = t + 1;
    end
    
    simulation_time_V(v_idx) = toc;
    fprintf('测试完毕！V=
    
    total_valid_requests = sum(chosen_mode(:) ~= 4);
    
    if total_valid_requests > 0
        avg_latency_V(v_idx) = sum(final_chosen_cost(chosen_mode ~= 4)) / total_valid_requests;
        avg_energy_V(v_idx) = sum(final_chosen_E(chosen_mode ~= 4)) / total_valid_requests;
    else
        avg_latency_V(v_idx) = 0;
        avg_energy_V(v_idx) = 0;
    end
    
    battery_history_V(:, v_idx) = mean(B, 2);
end



figure('Name', 'Average Latency vs V');
plot(V_list, avg_latency_V, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 8);
grid on;
xlabel('Control Parameter $V$', 'Interpreter', 'latex');
ylabel('Average Execution Cost/Latency (s)');
title('Impact of Control Parameter $V$ on Average Latency');

figure('Name', 'Average Energy vs V');
plot(V_list, avg_energy_V, '-s', 'LineWidth', 2, 'Color', '#D95319', 'MarkerFaceColor', '#D95319', 'MarkerSize', 8);
grid on;
xlabel('Control Parameter $V$', 'Interpreter', 'latex');
ylabel('Average Energy Consumption (J)');
title('Impact of Control Parameter $V$ on Average Energy');

figure('Name', 'Battery Evolution vs V');
hold on;
colors = lines(num_V);
for v_idx = 1:num_V
    plot(1:T, battery_history_V(:, v_idx), 'LineWidth', 1.5, 'Color', colors(v_idx, :));
end
grid on;
xlabel('Time Slot $t$', 'Interpreter', 'latex');
ylabel('Average Battery Energy Level (J)');
title('Battery Evolution under Different Control Parameter $V$');
legend_str = cellstr(num2str(V_list', 'V = %e'));
legend(legend_str, 'Location', 'best');

figure('Name', 'Simulation Time Cost');
b = bar(categorical(string(V_list)), simulation_time_V, 0.5);
b.FaceColor = '#EDB120';
grid on;
xlabel('Control Parameter $V$');
ylabel('Algorithm Run Time (Seconds)');
title('Simulation Computational Time Cost for Each $V$');
for i = 1:num_V
    text(i, simulation_time_V(i), sprintf('%.2fs', simulation_time_V(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
end

disp(' ');
disp('================ fsolve 求解误差分析 (全过程汇�? ================');
disp(['p_E_min 误差方差: ', num2str(var(err_fsolve_pEmin, 'omitnan')), '  | 平均误差: ', num2str(mean(err_fsolve_pEmin, 'omitnan'))]);
disp(['p_E_max 误差方差: ', num2str(var(err_fsolve_pEmax, 'omitnan')), '  | 平均误差: ', num2str(mean(err_fsolve_pEmax, 'omitnan'))]);
disp(['p_0     误差方差: ', num2str(var(err_fsolve_p0, 'omitnan')),    '  | 平均误差: ', num2str(mean(err_fsolve_p0, 'omitnan'))]);
disp('==================================================================');
