clc, clear
opt = optimset('Display', 'none');

err_fsolve_pEmin = [];
err_fsolve_pEmax = [];
err_fsolve_p0 = [];

%% basic parameter settings 
k = 1e-28;                
tau = 0.002;              
phi = 0.002;
omega = 1e6;              
sigma = 1e-13;            
p_tx_max = 1;             
f_max = 1.5e9;            
E_max = 0.002;            
L = 1000;                 
X = 737.5;                
W = L * X;                
E_H_max = 48e-6;          
p_H = E_H_max / (2*tau);  
g0 = power(10, -4);       
d0 = 1;                   

%% parameter control
N = 10;                   
M = 8;                    
T = 1000;
tau_d = 0.002;            
d = 50;                   
E_min = 0.02e-3;          
V = 1e-5;                 
max_connects = 4;         
min_distance = 10;        
max_distance = 50;        
eps_greedy = 0.8;         

rho_list = [0.2, 0.4, 0.6, 0.8, 1.0];
num_rho = length(rho_list);

compare_avg_cost = zeros(num_rho, 1);
compare_avg_battery = zeros(num_rho, 1);

ratio_local = zeros(num_rho, 1);
ratio_mec   = zeros(num_rho, 1);
ratio_drop  = zeros(num_rho, 1);

E_max_hat = min(max(k * W * (f_max)^2, p_tx_max * tau), E_max);
theta = E_max_hat + V * phi / E_min;

for r_idx = 1:num_rho
    rho = rho_list(r_idx);
    disp(['==================================================']);
    disp(['====> 开始仿�? 当前 rho = ', num2str(rho), ' (', num2str(r_idx), '/', num2str(num_rho), ') <====']);

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
        if mod(t, 200) == 0
            disp(['  -> 正在处理 Time slot #', num2str(t), ' / ', num2str(T)]);
        end
        
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
                % mobile execution
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
                
                % MEC execution
                for j = 1: M
                    h = h_mat(i, j);
                    E_tmp = sigma * L * log(2) / (omega * h);
                    p_L_taud = (power(2, L / (omega * tau_d)) - 1) * sigma / h;
                    
                    if E_tmp >= E_min
                        p_L = p_L_taud;
                    else
                        y1 = @(x) x * L - omega * log2(1 + h*x/sigma) * E_min;
                        [p_E_min, ~, ~] = fsolve(y1, 0.2, opt);
                        p_L = max(p_L_taud, p_E_min);
                    end
                    
                    if E_tmp >= E_max
                        p_U = 0;
                    else
                        y2 = @(x) x * L - omega * log2(1 + h*x/sigma) * E_max;
                        [p_E_max, ~, ~] = fsolve(y2, 100, opt);
                        p_U = min(p_tx_max, p_E_max);
                    end
                    
                    if p_L <= p_U
                        virtual_battery = B_hat(t, i);
                        y3 = @(x) virtual_battery * log2(1 + h*x/sigma) + ...
                            h * (V - virtual_battery*x) / log(2) / (sigma + h*x);
                        [p_0, ~, ~] = fsolve(y3, 0.5, opt);
                        
                        if (p_U < p_0 && B_hat(t, i) < 0) || B_hat(t, i) >= 0
                            p_mat(i, j) = p_U;
                        elseif p_0 < p_L && B_hat(t) < 0
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
            if rand() <= eps_greedy
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
                    if isempty(is), continue; end
                    
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
                                server_exe_cost(t, i) = server_cost_mat(i, j);
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
        
        %% step 3: update B
        B(t + 1, :) = B(t, :) - final_chosen_E(t, :) + e(t, :);
        t = t + 1;
    end
    
    total_requests = sum(sum(chosen_mode ~= 4));
    
    total_local = sum(sum(chosen_mode == 1));
    total_mec   = sum(sum(chosen_mode == 2));
    total_drop  = sum(sum(chosen_mode == 3));
    
    ratio_local(r_idx) = total_local / total_requests;
    ratio_mec(r_idx)   = total_mec / total_requests;
    ratio_drop(r_idx)  = total_drop / total_requests;
    
    compare_avg_cost(r_idx) = sum(sum(final_chosen_cost)) / total_requests;
    compare_avg_battery(r_idx) = mean(mean(B(2:end, :)));
end

disp(' ');
disp('================ 仿真结束，正在绘制对比视�?================');

x_labels = cell(1, num_rho);
for k = 1:num_rho
    x_labels{k} = ['\rho = ', num2str(rho_list(k))];
end

figure('Name', 'Mode Ratio vs \rho', 'Position', [100, 100, 500, 400]);
b = bar(1:num_rho, [ratio_local, ratio_mec, ratio_drop], 'stacked');
b(1).FaceColor = [0.4660 0.6740 0.1880];
b(2).FaceColor = [0.3010 0.7450 0.9330];
b(3).FaceColor = [0.8500 0.3250 0.0980];
set(gca, 'XTickLabel', x_labels);
title('Execution Mode Distribution vs. Task Request Probability (\rho)');
ylabel('Ratio of Total Requests');
legend('Local Execution', 'MEC Server Execution', 'Drop', 'Location', 'best');
grid on;

figure('Name', 'Avg Cost vs \rho', 'Position', [650, 100, 450, 400]);
bar(1:num_rho, compare_avg_cost, 0.5, 'FaceColor', [0.5 0.5 0.8]);
set(gca, 'XTickLabel', x_labels);
title('Overall Average Execution Cost');
ylabel('Average Delay / Penalty (s)');
grid on;

figure('Name', 'Avg Battery vs \rho', 'Position', [1150, 100, 450, 400]);
bar(1:num_rho, compare_avg_battery, 0.5, 'FaceColor', [0.8 0.6 0.2]);
set(gca, 'XTickLabel', x_labels);
title('Overall Average Battery Energy Level');
ylabel('Energy (J)');
grid on;
