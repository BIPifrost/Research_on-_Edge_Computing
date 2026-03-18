
opt = optimset('Display', 'none');


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


N = 10;
M = 8;
T = 200;
tau_d = 0.002;
d = 50;
E_min = 0.02e-3;
V = 1e-5;
rho = 0.6;
max_connects = 4;
min_distance = 10;
max_distance = 50;


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
    disp(['===> Time slot #', num2str(t), ' <==='])
    

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
        disp(['Mobile device #', num2str(i)])
        

        E_H_t = unifrnd(0, E_H_max);
        if B_hat(t, i) <= 0
            e(t, i) = E_H_t;
        end
        


        zeta = binornd(1, rho);
        if zeta == 0

            disp('no task request generated!')
            chosen_mode(t, i) = 4;
        else

            disp('task request generated!')
            


            f_L = max(sqrt(E_min / (k * W)), W / tau_d);
            f_U = min(sqrt(E_max / (k * W)), f_max);
            if f_L <= f_U

                disp('mobile execution ($\mathcal{P}_{ME}$) is feasible!')
                
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

                if f(t, i) == 0
                    disp('Something wrong! f is 0!')
                end
                

                mobile_exe_cost(t, i) = W / f(t, i);

                mobile_exe_E(t, i) = k * W * (f(t, i)^2);

                J_m(i) = -B_hat(t, i) * k * W * (f(t, i)^2 + V * W / f(t, i)); 
            else





                disp('mobile execution ($\mathcal{P}_{ME}$) is not feasible!')
                f(t, i) = 0;
                mobile_exe_cost(t, i) = 0;
                mobile_exe_E(t, i) = 0;

                J_m(i) = inf;
            end
            


            for j = 1: M
                disp(['MEC server #', num2str(j)])
                h = h_mat(i, j);
                
                E_tmp = sigma * L * log(2) / (omega * h);
                p_L_taud = (power(2, L / (omega * tau_d)) - 1) * sigma / h;

                if E_tmp >= E_min
                    p_L = p_L_taud;
                else

                    y = @(x) x * L - omega * log2(1 + h*x/sigma) * E_min;


                    p_E_min = fsolve(y, 0.2, opt);
                    p_L = max(p_L_taud, p_E_min);
                end

                if E_tmp >= E_max
                    p_U = 0;
                else

                    y = @(x) x * L - omega * log2(1 + h*x/sigma) * E_max;


                    p_E_max = fsolve(y, 100, opt);
                    p_U = min(p_tx_max, p_E_max);
                end
                
                if p_L <= p_U

                    disp('MEC server execution ($\mathcal{P}_{SE}$) is feasible!')

                    virtual_battery = B_hat(t, i);
                    y = @(x) virtual_battery * log2(1 + h*x/sigma) + ...
                        h * (V - virtual_battery*x) / log(2) / (sigma + h*x);
                    p_0 = fsolve(y, 0.5, opt);
                    
                    if (p_U < p_0 && B_hat(t, i) < 0) || B_hat(t, i) >= 0
                        p_mat(i, j) = p_U;
                    elseif p_0 < p_L && B_hat(t, i) < 0
                        p_mat(i, j) = p_L;
                    elseif p_0 >= p_L && p_0 <= p_U && B_hat(t, i) < 0
                        p_mat(i, j) = p_0;
                    end

                    if p_mat(i, j) == 0
                        disp('Something wrong! p is 0!')
                    end
                    

                    server_cost_mat(i, j) = L / (omega * log2(1 + h*p_mat(i, j)/sigma));

                    server_E_mat(i, j) = p_mat(i, j) * server_cost_mat(i, j);

                    J_s(i, j) = (-B_hat(t, i) * p_mat(i, j) + V) * server_cost_mat(i, j);
                    

                else





                    disp('MEC server execution ($\mathcal{P}_{SE}$) is not feasible!')
                    p_mat(i, j) = 0;
                    server_cost_mat(i, j) = 0;
                    server_E_mat(i, j) = 0;

                    J_s(i, j) = inf;
                    




                end
            end
            

            J_d = V * phi;
            disp(['J_m(i):', num2str(J_m(i))])
            disp(['J_s(i,:):', num2str(J_s(i, :))])
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
        for j = 1: M

            device_j_pairs = device_server_pairs(device_server_pairs(:, 2) == j, :);

            is = device_j_pairs(:, 1);
            if isempty(is)
                disp(['For current MEC server #', num2str(j), ', no mobile device choose it!'])

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
                        i = idxs(idx);
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
    

    B(t + 1, :) = B(t, :) - final_chosen_E(t, :) + e(t, :);
    t = t + 1;
    
end



figure
plot(1:T, B(1:T, :))
hold on
plot(1:T, repmat(theta + E_H_max, [T, 1]), '-')
title('Envolution of battery energy level')
xlabel('time slot')
ylabel('battery energy level $B_t$ of each mobile device', 'Interpreter','latex')



average_cost = zeros(T, N);
figure
for i = 1: N

    accumulated = 0;
    request_num = 0;
    for t = 1: T
        accumulated = accumulated + final_chosen_cost(t, i);
        if chosen_mode(t, i) ~= 4

            request_num = request_num + 1;
        end
        average_cost(t, i) = accumulated / request_num;
    end
    plot(1:T, average_cost(:, i));
    hold on
end
title('Envolution of average execution cost')
xlabel('time slot')
ylabel('average execution cost $\frac{1}{T} \sum_{t=0}^{T-1} cost^t$ of each mobile device', 'Interpreter','latex')


average_ratio = zeros(T, 3);
mobile_exe = 0; server_exe = 0; drop = 0;
request_num = 0;
i = 1;
for t = 1: T
    if final_chosen_cost(t, i) == 0
        continue
    else
        request_num = request_num + 1;
        if chosen_mode(t, i) == 1
            mobile_exe = mobile_exe + 1;
        elseif chosen_mode(t, i) == 2
            server_exe = server_exe + 1;
        else
            drop = drop + 1;
        end
    end
    average_ratio(t, :) = [mobile_exe, server_exe, drop] / request_num;
end
figure
plot(1:T, average_ratio(:, 1));
hold on
plot(1:T, average_ratio(:, 2));
hold on
plot(1:T, average_ratio(:, 3));
legend('mobile execution', 'MEC server execution', 'drop')
title('Envolution of average ratio of chosen modes')
xlabel('time slot')
ylabel('average  ratio of chosen modes $\frac{1}{T} \sum_{t=0}^{T-1} \{I_m^t, I_s^t, I_d^t\}$ of the i-th mobile device', 'Interpreter','latex')
