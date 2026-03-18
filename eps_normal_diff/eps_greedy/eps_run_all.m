clc;
close all;

scripts = {
    'eps_V.m', 'eps_eps.m', 'eps_N.m', 'eps_M.m', 'eps_max_connects.m', ...
    'eps_rho.m', 'eps_L.m', 'eps_tau_d.m', 'eps_E_H_max.m', ...
    'eps_max_distance.m', 'eps_min_distance.m', 'eps_f_max.m', 'eps_p_tx_max.m'
};

for k = 1:numel(scripts)
    fprintf('Running %s ...\n', scripts{k});
    run(scripts{k});
end

fprintf('Building overview report ...\n');
eps_build_overview_report('.');
fprintf('All parameter sweeps completed.\n');
