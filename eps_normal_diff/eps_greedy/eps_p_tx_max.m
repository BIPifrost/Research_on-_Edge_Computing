clc;
close all;

opts = struct();
opts.T = 800;
opts.seed_env = 20260307;
opts.seed_policy = 20260308;
opts.save_plot = true;
opts.save_mat = true;
opts.output_dir = '.';

values = [0.1, 0.2, 0.4, 0.6, 1.0];
sweep = eps_param_sweep_runner('p_tx_max', values, opts);

disp(['Completed sweep for p_tx_max.']);
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
