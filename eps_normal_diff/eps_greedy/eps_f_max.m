clc;
close all;

opts = struct();
opts.T = 800;
opts.seed_env = 20260307;
opts.seed_policy = 20260308;
opts.save_plot = true;
opts.save_mat = true;
opts.output_dir = '.';

values = [1e9, 2e9, 3e9, 4e9, 5e9];
sweep = eps_param_sweep_runner('f_max', values, opts);

disp(['Completed sweep for f_max.']);
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
