clc;
close all;

opts = struct();
opts.T = 800;
opts.seed_env = 20260307;
opts.seed_policy = 20260308;
opts.save_plot = true;
opts.save_mat = true;
opts.output_dir = '.';

values = [0.2, 0.35, 0.5, 0.7, 0.9];
sweep = eps_param_sweep_runner('rho', values, opts);

disp(['Completed sweep for rho.']);
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
