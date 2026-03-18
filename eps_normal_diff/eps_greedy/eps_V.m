clc;
close all;

opts = struct();
opts.T = 800;
opts.seed_env = 20260307;
opts.seed_policy = 20260308;
opts.save_plot = true;
opts.save_mat = true;
opts.output_dir = '.';

values = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3];
sweep = eps_param_sweep_runner('V', values, opts);

disp(['Completed sweep for V.']);
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
