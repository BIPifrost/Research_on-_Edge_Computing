clc;
close all;

opts = struct();
opts.T = 800;
opts.seed_env = 20260307;
opts.seed_policy = 20260308;
opts.save_plot = true;
opts.save_mat = true;
opts.output_dir = '.';

values = [0.002, 0.004, 0.006, 0.008, 0.010];
sweep = eps_param_sweep_runner('tau_d', values, opts);

disp(['Completed sweep for tau_d.']);
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
