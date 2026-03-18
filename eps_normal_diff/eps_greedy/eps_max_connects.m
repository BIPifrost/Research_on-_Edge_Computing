clc;
close all;

opts = struct();
opts.T = 800;
opts.seed_env = 20260307;
opts.seed_policy = 20260308;
opts.save_plot = true;
opts.save_mat = true;
opts.output_dir = '.';

values = [2, 4, 6, 8, 10];
sweep = eps_param_sweep_runner('max_connects', values, opts);

disp(['Completed sweep for max_connects.']);
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
