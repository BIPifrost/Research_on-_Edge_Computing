clc;
close all;

opts = struct();
opts.T = 1000;
opts.seed = 20260307;
opts.save_plot = true;
opts.save_mat = false;
opts.output_dir = '.';

values = [1e9, 2e9, 3e9, 4e9, 5e9];
sweep = greedy_param_sweep_runner('f_max', values, opts);

disp('Completed sweep for f_max.');
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
