clc;
close all;

opts = struct();
opts.T = 1000;
opts.seed = 20260307;
opts.save_plot = true;
opts.save_mat = false;
opts.output_dir = '.';

values = [1e-5, 3e-5, 2e-4, 1e-3, 3e-3];
sweep = greedy_param_sweep_runner('E_H_max', values, opts);

disp('Completed sweep for E_H_max.');
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
