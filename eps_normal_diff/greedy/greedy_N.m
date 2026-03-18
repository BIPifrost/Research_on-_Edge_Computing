clc;
close all;

opts = struct();
opts.T = 1000;
opts.seed = 20260307;
opts.save_plot = true;
opts.save_mat = false;
opts.output_dir = '.';

values = [10, 30, 50, 80, 100];
sweep = greedy_param_sweep_runner('N', values, opts);

disp('Completed sweep for N.');
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
