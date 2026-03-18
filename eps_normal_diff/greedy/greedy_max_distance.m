clc;
close all;

opts = struct();
opts.T = 1000;
opts.seed = 20260307;
opts.save_plot = true;
opts.save_mat = false;
opts.output_dir = '.';

values = [30, 50, 80, 120, 150];
sweep = greedy_param_sweep_runner('max_distance', values, opts);

disp('Completed sweep for max_distance.');
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
