clc;
close all;

opts = struct();
opts.T = 1000;
opts.seed = 20260307;
opts.save_plot = true;
opts.save_mat = false;
opts.output_dir = '.';

values = [0.1, 0.2, 0.4, 0.6, 1.0];
sweep = greedy_param_sweep_runner('p_tx_max', values, opts);

disp('Completed sweep for p_tx_max.');
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
