clc;
close all;

opts = struct();
opts.T = 800;
opts.seed_env = 20260307;
opts.seed_policy = 20260308;
opts.save_plot = true;
opts.save_mat = true;
opts.output_dir = '.';
opts.extra_cfg = struct('M', 3, 'max_connects', 2, 'rho', 0.8);
opts.title_suffix = '(M=3 C=2 rho=0.8)';

values = 0.1:0.1:0.9;
sweep = eps_param_sweep_runner('eps', values, opts);

disp(['Completed sweep for eps.']);
disp('Overall mode ratio rows correspond to values order:');
disp(values);
disp(sweep.overall);
