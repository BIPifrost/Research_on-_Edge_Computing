function sweep = greedy_param_sweep_runner(param_name, values, opts)
% GREEDY_PARAM_SWEEP_RUNNER
% Single-variable sweep for greedy_logco_basic.

if nargin < 3
    opts = struct();
end
if ~isfield(opts, 'T'), opts.T = 1000; end
if ~isfield(opts, 'seed'), opts.seed = 20260307; end
if ~isfield(opts, 'save_plot'), opts.save_plot = true; end
if ~isfield(opts, 'save_mat'), opts.save_mat = false; end
if ~isfield(opts, 'output_dir'), opts.output_dir = '.'; end
if ~isfield(opts, 'extra_cfg'), opts.extra_cfg = struct(); end
if ~isfield(opts, 'title_suffix'), opts.title_suffix = ''; end

if numel(values) < 2
    error('values must contain at least 2 points.');
end

labels = arrayfun(@(x) num2str(x, '%.4g'), values, 'UniformOutput', false);
param_label = format_param_label(param_name);
results = cell(numel(values), 1);
overall = zeros(numel(values), 3);

for k = 1:numel(values)
    cfg = struct();
    cfg.T = opts.T;
    cfg.seed = opts.seed;
    cfg.make_plot = false;
    cfg = apply_overrides(cfg, opts.extra_cfg);
    cfg.(param_name) = values(k);

    if strcmp(param_name, 'L')
        cfg.tau_d = min(max(0.002, 0.002 * (values(k) / 1000)), 0.010);
    end

    results{k} = greedy_logco_basic(cfg);
    overall(k, :) = results{k}.mode_ratio_overall;
end

figure('Name', ['greedy_', param_name], 'Position', [120 120 1380 520]);

subplot(1, 2, 1);
bar(overall, 'grouped');
set(gca, 'XTick', 1:numel(values), 'XTickLabel', labels);
legend({'Local', 'Remote', 'Drop'}, 'Location', 'best');
title(['Overall Mode Ratio vs ', param_label]);
xlabel(param_label);
ylabel('Ratio');
ylim([0, 1]);
grid on;

subplot(1, 2, 2);
for k = 1:numel(values)
    plot(results{k}.mode_ratio_cum(:, 2), 'LineWidth', 1.3); hold on;
end
for k = 1:numel(values)
    plot(results{k}.mode_ratio_cum(:, 3), '--', 'LineWidth', 1.1);
end
legend_items = cell(1, 2 * numel(values));
for k = 1:numel(values)
    legend_items{k} = ['remote, ', param_label, '=', labels{k}];
end
for k = 1:numel(values)
    legend_items{numel(values) + k} = ['drop, ', param_label, '=', labels{k}];
end
legend(legend_items, 'Location', 'best');
title('Cumulative Remote/Drop Ratio');
xlabel('Time Slot');
ylabel('Ratio');
ylim([0, 1]);
grid on;

sgtitle(compose_sgtitle(param_label, opts.title_suffix));

prefix = fullfile(opts.output_dir, ['greedy_', param_name]);
if opts.save_plot
    saveas(gcf, [prefix, '.png']);
end

sweep = struct();
sweep.param_name = param_name;
sweep.values = values;
sweep.labels = labels;
sweep.overall = overall;
sweep.results = results;

if opts.save_mat
    save([prefix, '.mat'], 'sweep');
end
end

function label = format_param_label(param_name)
if strcmp(param_name, 'E_H_max')
    label = 'E_{H}MAX';
elseif strcmp(param_name, 'f_max')
    label = 'F_{MAX}';
elseif strcmp(param_name, 'max_distance')
    label = 'MAX_{distance}';
elseif strcmp(param_name, 'p_tx_max')
    label = 'P_{TX}MAX';
else
    label = upper(param_name);
end
end

function out = apply_overrides(base, extra_cfg)
out = base;
if isempty(extra_cfg)
    return;
end

fields = fieldnames(extra_cfg);
for k = 1:numel(fields)
    out.(fields{k}) = extra_cfg.(fields{k});
end
end

function title_text = compose_sgtitle(param_label, title_suffix)
title_text = ['Single-Variable Sweep: ', param_label];
if ~isempty(title_suffix)
    title_text = sprintf('%s\n%s', title_text, title_suffix);
end
end
