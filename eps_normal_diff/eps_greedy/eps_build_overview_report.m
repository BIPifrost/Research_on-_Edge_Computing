function overview = eps_build_overview_report(output_dir)
% EPS_BUILD_OVERVIEW_REPORT
% Build one dashboard from all eps_*.mat sweep outputs.

if nargin < 1
    output_dir = '.';
end

params = {
    'V', 'eps', 'N', 'M', 'max_connects', 'rho', 'L', ...
    'tau_d', 'E_H_max', 'max_distance', 'min_distance', 'f_max', 'p_tx_max'
};

data = struct('param', {}, 'values', {}, 'overall', {});
for k = 1:numel(params)
    p = params{k};
    mat_path = fullfile(output_dir, ['eps_', p, '.mat']);
    if ~isfile(mat_path)
        warning('Missing sweep file: %s', mat_path);
        continue;
    end
    S = load(mat_path, 'sweep');
    if ~isfield(S, 'sweep')
        warning('Invalid sweep struct in %s', mat_path);
        continue;
    end
    row = struct();
    row.param = p;
    row.values = S.sweep.values(:);
    row.overall = S.sweep.overall;
    data(end + 1) = row; %#ok<AGROW>
end

if isempty(data)
    error('No sweep .mat files found. Run eps_*.m scripts first.');
end

fig = figure('Name', 'eps Overview Report', 'Position', [60 40 1750 980]);
lay = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(lay, 'eps-greedy Single-Variable Sweep Overview');

for k = 1:numel(data)
    nexttile;
    x = data(k).values;
    y = data(k).overall;
    plot(x, y(:, 1), '-o', 'LineWidth', 1.1, 'MarkerSize', 4); hold on;
    plot(x, y(:, 2), '-s', 'LineWidth', 1.1, 'MarkerSize', 4);
    plot(x, y(:, 3), '-^', 'LineWidth', 1.1, 'MarkerSize', 4);
    title(data(k).param, 'Interpreter', 'none');
    xlabel('value');
    ylabel('ratio');
    ylim([0, 1]);
    grid on;
    if k == 1
        legend({'Local', 'Remote', 'Drop'}, 'Location', 'best');
    end
end

for k = (numel(data) + 1):16
    nexttile;
    axis off;
end

saveas(fig, fullfile(output_dir, 'eps_overview_report.png'));
savefig(fig, fullfile(output_dir, 'eps_overview_report.fig'));

overview = struct();
overview.timestamp = char(datetime('now'));
overview.data = data;
save(fullfile(output_dir, 'eps_overview_report.mat'), 'overview');
end
