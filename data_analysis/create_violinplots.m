%% Create violin plots %%
clear; clc; close all
load("C:\Users\taeho\Desktop\output\variables\experiment.mat");
load("C:\Users\taeho\Desktop\output\variables\mean_fr.mat");
load("C:\Users\taeho\Desktop\output\results.mat");

%% Cycle through all significant units
stims = {'orient', 'spafreq', 'loc', 'natimg', 'cwimg'};
for idx = 1:size(wsr_forone,1)
    s = wsr_forone{idx,3};
    if s ~= 4
        e = wsr_forone{idx,1}; u = wsr_forone{idx,2};
        i = wsr_forone{idx,4}; r = experiment(e).UnitROI(u); r = r{:};
    
        output = 'C:\Users\taeho\Desktop\output\violinplots';
        filename = strjoin(['\exp' num2str(e) 'unit' num2str(u) stims{s} '_' r '.jpg'], '');
        output = fullfile(output, filename);

        plot_this(mean_fr, e, u, s, i, output);
    end
end

for idx = 1:size(kw_result,1)
    s = kw_result{idx,3};
    if s ~= 4
        e = kw_result{idx,1}; u = kw_result{idx,2};
        i = kw_result{idx,4}; r = experiment(e).UnitROI(u); r = r{:};
    
        output = 'C:\Users\taeho\Desktop\output\violinplots';
        filename = strjoin(['\exp' num2str(e) 'unit' num2str(u) stims{s} '_' r '.jpg'], '');
        output = fullfile(output, filename);

        plot_this(mean_fr, e, u, s, i, output);
    end
end

%% Create histogram plots to complement
stims = {'orient', 'spafreq', 'loc', 'natimg', 'cwimg'};
colors = {[0 0.4470 0.7410] [0.8500 0.3250 0.0980] [0.9290 0.6940 0.1250] ...
    [0.4940 0.1840 0.5560] [0.3010 0.7450 0.9330] [0.6350 0.0780 0.1840]};

for idx = 1:size(wsr_forone,1)
    s = wsr_forone{idx,3};
    if s ~= 4
        e = wsr_forone{idx,1}; u = wsr_forone{idx,2};
        i = wsr_forone{idx,4}; r = experiment(e).UnitROI(u); r = r{:};
    
        for id = 1:numel(i)
            output = 'C:\Users\taeho\Desktop\output\histogram';
            filename = strjoin(['\exp' num2str(e) 'unit' num2str(u) stims{s} '_' i(id) '_' r '.jpg'], '');
            output = fullfile(output, filename);

            f = figure('Position', [1 1 400 200]);
            histogram(mean_fr{e,s}{u,i(id)}(:,1), 'BinWidth', 5, 'FaceColor', [.5 .5 .5]); hold on
            histogram(mean_fr{e,s}{u,i(id)}(:,2), 'BinWidth', 5, 'FaceColor', colors{i(id)});
            xlabel('mean firing rate (sp/s)')
            ylabel('count')
            saveas(f, output);
            close(f)
        end
    end
end

for idx = 1:size(kw_result,1)
    s = kw_result{idx,3};
    if s ~= 4
        e = kw_result{idx,1}; u = kw_result{idx,2};
        i = kw_result{idx,4}; r = experiment(e).UnitROI(u); r = r{:};

        for id = 1:numel(i)
            output = 'C:\Users\taeho\Desktop\output\histogram';
            filename = strjoin(['\exp' num2str(e) 'unit' num2str(u) stims{s} '_' i(id) '_' r '.jpg'], '');
            output = fullfile(output, filename);

            f = figure('Position', [1 1 400 200]);
            histogram(mean_fr{e,s}{u,i(id)}(:,1), 'BinWidth', 5, 'FaceColor', [.5 .5 .5]); hold on
            histogram(mean_fr{e,s}{u,i(id)}(:,2), 'BinWidth', 5, 'FaceColor', colors{i(id)});
            xlabel('mean firing rate (sp/s)')
            ylabel('count')
            saveas(f, output);
            close(f)
        end
    end
end

%% Function
function plot_this(mean_fr, e, u, s, i, output)

colors = {[0 0.4470 0.7410] [0.8500 0.3250 0.0980] [0.9290 0.6940 0.1250] ...
    [0.4940 0.1840 0.5560] [0.3010 0.7450 0.9330] [0.6350 0.0780 0.1840]};
thiscolor = {[.5 .5 .5], [.5 .5 .5], [.5 .5 .5], [.5 .5 .5], [.5 .5 .5], ...
    [.5 .5 .5]};
for idx = 1:numel(i)
    thiscolor{i(idx)} = colors{i(idx)};
end

if s ~= 2
    n = max([size(mean_fr{e,s}{u,1},1), ...
        size(mean_fr{e,s}{u,2},1), ...
        size(mean_fr{e,s}{u,3},1), ...
        size(mean_fr{e,s}{u,4},1), ...
        size(mean_fr{e,s}{u,5},1), ...
        size(mean_fr{e,s}{u,6},1)]);
    data = nan(n, 6);
    for idx = 1:6
        data(1:numel(diff(mean_fr{e,s}{u,idx},1,2)),idx) = ...
            mean_fr{e,s}{u,idx}(:,2) - mean_fr{e,s}{u,idx}(:,1);
    end
else
    n = max([size(mean_fr{e,s}{u,1},1), ...
        size(mean_fr{e,s}{u,2},1), ...
        size(mean_fr{e,s}{u,3},1), ...
        size(mean_fr{e,s}{u,4},1)]);
    data = nan(n, 4);
    for idx = 1:4
        data(1:size(mean_fr{e,s}{u,idx},1),idx) = ...
            mean_fr{e,s}{u,idx}(:,2) - mean_fr{e,s}{u,idx}(:,1);
    end
end

if s == 1
    names = {'-60°', '-30°', '0°', '30°', '60°', '90°'};
elseif s == 2
    names = {'0.033 c/d', '0.067 c/d', '0.100 c/d', '0.133 c/d'};
elseif s == 3
    names = {'superior', 'central', 'inferior', 'left', 'front', 'right'};
elseif s == 5
    names = {'grating (f)', 'eagle (f)', 'fractal (f)', 'grating (l)', 'eagle (l)', 'fractal (l)'};
end

f = figure('Position', [1 1 400 600]);
distributionPlot(data, ...
    'color', thiscolor, 'histOpt', 2, 'showMM', 0, ...
    'xNames', names, ...
    'ylabel', 'pairwise difference in mean firing rate (sp/s)');
yline(0, 'r--')

if s == 2
    for idx = 1:4
        diffs = sort(mean_fr{e,s}{u,idx}(:,2) - mean_fr{e,s}{u,idx}(:,1));
        diffs = diffs(diffs ~= 0);
        [I,J] = ndgrid(diffs,diffs); 
        d = triu(I+J)./2; % Walsh averages triang matrix
        ld = sort(d(d~=0)); % Walsh averages vector
        hle = median(ld); % Hodge-Lehmann estimate
        if ismember(idx, i)
            line([idx-.5 idx+.5], [hle hle], 'Color', 'r', 'LineWidth', 2)
        else
            line([idx-.5 idx+.5], [hle hle], 'Color', 'k', 'LineWidth', 2)
        end
    end
else
    for idx = 1:6
        diffs = sort(mean_fr{e,s}{u,idx}(:,2) - mean_fr{e,s}{u,idx}(:,1));
        diffs = diffs(diffs ~= 0);
        [I,J] = ndgrid(diffs,diffs); 
        d = triu(I+J)./2; % Walsh averages triang matrix
        ld = sort(d(d~=0)); % Walsh averages vector
        hle = median(ld); % Hodge-Lehmann estimate
        if ismember(idx, i)
            line([idx-.5 idx+.5], [hle hle], 'Color', 'r', 'LineWidth', 2)
        else
            line([idx-.5 idx+.5], [hle hle], 'Color', 'k', 'LineWidth', 2)
        end
    end
end

saveas(f, output);
close(f)

end