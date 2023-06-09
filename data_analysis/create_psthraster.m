%% Create PSTH + raster plots %%
clear; clc; close all
load("C:\Users\taeho\Desktop\output\variables\experiment.mat");
load("C:\Users\taeho\Desktop\output\results.mat");

%% Cycle through each responsive unit-stimulus pairs 
stims = {'orient', 'spafreq', 'loc', 'natimg', 'cwimg'};

% Plot PSTH and raster together
for x = 1:size(wsr_forone, 1)
    e = wsr_forone{x,1}; u = wsr_forone{x,2};
    s = wsr_forone{x,3}; i = wsr_forone{x,4};
    r = experiment(e).UnitROI(u); r = r{:};
    if s ~= 4
        spike_times = experiment(e).Data(s).spike_times;
        spike_templates = experiment(e).Data(s).spike_templates;
        stimOn_t = experiment(e).Data(s).stimOn_t;
        stim_vals = experiment(e).Data(s).stim_vals;
        trange = experiment(e).Data(s).UnitTimeRange;
    
        spike_times = spike_times(spike_templates == u);
        idx = find(stimOn_t > trange(u,1) & stimOn_t < trange(u,2));
        stimOn_t = stimOn_t(idx);
        stim_vals = stim_vals(idx);
    
        stimOff_t = experiment(e).Data(s).stimOff_t;
        stimOff_t = stimOff_t(idx);
        if length(stimOn_t) ~= length(stimOff_t)
            this_l = min([length(stimOn_t) length(stimOff_t)]);
            stimOn_t = stimOn_t(1:this_l);
            stimOff_t = stimOff_t(1:this_l);
        end
        trial_l = median(stimOff_t - stimOn_t);
    
        output = 'C:\Users\taeho\Desktop\output\rasterpsth';
        filename = strjoin(['\exp' num2str(e) 'unit' num2str(u) stims{s} '_' r '.jpg'], '');
        output = fullfile(output, filename);
    
        plot_this(spike_times, stimOn_t, stim_vals, trial_l, s, i, output)
    end
end
for x = 1:size(kw_result, 1)
    e = kw_result{x,1}; u = kw_result{x,2};
    s = kw_result{x,3}; i = kw_result{x,4};
    r = experiment(e).UnitROI(u); r = r{:};
    if s ~= 4
        spike_times = experiment(e).Data(s).spike_times;
        spike_templates = experiment(e).Data(s).spike_templates;
        stimOn_t = experiment(e).Data(s).stimOn_t;
        stim_vals = experiment(e).Data(s).stim_vals;
        trange = experiment(e).Data(s).UnitTimeRange;
    
        spike_times = spike_times(spike_templates == u);
        idx = find(stimOn_t > trange(u,1) & stimOn_t < trange(u,2));
        stimOn_t = stimOn_t(idx);
        stim_vals = stim_vals(idx);
    
        stimOff_t = experiment(e).Data(s).stimOff_t;
        stimOff_t = stimOff_t(idx);
        if length(stimOn_t) ~= length(stimOff_t)
            this_l = min([length(stimOn_t) length(stimOff_t)]);
            stimOn_t = stimOn_t(1:this_l);
            stimOff_t = stimOff_t(1:this_l);
        end
        trial_l = median(stimOff_t - stimOn_t);
    
        output = 'C:\Users\taeho\Desktop\output\rasterpsth';
        filename = strjoin(['\exp' num2str(e) 'unit' num2str(u) stims{s} '_' r '.jpg'], '');
        output = fullfile(output, filename);
    
        plot_this(spike_times, stimOn_t, stim_vals, trial_l, s, i, output)
    end
end

clearvars -except experiment kw_result ph_result wsr_forone wsr_result

%% Plot PSTH of natural image responses as heatmaps
for x = 1:size(kw_result, 1)
    e = kw_result{x,1}; u = kw_result{x,2};
    s = kw_result{x,3}; i = kw_result{x,4};
    r = experiment(e).UnitROI(u); r = r{:};
    if s == 4
        spike_times = experiment(e).Data(s).spike_times;
        spike_templates = experiment(e).Data(s).spike_templates;
        stimOn_t = experiment(e).Data(s).stimOn_t;
        stim_vals = experiment(e).Data(s).stim_vals;
        trange = experiment(e).Data(s).UnitTimeRange;
    
        spike_times = spike_times(spike_templates == u);
        idx = find(stimOn_t > trange(u,1) & stimOn_t < trange(u,2));
        stimOn_t = stimOn_t(idx);
        stim_vals = stim_vals(idx);
    
        stimOff_t = experiment(e).Data(s).stimOff_t;
        stimOff_t = stimOff_t(idx);
        if length(stimOn_t) ~= length(stimOff_t)
            this_l = min([length(stimOn_t) length(stimOff_t)]);
            stimOn_t = stimOn_t(1:this_l);
            stimOff_t = stimOff_t(1:this_l);
        end
        trial_l = median(stimOff_t - stimOn_t);
    
        output = 'C:\Users\taeho\Desktop\output\rasterpsth';
        filename1 = strjoin(['\exp' num2str(e) 'unit' num2str(u) 'natimg_' r '.tif'], '');
        filename2 = strjoin(['\cbar1_exp' num2str(e) 'unit' num2str(u) 'natimg_' r '.tif'], '');
        filename3 = strjoin(['\cbar2_exp' num2str(e) 'unit' num2str(u) 'natimg_' r '.tif'], '');

        output1 = fullfile(output, filename1);
        output2 = fullfile(output, filename2);
        output3 = fullfile(output, filename3);
    
        get_heatmap(spike_times, stimOn_t, stim_vals, trial_l, i, output1, output2, output3)
    end
end

clearvars -except experiment kw_result ph_result wsr_forone wsr_result

%% Function to generate and save PSTH-raster plot
function plot_this(spike_times, stimOn_t, stim_vals, trial_l, stimType, items, output)

orients = {'-60°', '-30°', '0°', '30°', '60°', '90°'};
freqs = {'0.033 c/d', '0.067 c/d', '0.100 c/d', '0.133 c/d'};
loc = {'superior', 'central', 'inferior', 'left', 'front', 'right'};
cwimg = {'grating (f)', 'eagle (f)', 'fractal (f)', 'grating (l)', 'eagle (l)', 'fractal (l)'};
stims = unique(stim_vals);
colors = {[0 0.4470 0.7410]	[0.8500 0.3250 0.0980] [0.9290 0.6940 0.1250] ...
    [0.4940 0.1840 0.5560] [0.3010 0.7450 0.9330] [0.6350 0.0780 0.1840]};

% Initialize raster
raster_window = [-0.5, 1]; binsz = 0.001;
t_bin_edge = raster_window(1):binsz:raster_window(2);
t_offset = raster_window(1) - binsz;

% Initialize PSTH
psth_window = [-0.7 1.2]; 
t_edges = psth_window(1):binsz:psth_window(2);
t_centers = t_edges(1:end-1) + diff(t_edges);
actual_t = t_centers(201:1699);
% Generate the Gaussian filter
smooth_sz = 0.051;
gw = gausswin(round(smooth_sz/binsz), 3)';
gw = gw./sum(gw); % normalize amplitude

% Generate the rasters
raster = cell(length(stims), 2);
for i = 1:length(stims)
    this_stimOn_t = stimOn_t(stim_vals == stims(i));
    tempraster = zeros(length(this_stimOn_t), length(t_bin_edge)-1);
    for trial = 1:length(this_stimOn_t)
        tempraster(trial,:) = histcounts(spike_times, ...
            t_bin_edge + this_stimOn_t(trial));
    end
    [raster{i,2}, raster{i,1}] = find(tempraster);
    raster{i,1} = raster{i,1}*binsz + t_offset;
end

% Generate the PSTHs
psth = cell(length(stims), 2);
for i = 1:length(stims)
    this_stimOn_t = stimOn_t(stim_vals == stims(i));
    tempraster2 = zeros(length(this_stimOn_t), length(t_centers));
    for trial = 1:length(this_stimOn_t)
        tempraster2(trial,:) = histcounts(spike_times, ...
            t_edges + this_stimOn_t(trial));
    end
    % Smooth the raster with the Gaussian filter
    smooth_raster = conv2(gw, 1, tempraster2', 'same');
    smooth_raster = smooth_raster'./binsz; % convert to FR
    % Get the PSTH average and standard error of the raster
    psth_avg = mean(smooth_raster);
    psth_sem = std(smooth_raster)./sqrt(sum(smooth_raster));
    fixthese = isnan(psth_sem); psth_sem(fixthese) = 0;
    % Trim the ends of the PSTH to account for convolution
    psth{i,1} = psth_avg(201:1699);
    psth{i,2} = psth_sem(201:1699);
end

% Plot rasters and PSTHs in tiledlayout
f = figure('Position', [1 1 500 700]); 
tiledlayout(length(stims)*3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:length(stims)
    nexttile([2 1])
    if ismember(i, items)
        scatter(raster{i,1}, raster{i,2}, 2, colors{i}, 'filled')
    else
        scatter(raster{i,1}, raster{i,2}, 2, [.5 .5 .5], 'filled')
    end
    xlim(raster_window)
    xline(0, 'g--', 'LineWidth', 1.5)
    xline(trial_l, 'r--', 'LineWidth', 1.5)
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    if stimType == 1
        ylabel(orients{i})
    elseif stimType == 2
        ylabel(freqs{i})
    elseif stimType == 3
        ylabel(loc{i})
    elseif stimType == 5
        ylabel(cwimg{i})
    end
end
% Plot PSTHs
nexttile([length(stims) 1])
labels = cell(length(stims), 1);
for i = 1:length(stims)
    if ismember(i, items)
        fill([actual_t, fliplr(actual_t)], [psth{i,1}+psth{i,2}, ...
            fliplr(psth{i,1}-psth{i,2})], colors{i}, 'FaceAlpha', .5, ...
            'EdgeColor', 'none'); hold on
        plot(actual_t, psth{i,1}, 'Color', colors{i}); hold on
        labels{2*i-1} = '';
        if stimType == 1
            labels{2*i} = orients{i};
        elseif stimType == 2
            labels{2*i} = freqs{i};
        elseif stimType == 3
            labels{2*i} = loc{i};
        elseif stimType == 5
            labels{2*i} = cwimg{i};
        end
    else
        fill([actual_t, fliplr(actual_t)], [psth{i,1}+psth{i,2}, ...
            fliplr(psth{i,1}-psth{i,2})], [.5 .5 .5], 'FaceAlpha', .3, ...
            'EdgeColor', 'none'); hold on
        plot(actual_t, psth{i,1}, 'Color', [.5 .5 .5]); hold on
        labels{2*i-1} = ''; labels{2*i} = '';
    end
end

xline(0, 'g--', 'LineWidth', 1.5)
xline(trial_l, 'r--', 'LineWidth', 1.5)
xlim([actual_t(1) actual_t(end)])
xlabel('time from stimulus onset (s)')
ylabel('firing rate (sp/s)')
legend(labels)

% Save figure in output
saveas(f, output)
close(f)

end

function get_heatmap(spike_times, stimOn_t, stim_vals, trial_l, items, output1, output2, output3)

stims = unique(stim_vals);

% Initialize PSTH
psth_window = [-0.7 1.2]; binsz = 0.001;
t_edges = psth_window(1):binsz:psth_window(2);
t_centers = t_edges(1:end-1) + diff(t_edges);
actual_t = round(t_centers(201:1699),4);
% Generate the Gaussian filter
smooth_sz = 0.051;
gw = gausswin(round(smooth_sz/binsz), 3)';
gw = gw./sum(gw); % normalize amplitude

% Generate the PSTHs
psth = cell(length(stims), 1);
for i = 1:length(stims)
    this_stimOn_t = stimOn_t(stim_vals == stims(i));
    tempraster2 = zeros(length(this_stimOn_t), length(t_centers));
    for trial = 1:length(this_stimOn_t)
        tempraster2(trial,:) = histcounts(spike_times, ...
            t_edges + this_stimOn_t(trial));
    end
    % Smooth the raster with the Gaussian filter
    smooth_raster = conv2(gw, 1, tempraster2', 'same');
    smooth_raster = smooth_raster'./binsz; % convert to FR
    % Get the PSTH average and standard error of the raster
    psth_avg = mean(smooth_raster);
    % Trim the ends of the PSTH to account for convolution
    psth{i} = psth_avg(201:1699);
end

% Figure out caxis
thismax = 0;
for i = 1:size(psth,1)
    if i == 1
        thismin = floor(min(psth{i}));
    else
        thismin = min([thismin floor(min(psth{i}))]);
    end
    thismax = max([thismax ceil(max(psth{i}))]);
end
[~,thiszero] = find(actual_t == 0);
[~,thisend] = min(abs(actual_t-trial_l));

% 
f = figure("Position", [1 1 1250 700]); 
tlo = tiledlayout(30, 46, 'TileSpacing', 'compact');
for i = 1:size(psth,1)
    nexttile([1 45])
    psth{i}(thiszero) = thismax;
    psth{i}(thisend) = thismax;
    if sum(ismember(items, i)) == 1
        imagesc(psth{i}); colormap(gca, 'turbo'); 
        ylabel(num2str(i))
    else
        imagesc(psth{i}); colormap(gca, 'gray'); 
    end
    set(gca, 'YTick', []);
    if i ~= 30
        set(gca,'XTick',[]);
    else
        set(gca,'XTick',100:100:1400)
        set(gca,'TickLength',[0 0])
        set(gca,'XTickLabel',actual_t(100:100:1400));
    end
end
caxis([thismin thismax])
xlabel(tlo, 'time from stimulus onset (s)')
ylabel(tlo, 'natural image ID')

% Quickly create this colorbar
a = thismin:(thismax-thismin)/(length(actual_t)-1):thismax;
f2 = figure("Position", [1 1 1250 350]);
imagesc(a); colormap('turbo'); colorbar;
f3 = figure("Position", [1 1 1250 350]);
imagesc(a); colormap('gray'); colorbar;

% Save figures in output
saveas(f, output1)
close(f)
saveas(f2, output2)
close(f2)
saveas(f3, output3)
close(f3)

end