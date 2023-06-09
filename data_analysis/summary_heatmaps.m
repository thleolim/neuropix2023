%% Summary Heatmap Figures %%
clear; clc; close all
load("C:\Users\taeho\Desktop\output\results.mat");

%% Summary heatmap of CP units
cp_heatmap = zeros(40, 5); % 0 means non-responsive
cp_units = cell(40, 1);

nsel_units = zeros(0, 1); id = 1;
% Find the non-selective units among kw_result
for idx = 1:size(kw_result, 1)
    if kw_result{idx,7} >= 0.05
        if (kw_result{idx,3} == 1 || kw_result{idx,3} == 3 || ...
                kw_result{idx,3} == 5) && numel(kw_result{idx,4}) == 6
            nsel_units(id) = 1000*kw_result{idx,1} + 10*kw_result{idx,2} ...
                + kw_result{idx,3}; id = id + 1;
        elseif kw_result{idx,3} == 2 && numel(kw_result{idx,4}) == 4 
            nsel_units(id) = 1000*kw_result{idx,1} + 10*kw_result{idx,2} ...
                + kw_result{idx,3}; id = id + 1;
        elseif kw_result{idx,3} == 4 && numel(kw_result{idx,4}) == 30
            nsel_units(id) = 1000*kw_result{idx,1} + 10*kw_result{idx,2} ...
                + kw_result{idx,3}; id = id + 1;
        end
    elseif idx == 2 % failed post-hoc
        nsel_units(id) = 4032;
        id = id + 1;
    end
end

for idx = 1:size(cp_heatmap, 1)
    e = wsr_result{idx,1}; u = wsr_result{idx,2};
    cp_units{idx} = ['Unit ' num2str(e) '.' num2str(u)];
    for s = 1:5
        if ~isnan(wsr_result{idx,6}{s}) 
            if mean(wsr_result{idx,6}{s}) ~= 0
                if mean(wsr_result{idx,6}{s}) ~= 1
                    cp_heatmap(idx,s) = 1; % selective
                else
                    id = 1000*e + 10*u + s;
                    if ismember(id, nsel_units) 
                        cp_heatmap(idx,s) = 2; % non-selective
                    else
                        cp_heatmap(idx,s) = 1; % selective
                    end
                end
            end
        else
            cp_heatmap(idx,s) = NaN; % excluded by compute_timechunks
        end
    end
end

clearvars -except cp_units cp_heatmap nsel_units kw_result wsr_result

%% Summary heatmap of SNr units
snr_heatmap = zeros(52, 5); % 0 means non-responsive
snr_units = cell(52, 1);

for idx = 41:92
    e = wsr_result{idx,1}; u = wsr_result{idx,2};
    snr_units{idx-40} = ['Unit ' num2str(e) '.' num2str(u)];
    for s = 1:5
        if ~isnan(wsr_result{idx,6}{s}) 
            if mean(wsr_result{idx,6}{s}) ~= 0
                if mean(wsr_result{idx,6}{s}) ~= 1
                    snr_heatmap(idx-40,s) = 1; % selective
                else
                    id = 1000*e + 10*u + s;
                    if ismember(id, nsel_units) 
                        snr_heatmap(idx-40,s) = 2; % non-selective
                    else
                        snr_heatmap(idx-40,s) = 1; % selective
                    end
                end
            end
        else
            snr_heatmap(idx-40,s) = NaN; % excluded by compute_timechunks
        end
    end
end

clearvars -except cp_units cp_heatmap snr_units snr_heatmap nsel_units kw_result wsr_result

%% Summary heatmaps of all responsive units to items
% Get the ID and name of all responsive CP and SNr units
r_cpid = zeros(0,1); r_snrid = zeros(0,1);
r_cpnames = cell(0,1); r_snrnames = cell(0,1);
id = 1;
for idx = 1:size(cp_heatmap, 1)
    if any(cp_heatmap(idx, :) == 1 | cp_heatmap(idx, :) == 2)
        r_cpid(id) = idx;
        r_cpnames{id} = cp_units{idx};
        id = id + 1;
    end
end
id = 1;
for idx = 1:size(snr_heatmap, 1)
    if any(snr_heatmap(idx, :) == 1 | snr_heatmap(idx, :) == 2)
        r_snrid(id) = idx+40;
        r_snrnames{id} = snr_units{idx};
        id = id + 1;
    end
end

% Make heatmaps of each responsive CP/SNr units with each stimulus types
heatmap = cell(2, 5); stim_l = [6 4 6 30 6];
for s = 1:5
    heatmap{1,s} = zeros(length(r_cpid), stim_l(s));
    for idx = 1:length(r_cpid)
        id = r_cpid(idx);
        if ~isnan(wsr_result{id,6}{s})
            heatmap{1,s}(idx,:) = wsr_result{id,6}{s};
            for l = 1:stim_l(s)
                if heatmap{1,s}(idx,l) ~= 0
                    heatmap{1,s}(idx,l) = wsr_result{id,5}{s}(l);
                end
            end
        else
            heatmap{1,s}(idx,:) = NaN;
        end
    end
    heatmap{2,s} = zeros(length(r_snrid), stim_l(s));
    for idx = 1:length(r_snrid)
        id = r_snrid(idx);
        if ~isnan(wsr_result{id,6}{s})
            heatmap{2,s}(idx,:) = wsr_result{id,6}{s};
            for l = 1:stim_l(s)
                if heatmap{2,s}(idx,l) ~= 0
                    heatmap{2,s}(idx,l) = wsr_result{id,5}{s}(l);
                end
            end
        else
            heatmap{2,s}(idx,:) = NaN;
        end
    end
end

% Make heatmaps of each responsive CP/SNr units with each stimulus types
heatmap2 = cell(2, 5); stim_l = [6 4 6 30 6];
for s = 1:5
    heatmap2{1,s} = zeros(length(r_cpid), stim_l(s));
    for idx = 1:length(r_cpid)
        id = r_cpid(idx);
        if ~isnan(wsr_result{id,6}{s})
            heatmap2{1,s}(idx,:) = wsr_result{id,6}{s};
            for l = 1:stim_l(s)
                if heatmap2{1,s}(idx,l) ~= 0
                    heatmap2{1,s}(idx,l) = wsr_result{id,4}{s}(l);
                else
                    heatmap2{1,s}(idx,l) = -13;
                end
            end
        else
            heatmap2{1,s}(idx,:) = NaN;
        end
    end
    heatmap2{2,s} = zeros(length(r_snrid), stim_l(s));
    for idx = 1:length(r_snrid)
        id = r_snrid(idx);
        if ~isnan(wsr_result{id,6}{s})
            heatmap2{2,s}(idx,:) = wsr_result{id,6}{s};
            for l = 1:stim_l(s)
                if heatmap2{2,s}(idx,l) ~= 0
                    heatmap2{2,s}(idx,l) = wsr_result{id,4}{s}(l);
                else
                    heatmap2{2,s}(idx,l) = -13;
                end
            end
        else
            heatmap2{2,s}(idx,:) = NaN;
        end
    end
end

clearvars -except cp_units cp_heatmap snr_units snr_heatmap nsel_units kw_result wsr_result r_cpid r_snrid r_cpnames r_snrnames heatmap heatmap2

%% Plot and save the figures
stims = {'Orientation', 'Frequency', 'Location', 'Natural Images', 'Choice Images'};
map = [0 0 0
    0 0 1
    1 0 0];
map2 = turbo; map2(1,:) = [0 0 0];

f1 = figure(1);
set(gcf, "Position", [1 1 500 700*40/52])
imagesc(cp_heatmap, 'AlphaData', ~isnan(cp_heatmap))
yticks(1:size(cp_heatmap,1))
yticklabels(cp_units)
xticks(1:size(cp_heatmap,2))
xticklabels(stims)
caxis([0 2])
colormap(map)
set(gca, "TickLength", [0 0])
saveas(f1, "C:\Users\taeho\Desktop\output\cp_heatmap.tif")

f2 = figure(2);
set(gcf, "Position", [1 1 500 700])
imagesc(snr_heatmap, 'AlphaData', ~isnan(snr_heatmap))
yticks(1:size(snr_heatmap,1))
yticklabels(snr_units)
xticks(1:size(snr_heatmap,2))
xticklabels(stims)
caxis([0 2])
colormap(map)
set(gca, "TickLength", [0 0])
saveas(f2, "C:\Users\taeho\Desktop\output\snr_heatmap.tif")

stim_labels = {'Orientation', 'Frequency (c/d)', 'Grating Location', 'Natural Image ID', 'Choice Image'};
stim_l = [6 4 6 30 6]; orient = {'-60°', '-30°', '0°', '30°', '60°', '90°'};
spafreq = {'0.033', '0.067', '0.100', '0.133'};
loc = {'superior', 'central', 'inferior', 'left', 'front', 'right'};
choice = {'grating (f)', 'eagle (f)', 'fractal (f)', 'grating (l)', 'eagle (l)', 'fractal (l)'};
f3 = figure(3); 
tiledlayout(length(r_cpid)+length(r_snrid), sum(stim_l), 'TileSpacing', 'compact', 'Padding', 'compact')
set(gcf, "Position", [1 1 1000 700*25/52])
for s = 1:5
    nexttile([length(r_cpid) stim_l(s)])
    imagesc(heatmap{1,s}, 'AlphaData', ~isnan(heatmap{1,s}))
    set(gca, "TickLength", [0 0])
    caxis([0 1])
    colormap(map2)
    xticks([])
    if s == 1
        yticks(1:length(r_cpid))
        yticklabels(r_cpnames)
        ylabel('CP Units')
    else
        yticks([])
    end
end
for s = 1:5
    nexttile([length(r_snrid) stim_l(s)])
    imagesc(heatmap{2,s}, 'AlphaData', ~isnan(heatmap{2,s}))
    set(gca, "TickLength", [0 0])
    caxis([0 1])
    colormap(map2)
    xlabel(stim_labels{s})
    xticks(1:stim_l(s))
    if s == 1
        ylabel('SNr Units')
        yticks(1:length(r_snrid))
        yticklabels(r_snrnames)
        xticklabels(orient)
    else
        yticks([])
        if s == 2
            xticklabels(spafreq)
        elseif s == 3
            xticklabels(loc)
        elseif s == 4
            xticklabels(1:30)
        elseif s == 5
            xticklabels(choice)
        end
    end
end
saveas(f3, "C:\Users\taeho\Desktop\output\heatmap.tif")

stim_labels = {'Orientation', 'Frequency (c/d)', 'Grating Location', 'Natural Image ID', 'Choice Image'};
stim_l = [6 4 6 30 6]; orient = {'-60°', '-30°', '0°', '30°', '60°', '90°'};
spafreq = {'0.033', '0.067', '0.100', '0.133'};
loc = {'superior', 'central', 'inferior', 'left', 'front', 'right'};
choice = {'grating (f)', 'eagle (f)', 'fractal (f)', 'grating (l)', 'eagle (l)', 'fractal (l)'};
f4 = figure(4); 
tiledlayout(length(r_cpid)+length(r_snrid), sum(stim_l), 'TileSpacing', 'compact', 'Padding', 'compact')
set(gcf, "Position", [1 1 1000 700*25/52])
for s = 1:5
    nexttile([length(r_cpid) stim_l(s)])
    imagesc(heatmap2{1,s}, 'AlphaData', ~isnan(heatmap2{1,s}))
    set(gca, "TickLength", [0 0])
    caxis([-13 13])
    colormap(map2)
    xticks([])
    if s == 1
        yticks(1:length(r_cpid))
        yticklabels(r_cpnames)
        ylabel('CP Units')
    else
        yticks([])
    end
end
for s = 1:5
    nexttile([length(r_snrid) stim_l(s)])
    imagesc(heatmap2{2,s}, 'AlphaData', ~isnan(heatmap2{2,s}))
    set(gca, "TickLength", [0 0])
    caxis([-13 13])
    colormap(map2)
    xlabel(stim_labels{s})
    xticks(1:stim_l(s))
    if s == 1
        ylabel('SNr Units')
        yticks(1:length(r_snrid))
        yticklabels(r_snrnames)
        xticklabels(orient)
    else
        yticks([])
        if s == 2
            xticklabels(spafreq)
        elseif s == 3
            xticklabels(loc)
        elseif s == 4
            xticklabels(1:30)
        elseif s == 5
            xticklabels(choice)
        end
    end
end
saveas(f4, "C:\Users\taeho\Desktop\output\heatmap2.tif")

clearvars -except cp_units cp_heatmap snr_units snr_heatmap nsel_units kw_result wsr_result