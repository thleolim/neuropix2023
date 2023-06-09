%% Create tables of summary statistics %% 
clear; clc; close all
load("C:\Users\taeho\Desktop\output\variables\experiment.mat");
load("C:\Users\taeho\Desktop\output\variables\mean_fr.mat");
load("C:\Users\taeho\Desktop\output\results.mat");

%% Create a cell array storing all significant signed-rank result
orients = {'-60°', '-30°', '0°', '30°', '60°', '90°'};
freqs = {'0.033 c/d', '0.067 c/d', '0.100 c/d','0.133 c/d'};
locs = {'superior', 'central', 'inferior', ...
    'left', 'front', 'right'};
nats = cellstr(strcat('Natural #', num2str((1:30)')));
choice = cellstr(strcat('Choice #', num2str((1:6)')));
stimlabels = {orients, freqs, locs, nats, choice};

sig_wsr = cell(154, 8);      % column 1: unit ID
                             % column 2: unit ROI
                             % column 3: stimulus name
                             % column 4: W test statistic
                             % column 5: n sample size
                             % column 6: p-value
                             % column 7: Hodges-Lehmann estimator
                             % column 8: rank biserial correlation

id = 1;
for u = 1:92
    for s = 1:5
        if ~isnan(wsr_result{u,6}{s})
            if s == 1 || s == 3 || s == 5
                n_items = 6;
            elseif s == 2
                n_items = 4;
            elseif s == 4
                n_items = 30;
            end
            for i = 1:n_items
                if wsr_result{u,6}{s}(i) == 1
                    sig_wsr{id,1} = [num2str(wsr_result{u,1}) '.' num2str(wsr_result{u,2})];
                    if wsr_result{u,1} == 4
                        sig_wsr{id,2} = 'CP';
                    else
                        sig_wsr{id,2} = 'SNr';
                    end
                    sig_wsr{id,3} = string(stimlabels{s}{i});
                    sig_wsr{id,4} = sprintf('%0.1f',wsr_result{u,7}{s}(i,1));
                    sig_wsr{id,5} = wsr_result{u,7}{s}(i,2);
                    sig_wsr{id,6} = sprintf('%0.2e',wsr_result{u,3}{s}(i));
                    sig_wsr{id,7} = sprintf('%0.2f',wsr_result{u,4}{s}(i));
                    sig_wsr{id,8} = sprintf('%0.2f',wsr_result{u,5}{s}(i));
                    id = id + 1;
                end
            end
        end
    end
end

%% Table 1: Sign-rank results for unit-stimulus responsive to one item only
orients = {'Grating: -60°', 'Grating: -30°', 'Grating: 0°', 'Grating: 30°',...
    'Grating: 60°', 'Grating: 90°'};
freqs = {'Grating: 0.033 c/d', 'Grating: 0.067 c/d', 'Grating: 0.100 c/d',...
    'Grating: 0.133 c/d'};
locs = {'Location: Superior', 'Location: Central', 'Location: Inferior', ...
    'Location: Left', 'Location: Front', 'Location: Right'};
nats = cellstr(strcat('Natural Image: ', num2str((1:30)')));
choice = {'Choice: Grating (Front)', 'Choice: Eagle (Front)', 'Choice: Fractal (Front)', ...
    'Choice: Grating (Left)', 'Choice: Eagle (Left)', 'Choice: Fractal (Left)'};
stimlabels = {orients, freqs, locs, nats, choice};

unit_id = cell(size(wsr_forone, 1), 1);
roi = cell(size(wsr_forone, 1), 1);
stim = cell(size(wsr_forone, 1), 1);
wdesc = cell(size(wsr_forone, 1), 1);
ndesc = cell(size(wsr_forone, 1), 1);
pdesc = cell(size(wsr_forone, 1), 1);
rdesc = cell(size(wsr_forone, 1), 1);

for idx = 1:size(wsr_forone, 1)
    e = wsr_forone{idx,1}; u = wsr_forone{idx,2};
    s = wsr_forone{idx,3}; i = wsr_forone{idx,4};

    unit_id{idx} = ['Unit ' num2str(e) '.' num2str(u)];
    roi{idx} = experiment(e).UnitROI(u); roi{idx} = roi{idx}{:};
    stim{idx} = stimlabels{s}{i};

    w = num2str(wsr_forone{idx,8}(1));
    n = num2str(wsr_forone{idx,8}(2));
    p = sprintf('%0.2e', wsr_forone{idx,5});
    r = sprintf('%0.4f', wsr_forone{idx,7});
    wdesc{idx} = ['W(' n ') = ' w];
    pdesc{idx} = ['p = ' p];
    rdesc{idx} = ['r = ' r];
end

wsr_tbl = table(unit_id, roi, stim, wdesc, pdesc, rdesc, 'VariableNames', ...
    {'Unit ID', 'Putative ROI', 'Stimulus', 'W', 'p', 'r'});

clearvars -except experiment mean_fr wsr_forone wsr_tbl kw_result ph_result

%% Table 2: Kruskal-Wallis results for unit-stimulus responsive to 2+ items
stimcats = {'Orientation', 'Frequency', 'Location', 'Natural Image', 'Choice Image'};
orients = {'-60°', '-30°', '0°', '30°', '60°', '90°'};
freqs = {'0.033 c/d', '0.067 c/d', '0.100 c/d', '0.133 c/d'};
locs = {'Superior', 'Central', 'Inferior', 'Left', 'Front', 'Right'};
nats = cellstr(num2str((1:30)'));
% There's no Kruskal-Wallis test for choice images
stimitems = {orients, freqs, locs, nats};

unit_id = cell(size(kw_result, 1), 1);
roi = cell(size(kw_result, 1), 1);
stimc = cell(size(kw_result, 1), 1);
stimi = cell(size(kw_result, 1), 1);
kwchi = cell(size(kw_result, 1), 1);
kwp = cell(size(kw_result, 1), 1);
kweta = cell(size(kw_result, 1), 1);

for idx = 1:size(kw_result, 1)
    e = kw_result{idx,1}; u = kw_result{idx,2};
    s = kw_result{idx,3}; i = kw_result{idx,4};

    unit_id{idx} = [num2str(e) '.' num2str(u)];
    roi{idx} = experiment(e).UnitROI(u); roi{idx} = roi{idx}{:};

    stimc{idx} = stimcats{s};
    thisitems = cell(length(i), 1);
    for id = 1:length(i)
        thisitems{id} = stimitems{s}{i(id)};
    end
    stimi{idx} = [strjoin(thisitems, ', ')];

    df = sprintf('%d', kw_result{idx,6}(1));
    chi = sprintf('%0.3f', kw_result{idx,6}(2));
    if kw_result{idx,7} < 0.0001
        p = '< 0.0001';
    else
        p = sprintf('= %0.4f', kw_result{idx,7});
    end
    if kw_result{idx,8} < 0.01
        eta = '< 0.010';
    else
        eta = sprintf('= %0.3f', kw_result{idx,8});
    end
    kwchi{idx} = ['χ2(' df ') = ' chi];
    kwp{idx} = ['p ' p];
    kweta{idx} = ['η2 ' eta];
end

kw_tbl = table(unit_id, roi, stimc, stimi, kwchi, kwp, kweta, 'VariableNames', ...
    {'Unit ID', 'ROI', 'Stimulus', 'Items', 'Chi', 'p', 'eta'});

clearvars -except experiment mean_fr wsr_forone wsr_tbl kw_result kw_tbl ph_result

%% Table 3: Unit Classification with Parameters
load("C:\Users\taeho\Desktop\output\unitclass.mat")

unit_id = cell(size(unit_char, 1), 1);
meanfr = cell(size(unit_char, 1), 1);
wavedur = cell(size(unit_char, 1), 1);
cv2 = cell(size(unit_char, 1), 1);
prop = cell(size(unit_char, 1), 1);
class = cell(size(unit_char, 1), 1);

for idx = 1:size(unit_char, 1)
    unit_id{idx} = ['Unit ' num2str(unit_char{idx,1}(1)) '.' num2str(unit_char{idx,1}(2))];
    meanfr{idx} = [num2str(round(unit_char{idx,2}(1),2)) ' sp/s ± ' ...
        num2str(round(unit_char{idx,2}(2),3)) ' sp/s'];
    wavedur{idx} = [num2str(round(unit_char{idx,3})) ' μs'];
    cv2{idx} = num2str(round(unit_char{idx,4},3));
    prop{idx} = [num2str(round(unit_char{idx,5},4))];
    class{idx} = unit_char{idx,6};
end

char_tbl = table(unit_id, meanfr, wavedur, cv2, prop, class, 'VariableNames', ...
    {'Unit ID', 'Baseline Mean Firing Rate', 'Waveform Duration', 'CV2', 'Proportion of ISI > 2 s', 'Putative Cell-Type'});

%% Table 4: Post-hoc table
stimcats = {'Orientation', 'Frequency'};
orients = {'-60°', '-30°', '0°', '30°', '60°', '90°'};
freqs = {'0.033 c/d', '0.067 c/d', '0.100 c/d', '0.133 c/d'};
stimitems = {orients, freqs};

unit_id = cell(size(ph_result, 1), 1);
roi = cell(size(ph_result, 1), 1);
stimc = cell(size(ph_result, 1), 1);
stimi = cell(size(ph_result, 1), 1);
phres = cell(size(ph_result, 1), 1);

for idx = 1:size(ph_result, 1)
    e = ph_result{idx,1}; u = ph_result{idx,2};
    s = ph_result{idx,3}; i = reshape(ph_result{idx,4}', 1, []);

    unit_id{idx} = ['Unit ' num2str(e) '.' num2str(u)];
    roi{idx} = experiment(e).UnitROI(u); roi{idx} = roi{idx}{:};

    stimc{idx} = stimcats{s};
    thisitems = cell(numel(i), 1);
    for id = 1:length(i)
        thisitems{id} = stimitems{s}{i(id)};
    end
    for id = 1:2:numel(i)-1
        theseitems = {thisitems{id}, thisitems{id+1}};
        if ph_result{idx,5} < 0
            stimi{idx} = [stimi{idx} ' ' strjoin(theseitems, ' < ')];
        else
            stimi{idx} = [stimi{idx} ' ' strjoin(theseitems, ' > ')];
        end
    end
    for id = 1:length(ph_result{idx,5})
        phres{idx} = [phres{idx} ' ' sprintf('z = %0.1f, p = %0.4f', ...
            ph_result{idx,5}(id), ph_result{idx,6}(id))];
    end
end

ph_tbl = table(unit_id, roi, stimc, stimi, phres, 'VariableNames', ...
    {'Unit ID', 'ROI', 'Stimulus', 'Items', 'Post-hoc'});

clearvars -except experiment mean_fr wsr_forone wsr_tbl kw_result kw_tbl ph_result ph_tbl
