%% Neuronal Classification %%
clear; clc; close all
load("C:\Users\taeho\Desktop\output\variables\experiment.mat");
load("C:\Users\taeho\Desktop\output\variables\mean_fr.mat");
load("C:\Users\taeho\Desktop\output\results.mat");

%% Get ID of all responsive units
unit_id = zeros(0, 3); % experiment ID, unit ID, unique ID
id = 1;
for idx = 1:size(wsr_forone,1)
    unit_id(id,1) = wsr_forone{idx,1};
    unit_id(id,2) = wsr_forone{idx,2};
    unit_id(id,3) = 100*wsr_forone{idx,1} + wsr_forone{idx,2};
    id = id + 1;
end
for idx = 1:size(kw_result,1)
    unit_id(id,1) = kw_result{idx,1};
    unit_id(id,2) = kw_result{idx,2};
    unit_id(id,3) = 100*kw_result{idx,1} + kw_result{idx,2};
    id = id + 1;
end

[~, this] = unique(unit_id(:,3));
unit_id = unit_id(this, 1:2);
unit_char = cell(size(unit_id,1), 6); % column 1: unit ID [e u]
                                      % column 2: baseline firing rate [avg sem]
                                      % column 3: waveform duration us
                                      % column 4: cv2
                                      % column 5: prop ISI > 2s
                                      % column 6: putative classification
for idx = 1:size(unit_id,1)
    unit_char{idx,1} = unit_id(idx,:);
end

% Now for all BG units
unit_id = zeros(0, 3); % experiment ID, unit ID, unique ID
id = 1;
for idx = 1:size(wsr_result,1)
    unit_id(id,1) = wsr_result{idx,1};
    unit_id(id,2) = wsr_result{idx,2};
    unit_id(id,3) = 100*wsr_result{idx,1} + wsr_result{idx,2};
    id = id + 1;
end

[~, this] = unique(unit_id(:,3));
unit_id = unit_id(this, 1:2);
all_char = cell(size(unit_id,1), 6); % column 1: unit ID [e u]
                                      % column 2: baseline firing rate [avg sem]
                                      % column 3: waveform duration us
                                      % column 4: cv2
                                      % column 5: prop ISI > 2s
                                      % column 6: putative classification
for idx = 1:size(unit_id,1)
    all_char{idx,1} = unit_id(idx,:);
end

clearvars -except experiment mean_fr unit_char all_char

%% Get baseline firing rate for each unit
for idx = 1:size(unit_char,1)
    e = unit_char{idx,1}(1); u = unit_char{idx,1}(2);
    base = [];
    for s = 1:5
        if ~isnan(mean_fr{e,s}{u,1})
            n = size(mean_fr{e,s},2);
            for i = 1:n
                base = [base; mean_fr{e,s}{u,i}(:,2)];
            end
        end
    end
    base = base(~isnan(base));
    unit_char{idx,2} = [mean(base) std(base)/sqrt(length(base))];
end

% Now for all BG units
for idx = 1:size(all_char,1)
    e = all_char{idx,1}(1); u = all_char{idx,1}(2);
    base = [];
    for s = 1:5
        if ~isnan(mean_fr{e,s}{u,1})
            n = size(mean_fr{e,s},2);
            for i = 1:n
                base = [base; mean_fr{e,s}{u,i}(:,2)];
            end
        end
    end
    base = base(~isnan(base));
    all_char{idx,2} = [mean(base) std(base)/sqrt(length(base))];
end

clearvars -except experiment mean_fr unit_char all_char

%% Get template duration for each unit
t = 1e3 * ((0:81) / 30000);

for idx = 1:size(unit_char,1)
    e = unit_char{idx,1}(1); u = unit_char{idx,1}(2);
    templates = squeeze(experiment(e).Templates(u,:,:));
    [~, max_site] = max(max(abs(templates), [], 1));
    waveform = templates(:,max_site);
    [~, trough] = min(waveform);
    waveform_temp = waveform; waveform_temp(1:trough) = 0;
    [~, peak] = max(waveform_temp);
    unit_char{idx,3} = (peak - trough)/30000 * 1e6;

    f = figure;
    plot(t, waveform, 'k');
    xlabel('time (ms)')
    output = fullfile('C:\Users\taeho\Desktop\output\waveforms', ...
        ['exp' num2str(e) '_unit' num2str(u) '.tif']);
    saveas(f, output);
    close(f)
end

% Now for all units
for idx = 1:size(all_char,1)
    e = all_char{idx,1}(1); u = all_char{idx,1}(2);
    templates = squeeze(experiment(e).Templates(u,:,:));
    [~, max_site] = max(max(abs(templates), [], 1));
    waveform = templates(:,max_site);
    [~, trough] = min(waveform);
    waveform_temp = waveform; waveform_temp(1:trough) = 0;
    [~, peak] = max(waveform_temp);
    all_char{idx,3} = (peak - trough)/30000 * 1e6;
end

clearvars -except experiment mean_fr unit_char all_char

%% Get CV2 and ISI proportion of each unit
for idx = 1:size(unit_char,1)
    e = unit_char{idx,1}(1); u = unit_char{idx,1}(2);
    this_trange = experiment(e).Data(1).UnitTimeRange(u,:);
    this_spikes = experiment(e).Data(1).spike_templates == u;
    this_spikes = experiment(e).Data(1).spike_times(this_spikes);
    this_spikes = this_spikes(this_spikes >= this_trange(1) ...
        & this_spikes <= this_trange(2));
    isi = diff(this_spikes);
    
    isi_ratio = zeros(length(isi)-1, 1);
    for i = 1:length(isi_ratio)
        isi_ratio(i) = abs(isi(i)-isi(i+1))/mean([isi(i) isi(i+1)]);
    end
    cv2 = nanmean(isi_ratio);

    prop_isi = sum(isi(isi>2))/diff(this_trange);

    unit_char{idx,4} = cv2;
    unit_char{idx,5} = prop_isi;
end

% Now for all BG units
for idx = 1:size(all_char,1)
    e = all_char{idx,1}(1); u = all_char{idx,1}(2);
    this_trange = experiment(e).Data(1).UnitTimeRange(u,:);
    this_spikes = experiment(e).Data(1).spike_templates == u;
    this_spikes = experiment(e).Data(1).spike_times(this_spikes);
    this_spikes = this_spikes(this_spikes >= this_trange(1) ...
        & this_spikes <= this_trange(2));
    isi = diff(this_spikes);
    
    isi_ratio = zeros(length(isi)-1, 1);
    for i = 1:length(isi_ratio)
        isi_ratio(i) = abs(isi(i)-isi(i+1))/mean([isi(i) isi(i+1)]);
    end
    cv2 = nanmean(isi_ratio);

    prop_isi = sum(isi(isi>2))/diff(this_trange);

    all_char{idx,4} = cv2;
    all_char{idx,5} = prop_isi;
end

clearvars -except experiment mean_fr unit_char all_char

%% Putative classification
msn_n = 0; tan_n = 0; fsi_n = 0; uin_n = 0; cp_n = 0; snr_n = 0;
msn_n2 = 0; tan_n2 = 0; fsi_n2 = 0; uin_n2 = 0; cp_n2 = 0; snr_n2 = 0;

for idx = 1:size(unit_char,1)
    e = unit_char{idx,1}(1); u = unit_char{idx,1}(2);
    if strcmp(experiment(e).UnitROI{u}, 'CP')
        cp_n = cp_n + 1;
        if unit_char{idx,3} > 400
            if unit_char{idx,5} >= 0.35
                unit_char{idx,6} = 'MSN';
                msn_n = msn_n + 1;
            else
                unit_char{idx,6} = 'TAN';
                tan_n = tan_n + 1;
            end
        else
            if unit_char{idx,5} < 0.35
                unit_char{idx,6} = 'FSI';
                fsi_n = fsi_n + 1;
            else
                unit_char{idx,6} = 'UIN';
                uin_n = uin_n + 1;
            end
        end
    else
        unit_char{idx,6} = 'SNr';
        snr_n = snr_n + 1;
    end
end

% Now for all BG units
for idx = 1:size(all_char,1)
    e = all_char{idx,1}(1); u = all_char{idx,1}(2);
    if strcmp(experiment(e).UnitROI{u}, 'CP')
        cp_n2 = cp_n2 + 1;
        if all_char{idx,3} > 400
            if all_char{idx,5} >= 0.35
                all_char{idx,6} = 'MSN';
                msn_n2 = msn_n2 + 1;
            else
                all_char{idx,6} = 'TAN';
                tan_n2 = tan_n2 + 1;
            end
        else
            if all_char{idx,5} < 0.35
                all_char{idx,6} = 'FSI';
                fsi_n2 = fsi_n2 + 1;
            else
                all_char{idx,6} = 'UIN';
                uin_n2 = uin_n2 + 1;
            end
        end
    else
        all_char{idx,6} = 'SNr';
        snr_n2 = snr_n2 + 1;
    end
end

disp([num2str(msn_n) '/' num2str(msn_n2) ' MSN units responsive.'])
disp([num2str(tan_n) '/' num2str(tan_n2) ' TAN units responsive.'])
disp([num2str(fsi_n) '/' num2str(fsi_n2) ' FSI units responsive.'])
disp([num2str(uin_n) '/' num2str(uin_n2) ' UIN units responsive.'])
disp([num2str(cp_n) '/' num2str(cp_n2) ' CP units responsive.'])
disp([num2str(snr_n) '/' num2str(snr_n2) ' SNr units responsive.'])

clearvars -except unit_char all_char

%% Save
save("C:\Users\taeho\Desktop\output\unitclass.mat")