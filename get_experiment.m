%% Get all my experiment data
function experiment = get_experiment

% Below are all experiments with successful pre-processing
exp = {'JF088_2022-11-29_site1', 'JF088_2022-11-30_site1', ...
    'JF089_2022-11-15_site1', 'JF089_2022-11-15_site2', ...
    'JF090_2022-11-15_site1', 'JF090_2022-11-16_site1'};
% Initialize experiment as struct to store each experiment data
experiment = struct('experiment', exp); % keep it simple
[~, ~, st, ~] = bd_loadAllenAtlas(fullfile('C:\Users\taeho\Desktop\neuropix2023', 'files', 'allen_mouse_25um_v1.2/'));

for e = 1:6
    % Decompose experiment to figure out defintions
    % experiment must be a string in the format e.g. 'JF088_2022-11-29_site1'
    animal = exp{e}(1:5); day = exp{e}(7:16); site = exp{e}(18:22);
    
    % Define the bombcell and kilosort folders
    bcpath = fullfile('C:\Users\taeho\Desktop\output\bombcell_output', [animal '_' day '_' site]);
    kspath = fullfile('D:\recordings', [animal, '_', day, '_', site]);
    load(fullfile(kspath, 'sync.mat'));
    % Extract the good units from the kilosort output folder
    [spike_times, spike_templates, template_depths, ~, ~, ~, template_xdepths, ...
        ~, templates, trange] = load_spikes(bcpath, kspath, site, 0);
    % Store an array of good unit indices in the experiment structf
    experiment(e).UnitID = unique(spike_templates);
    % Get all templates
    experiment(e).Templates = templates;

    % Load the probe_ccf.mat file
    load(fullfile('C:\Users\taeho\Desktop\output\braindraw_output\', animal, day, 'probe_ccf.mat'));
    experiment(e).UnitROI = cell(numel(experiment(e).UnitID), 1);
    % Loop across all units
    for unit = 1:numel(experiment(e).UnitID)
        % Find the shank at which the unit is located at
        if mod(template_xdepths(unit), 250) == 0
            shank = template_xdepths(unit)/250 + 1;
        elseif mod(template_xdepths(unit), 250) == 32
            shank = (template_xdepths(unit)-32)/250 + 1;
        end
        % Normalize the unit depth to probe trajectory coordinates
        if e == 3
            depth = template_depths(unit) - 1440;
        elseif e == 4
            depth = template_depths(unit) - 1440;
        else
            depth = template_depths(unit) - 2160;
        end
        depth = (probe_ccf(shank).probe_depths(1) + depth)./25;
        % Find the ROI ID that corresponds to the unit
        probe_trajectory_depths = pdist2(probe_ccf(shank).trajectory_coords, ...
            probe_ccf(shank).trajectory_coords(1, :));
        coord = abs(probe_trajectory_depths - depth);
        coord = coord - min(coord);
        coord = find(~coord, 1);
        roiID = probe_ccf(shank).trajectory_areas(coord);
        % Store the ROI name to the structure finally
        experiment(e).UnitROI{unit} = string(st.name(st.id == roiID));
    end
    
    % Get all session ID and types for the specified recording
    expname = {'orient', 'spafreq', 'loc', 'natimg', 'cwimg'};
    if strcmp(animal, 'JF088') && strcmp(day, '2022-11-29')
        expid = {'3', '3', '4', '5', '6'};
    elseif strcmp(animal, 'JF088') && strcmp(day, '2022-11-30')
        expid = {'2', '2', '3', '4', '1'};
    elseif strcmp(animal, 'JF089') && strcmp(site, 'site2')
        expid = {'5', '5', '6', '7', '8'};
    else
        expid = {'1', '1', '2', '3', '4'};
    end
    
    % Initialize the experiment(e).Data structure
    experiment(e).Data = struct('stimType', expname, 'spike_times', cell(1,5), ...
        'spike_templates', spike_templates, 'stimOn_t', cell(1,5), 'stim_vals', cell(1,5));
    
    % Loop across sessions and assign spike_times, stimOn_t, and stim_vals
    for s = 1:5
        % Load Timeline.mat and Block.mat for the specified session
        load(fullfile('D:\', animal, day, expid{s}, [day '_' expid{s} '_' animal '_Timeline.mat']));
        load(fullfile('D:\', animal, day, expid{s}, [day '_' expid{s} '_' animal '_Block.mat']));
        % Align spike times to the experiment
        experiment(e).Data(s).spike_times = align_ephys2timeline(spike_times, Timeline, sync);
        % Align timechunks to keep to the start of each experiment
        t_diff = experiment(e).Data(s).spike_times(1) - spike_times(1);
        experiment(e).Data(s).UnitTimeRange = trange + t_diff;
        % Get the stimulus values
        if strcmp(expname{s}, 'orient')
            experiment(e).Data(s).stim_vals = block.events.stimOrientationValues';
        elseif strcmp(expname{s}, 'spafreq')
            experiment(e).Data(s).stim_vals = block.events.stimSpatialFreqValues';
        else
            experiment(e).Data(s).stim_vals = block.events.stim_idValues';
        end
        % Get stimOn and stimOff times
        [experiment(e).Data(s).stimOn_t, experiment(e).Data(s).stimOff_t] = ...
            get_stimOn_t(Timeline);
        % Correct for stimOn_t and stim_vals mismatch
        if e == 1 && s == 4
            [~, ~, stimOn_t, stimOff_t] = get_stimOn_t(Timeline);
            experiment(e).Data(s).stimOn_t = stimOn_t(2:end);
            experiment(e).Data(s).stimOff_t = stimOff_t(2:end);
            experiment(e).Data(s).stim_vals = experiment(e).Data(s).stim_vals(2:end);
        elseif e == 3 && s == 3
            [~, ~, stimOn_t, stimOff_t] = get_stimOn_t(Timeline);
            experiment(e).Data(s).stimOn_t = stimOn_t(2:end);
            experiment(e).Data(s).stimOff_t = stimOff_t(2:end);
            experiment(e).Data(s).stim_vals = experiment(e).Data(s).stim_vals(2:end);
        elseif e == 6 && s == 4
            [~, ~, stimOn_t, stimOff_t] = get_stimOn_t(Timeline);
            experiment(e).Data(s).stimOn_t = stimOn_t(2:end);
            experiment(e).Data(s).stimOff_t = stimOff_t(2:end);
            experiment(e).Data(s).stim_vals = experiment(e).Data(s).stim_vals(2:end);
        end
    end
end

end