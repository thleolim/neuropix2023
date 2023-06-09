%% Load Unit Spiking Activity in a Defined Experiment %%
% This function is used to extract spiking information of defined units of
% an experiment in indices or time. Channel positions are adjusted to make
% it understandable and good/+mua units are extracted.
% Input:
%   bcpath: path to bombcell output folder
%   kspath: path to kilosort output folder
%   site: 'site1' or 'site2'
%   flag: 1 for good + mua units, 0 for good units only
% Output: 
%   spike_times: vector of all spiking time in an experiment
%   spike_templates: vector of corresponding template index to spike_times
%   template_depths: vector of depth to corresponding template index
%   template_amplitudes: vector of amplitude to corresponding template index
%   spike_depths: vector of corresponding template depths to spike_times
%   spike_xdepths: vector of corresponding template x-coord to spike_times
%   template_xdepths: vector of x-coord to corresponding template index
%   channel_positions: array of all channel-matched x-coord and depth
%   templates: array of template amplitudes across channels
%   trange: two column array of start and stop times to use

function [spike_times, spike_templates, template_depths, template_amplitudes, ...
    spike_depths, spike_xdepths, template_xdepths, channel_positions, templates, ...
    trange] = load_spikes(bcpath, kspath, site, flag)

% Load bombcell data
qMetric = parquetread(fullfile(bcpath, 'templates._bc_qMetrics.parquet'));
param = parquetread(fullfile(bcpath, '_bc_parameters._bc_qMetrics.parquet'));
trange = [qMetric.useTheseTimesStart qMetric.useTheseTimesStop];
unitType = bc_getQualityUnitType(param, qMetric);

% Load kilosort data
spike_times = double(readNPY(fullfile(kspath, 'spike_times.npy'))) ./ 30000;
spike_templates_0idx = readNPY(fullfile(kspath, 'spike_templates.npy'));
templates_whitened = readNPY(fullfile(kspath, 'templates.npy'));
channel_positions = readNPY(fullfile(kspath, 'channel_positions.npy'));
winv = readNPY(fullfile(kspath, 'whitening_mat_inv.npy'));
template_amplitudes = readNPY(fullfile(kspath, 'amplitudes.npy'));

% Flush channel depths to 0, invert, then reposition channel depths
% As a result, channel depth increases down the shank
if strcmp(site, 'site1')
    channel_positions(:, 2) = channel_positions(:, 2) - min(channel_positions(:,2));
    channel_positions(:, 2) = 2865 - channel_positions(:, 2);
elseif strcmp(site, 'site2')
    channel_positions(:, 2) = channel_positions(:, 2) - min(channel_positions(:,2));
    channel_positions(:, 2) = 2145 - channel_positions(:, 2);
end

% Unwhiten templates
templates = zeros(size(templates_whitened));
for t = 1:size(templates_whitened, 1)
    templates(t, :, :) = squeeze(templates_whitened(t, :, :)) * winv;
end

% Get channel_position idx for each template
[~, max_site] = max(max(abs(templates), [], 2), [], 3);
template_xdepths = channel_positions(max_site, 1);
% Get min-max range of each template in each template
template_chan_amp = squeeze(range(templates, 2));
% Filter out low amplitude channels for each template
template_chan_amp_thresh = max(template_chan_amp, [], 2) * 0.5;
template_chan_amp_overthresh = template_chan_amp .* (template_chan_amp >= ...
    template_chan_amp_thresh);
% Get center-of-mass on thresholded channel amplitudes
template_depths = sum(template_chan_amp_overthresh.*channel_positions(:, 2)', 2) ...
    ./ sum(template_chan_amp_overthresh, 2);
% Get depth of each spike
spike_depths = template_depths(spike_templates_0idx+1);
spike_xdepths = template_xdepths(spike_templates_0idx+1);

% Find all the good units defined by quality metrics
if flag % both mua and good units
    good_templates = unitType == 1 | 2;
    good_templates_idx = find(unitType == 1 | 2) - 1;
else % only good units
    good_templates = unitType == 1;
    good_templates_idx = find(unitType == 1) - 1;
end

% Throw out all non-good template data
template_depths = template_depths(good_templates);
template_xdepths = template_xdepths(good_templates);
templates = templates(good_templates,:,:);
trange = trange(good_templates, :);

% Throw out all non-good spike data
good_spike_idx = ismember(spike_templates_0idx, good_templates_idx);
spike_times = spike_times(good_spike_idx);
template_amplitudes = template_amplitudes(good_spike_idx);
spike_depths = spike_depths(good_spike_idx);
spike_xdepths = spike_xdepths(good_spike_idx);
spike_templates_0idx = spike_templates_0idx(good_spike_idx);

% Rename the spike templates according to the remaining templates
new_spike_idx = nan(max(spike_templates_0idx)+1, 1);
new_spike_idx(good_templates_idx+1) = 1:length(good_templates_idx);
spike_templates = new_spike_idx(spike_templates_0idx+1);

end