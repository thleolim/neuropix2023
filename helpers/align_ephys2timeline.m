function spike_times_aligned = align_ephys2timeline(spike_times, Timeline, sync, exp)

% Create electrophysiology timestamps based on sync
ephys_sample_rate = 30000;
ephys_timestamps = [0:size(sync, 2) - 1] / ephys_sample_rate;
% Binarize sync into flipper ON/OFF
flip_threshold = range(sync)/2;
ephys_binary_values = sync > flip_threshold;
% Get electrophysiology flipper flip times
ephys_sync_values = find((~ephys_binary_values(1:end-1) & ephys_binary_values(2:end)) | ...
    (ephys_binary_values(1:end-1) & ~ephys_binary_values(2:end))) + 1;
ephys_sync_timestamps = ephys_timestamps(ephys_sync_values)';

% Get flipper flip times from Timeline now
flipper_idx = strcmp({Timeline.hw.inputs.name}, 'flipper');
flipper_thresh = 2; % TTL threshold
flipper_trace = Timeline.rawDAQData(:, flipper_idx) > flipper_thresh;
flipper_flip = find((~flipper_trace(1:end-1) & flipper_trace(2:end)) | ...
        (flipper_trace(1:end-1) & ~flipper_trace(2:end))) + 1;
timeline_flip_times = Timeline.rawDAQTimestamps(flipper_flip)';

% Infer experiment # by comparing each sync blocks with Timeline length
flip_diff_thresh = 1; % flipper flip length < 1 sec
flipper_expt_idx = [1; find(abs(diff(ephys_sync_timestamps)) > ...
    flip_diff_thresh) + 1; length(ephys_sync_timestamps)+1];
[~, infer_expid] = min(abs(diff(flipper_expt_idx) - length(timeline_flip_times)));

% Get electrophysiology flipper flip times of inferred/provided experiment block 
if ~exist('exp', 'var')
    exp = infer_expid;
end
if exp == infer_expid % use provided experiment #
    ephys_flip_times = ephys_sync_timestamps(...
        flipper_expt_idx(exp):flipper_expt_idx(exp+1)-1);
else                  % use inferred experiment #
    ephys_flip_times = ephys_sync_timestamps(...
        flipper_expt_idx(infer_expid):flipper_expt_idx(infer_expid+1)-1);
    disp('Provided experiment ID was incorrect!');
end

% Align spike times to Timeline by interpolating with flipper flip times
spike_times_aligned = interp1(ephys_flip_times, timeline_flip_times, spike_times, 'linear', 'extrap');

end