function [stimOn_t, stimOff_t, stimOn_t2, stimOff_t2] = get_stimOn_t(Timeline)

photodiode_idx = strcmp({Timeline.hw.inputs.name}, 'photoDiode');

% Get Timeline timestamps for when photodiode is ON 
stimScreen_on = Timeline.rawDAQData(:, photodiode_idx) > 0.2;
stimScreen_on_t = Timeline.rawDAQTimestamps(stimScreen_on);

% Downsample (with median) the photodiode trace
photodiode_trace_medfilt = medfilt1(Timeline.rawDAQData(stimScreen_on, photodiode_idx), 3);
% Convolve the photodiode trace to reveal its delayed differential
photodiode_diff_t = 20; % in ms
photodiode_diff_samples = round(Timeline.hw.daqSampleRate/1000*photodiode_diff_t);
photodiode_diff_filt = [1, zeros(1, photodiode_diff_samples), -1];
photodiode_diff_conv = abs(conv(photodiode_trace_medfilt, photodiode_diff_filt, 'valid'));

% Binarize the photodiode trace into stim ON/OFF
photodiode_diff_thresh = range(Timeline.rawDAQData(:, photodiode_idx)) * 0.2;
photodiode_trace_diff = photodiode_diff_conv > photodiode_diff_thresh;

% Get the timestamps when the photodiode detected a flip in ON/OFF
photodiode_flip = find(~photodiode_trace_diff(1:end-1) & ...
    photodiode_trace_diff(2:end)) + photodiode_diff_samples + 1;
photodiode_flip_times = stimScreen_on_t(photodiode_flip)';

% Get the timestamps for stimOn
stimOn_t = photodiode_flip_times(2:2:end-1);
stimOn_t2 = photodiode_flip_times(1:2:end);
% Additional: stimOff times
stimOff_t = photodiode_flip_times(3:2:end);
stimOff_t2 = photodiode_flip_times(2:2:end);

end