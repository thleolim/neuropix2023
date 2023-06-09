% Get mean_fr cell array
load("C:\Users\taeho\Desktop\output\variables\experiment.mat");

mean_fr = cell(6, 5); % nExperiment X nStimulusType

% Define the baseline and response windows
base_win = [-0.4 0]; resp_win = [0 0.2];
base_dur = diff(base_win); resp_dur = diff(resp_win);

for e = 1:6 % experiments
    tic
    % Get experiment data
    expdata = experiment(e).Data;
    n_units = length(experiment(e).UnitID);

    for s = 1:5 % stimulus Type session
        stims = unique(expdata(s).stim_vals); n_stims = length(stims);
        mean_fr{e,s} = cell(n_units, n_stims);

        for i = 1:n_stims % stimulus item
            % Extract the stimulus item's stimOn times
            this_stimTimes = (expdata(s).stim_vals == stims(i));
            this_stimTimes = expdata(s).stimOn_t(this_stimTimes);
            n_trial = length(this_stimTimes);
            % Get the basleine and response time edges
            base_edge = base_win + this_stimTimes;
            resp_edge = this_stimTimes + resp_win;

            for u = 1:n_units % units
                if this_stimTimes(end) > experiment(e).Data(s).UnitTimeRange(u,1) && ...
                        this_stimTimes(1) < experiment(e).Data(s).UnitTimeRange(u,2)
                    % Ultimately, each entry will be an array of two columns
                    % across trials: first being baseline firing rate and
                    % second being response firing rate
                    mean_fr{e,s}{u,i} = zeros(n_trial, 2);
                    % Extract the unit spike times
                    this_spikeTimes = (expdata(s).spike_templates == u);
                    this_spikeTimes = expdata(s).spike_times(this_spikeTimes);
                    % Count window spike density and divide by window length
                    mean_fr{e,s}{u,i}(:,1) = arrayfun(@(t) histcounts(...
                        this_spikeTimes, [base_edge(t,1) base_edge(t,2)]), ...
                        1:n_trial)/base_dur;
                    mean_fr{e,s}{u,i}(:,2) = arrayfun(@(t) histcounts(...
                        this_spikeTimes, [resp_edge(t,1) resp_edge(t,2)]), ...
                        1:n_trial)/resp_dur;
                else
                    mean_fr{e,s}{u,i} = NaN;
                end
            end
        end
    end
    toc
end
clearvars -except mean_fr
save("C:\Users\taeho\Desktop\output\variables\mean_fr.mat");