%% Get Multiunit Correlation Matrices %%
% This function bins unit spike density by depth and cross-correlates each
% bin. This should provide clusters of high correlation. 
% Input:
%   site1_bc: path to folder containing site1 bombcell output
%   site1_ks: path to folder containing site1 kilosort output
%   site2_bc: path to folder containing site2 bombcell output
%   site2_ks: path to folder containing site2 kilosort output
%   (site1_bc/site1_ks/site2_bc/site2_ks = [] if it doesn't exist)
%   shank: which shank unit information to extract (1, 2, 3, 4)
% Output:
%   depth_centres: depth coordinates of the bin centres
%   mua_corr: an nxn correlation matrix of n bins

function [mua_corr, depth_centres] = get_mua(site1_bc, site1_ks, site2_bc, site2_ks, shank)

if ~isempty(site1_bc)
    [spike_times1, spike_templates1, template_depths1, ~, ~, spike_xdepths1, ...
        template_xdepths1] = load_spikes(site1_bc, site1_ks, 'site1', 1);
    % Define channel position/ template index/ spike times of defined shank
    theseChannelPositions1 = [(shank-1) * 250, (shank-1)*250 + 32];
    theseTemplates1 = ismember(template_xdepths1, theseChannelPositions1);
    theseSpikes1 = ismember(spike_xdepths1, theseChannelPositions1);
    % Redefine spike times/ template index of defined shank
    spike_times1 = spike_times1(theseSpikes1);
    spike_templates1 = spike_templates1(theseSpikes1);
    good_templates_idx1 = unique(spike_templates1); 
    new_spike_idx1 = nan(max(spike_templates1), 1); 
    new_spike_idx1(good_templates_idx1) = 1:length(good_templates_idx1); 
    spike_templates1 = new_spike_idx1(spike_templates1); 
    % Get array of corresponding unit depth
    template_depths1 = template_depths1(theseTemplates1);
end

if ~isempty(site2_bc)
    [spike_times2, spike_templates2, template_depths2, ~, ~, spike_xdepths2, ...
        template_xdepths2] = load_spikes(site2_bc, site2_ks, 'site2', 1);
    % Define channel position/ template index/ spike times of defined shank
    theseChannelPositions2 = [(shank-1) * 250, (shank-1)*250 + 32];
    theseTemplates2 = ismember(template_xdepths2, theseChannelPositions2);
    theseSpikes2 = ismember(spike_xdepths2, theseChannelPositions2);
    % Redefine spike times/ template index of defined shank
    spike_times2 = spike_times2(theseSpikes2);
    spike_templates2 = spike_templates2(theseSpikes2);
    good_templates_idx2 = unique(spike_templates2); 
    new_spike_idx2 = nan(max(spike_templates2), 1); 
    new_spike_idx2(good_templates_idx2) = 1:length(good_templates_idx2); 
    spike_templates2 = new_spike_idx2(spike_templates2); 
    % Get array of corresponding unit depth
    template_depths2 = template_depths2(theseTemplates2);
end

% Assign spike_times, spike_templates, template_depths according to site1
% and site2 availability
if ~isempty(site1_bc) && isempty(site2_bc) % only site1 avliable
    spike_times = spike_times1;
    spike_templates = spike_templates1;
    template_depths = template_depths1;
elseif isempty(site1_bc) && ~isempty(site2_bc) % only site2 available
    spike_times = spike_times2;
    spike_templates = spike_templates2;
    template_depths = template_depths2;
elseif ~isempty(site1_bc) && ~isempty(site2_bc) % both site1 site2 available
    spike_times = [spike_times2; spike_times1];
    spike_templates = [spike_templates2; spike_templates1 + length(template_depths2)];
    template_depths = [template_depths2; template_depths1];
end

n_corr_groups = 20; % # of bins in correlation matrix (should be even)
spike_binning = 0.01; % in sec

x = range(template_depths) / (n_corr_groups + 1 - 1);
min_depth = min(template_depths) - .5*x; max_depth = max(template_depths) + .5*x;
depth_group_edges = linspace(min_depth, max_depth, n_corr_groups + 1);

depth_group = discretize(template_depths, depth_group_edges);
depth_centres = depth_group_edges(1:end-1) + (diff(depth_group_edges)/2);
unique_depths = 1:length(depth_group_edges) - 1;

corr_edges = nanmin(spike_times):spike_binning:nanmax(spike_times);

binned_spikes_depth = zeros(length(unique_depths),length(corr_edges)-1);
for curr_depth = 1:length(unique_depths)
    binned_spikes_depth(curr_depth,:) = histcounts(spike_times( ...
        ismember(spike_templates,find(depth_group == unique_depths(curr_depth)))), ...
        corr_edges);
end

mua_corr = corrcoef(binned_spikes_depth');

if ~isempty(site1_bc) && ~isempty(site2_bc)
    mua_corr_temp = mua_corr;
    mua_corr = zeros(n_corr_groups, n_corr_groups/2);
    mua_corr(1:n_corr_groups/2, :) = mua_corr_temp(1:n_corr_groups/2, 1:n_corr_groups/2);
    mua_corr(n_corr_groups/2+1:end, :) = mua_corr_temp(n_corr_groups/2+1:end, n_corr_groups/2+1:end);
end

end