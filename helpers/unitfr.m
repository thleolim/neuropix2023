%% Get Unit Normalized Firing Rate and Their Depths %%
% This function retrieves the normalized firing rate of all good units and
% their corresponding depth for the ephys2hist GUI.
% Input:
%   site1_bc: path to folder containing site1 bombcell output
%   site1_ks: path to folder containing site1 kilosort output
%   site2_bc: path to folder containing site2 bombcell ouptut
%   site2_ks: path to folder containing site2 kilosort ouptut
%   (site1_bc/site1_ks/site2_bc/site2_ks = [] if it doesn't exist)
%   shank: which shank unit information to extract (1, 2, 3, 4)
% Output:
%   normfr: array of normalized firing rate of each unit
%   depths: array of depths of each corresponding unit

function [normfr, depths] = unitfr(site1_bc, site1_ks, site2_bc, site2_ks, shank)

% Extract site1 good unit spikes of defined shank
if ~isempty(site1_bc)
    [~, spike_templates1, template_depths1, ~, ~, spike_xdepths1, ...
        template_xdepths1] = load_spikes(site1_bc, site1_ks, 'site1', 0);
    % Define channel position/ template index/ spike times of defined shank
    theseChannelPositions1 = [(shank-1) * 250, (shank-1)*250 + 32];
    theseTemplates1 = ismember(template_xdepths1, theseChannelPositions1);
    theseSpikes1 = ismember(spike_xdepths1, theseChannelPositions1);
    % Redefine spike times/ template index of defined shank
    spike_templates1 = spike_templates1(theseSpikes1);
    good_templates_idx1 = unique(spike_templates1); 
    new_spike_idx1 = nan(max(spike_templates1), 1); 
    new_spike_idx1(good_templates_idx1) = 1:length(good_templates_idx1); 
    spike_templates1 = new_spike_idx1(spike_templates1); 
    % Get array of spike count of each unit index
    [~,~,spike_templates_reidx1] = unique(spike_templates1);
    template_spike_n1 = accumarray(spike_templates_reidx1,1);
    % Get array of corresponding unit depth
    template_depths1 = template_depths1(theseTemplates1);
end

% Extract site2 good unit spikes of defined shank
if ~isempty(site2_bc)
    [~, spike_templates2, template_depths2, ~, ~, spike_xdepths2, ...
        template_xdepths2] = load_spikes(site2_bc, site2_ks, 'site2', 0);
    % Define channel position/ template index/ spike times of defined shank
    theseChannelPositions2 = [(shank-1) * 250, (shank-1)*250 + 32];
    theseTemplates2 = ismember(template_xdepths2, theseChannelPositions2);
    theseSpikes2 = ismember(spike_xdepths2, theseChannelPositions2);
    % Redefine spike times/ template index of defined shank
    spike_templates2 = spike_templates2(theseSpikes2);
    good_templates_idx2 = unique(spike_templates2); 
    new_spike_idx2 = nan(max(spike_templates2), 1); 
    new_spike_idx2(good_templates_idx2) = 1:length(good_templates_idx2); 
    spike_templates2 = new_spike_idx2(spike_templates2); 
    % Get array of spike count of each unit index
    [~,~,spike_templates_reidx2] = unique(spike_templates2);
    template_spike_n2 = accumarray(spike_templates_reidx2,1);
    % Get array of corresponding unit depth
    template_depths2 = template_depths2(theseTemplates2);
end

% Find the log-normalized firing rate of each unit
if ~isempty(site1_bc) && isempty(site2_bc) % only site1 avliable
    normfr = mat2gray(log10(template_spike_n1+1));
elseif isempty(site1_bc) && ~isempty(site2_bc) % only site2 available
    normfr = mat2gray(log10(template_spike_n2+1));
elseif ~isempty(site1_bc) && ~isempty(site2_bc) % both site1 site2 available
    ntot = sum(template_spike_n1); ntot2 = sum(template_spike_n2);
    if ntot2 > ntot
        template_spike_n2 = template_spike_n2 * (ntot/ntot2);
    elseif ntot > ntot2
        template_spike_n1 = template_spike_n1 * (ntot2/ntot);
    end
    template_spike_n = [template_spike_n2; template_spike_n1];
    normfr = mat2gray(log10(template_spike_n+1));
end

% Find the depth of each corresponding unit 
if ~isempty(site1_bc) && isempty(site2_bc) % only site1 avliable
    depths = template_depths1;
elseif isempty(site1_bc) && ~isempty(site2_bc) % only site2 available
    depths = template_depths2;
elseif ~isempty(site1_bc) && ~isempty(site2_bc) % both site1 site2 available
    depths = [template_depths2; template_depths1];
end

end