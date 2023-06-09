%% Count Spike Density for Each Channel for a Specified Shank %%
% This function counts the number of spikes of a given experiment across
% all channels of a specified shank. If a shank0 recording exists and shank
% 1 is specified, only the shank0 recording is used. Otherwise, a
% combination of available site1/site2 recordings are used.
% Input:
%   site1_ks: path to site1 kilosort output
%   site2_ks: path to site2 kilosort ouptut
%   shank0: path to site0 SpikeGLX output.
%   (site1_ks/site2_ks/shank0 = [] if it doesn't exist)
%   shank: which shank unit information to extract (1, 2, 3, 4)
%   neuropix_path: path to the leo_neuropix package
% Output:
%   chan_scounts: count of spikes in each channel
%   scount_depths: depth of each corresponding channel in chan_scounts

function [chan_scounts, scount_depths] = count_spikes(site1_ks, site2_ks, shank0, shank, neuropix_path)

config_path = fullfile(neuropix_path, 'files');

if shank == 1 && ~isempty(shank0)
    chan_scounts = readNPY(fullfile(shank0, 'channels.spike_count.npy'));
    channel_id = readNPY(fullfile(shank0, 'channels.id.npy'));
    load(fullfile(config_path, "shank0_chanmap.mat"), 'ycoords');
    ycoords = 2865 - ycoords(1:end-1,1); scount_depths = ycoords(channel_id);
else
    if ~isempty(site1_ks)
        chan_scounts1 = readNPY(fullfile(site1_ks, 'channels.spike_count.npy'));
        channel_id1 = readNPY(fullfile(site1_ks, 'channels.id.npy'));
        load(fullfile(config_path, "bottrow0_chanmap.mat"), 'ycoords');
        load(fullfile(config_path, "bottrow0_chanmap.mat"), 'xcoords');
        ycoords = 2865 - ycoords(1:end-1,1);
        ycoords = ycoords(channel_id1); xcoords = xcoords(channel_id1);
        thisShank = [(shank - 1)*250 (shank - 1)*250 + 32];
        shankPos = logical(sum(xcoords == thisShank, 2));
        ycoords1 = ycoords(shankPos); 
        chan_scounts1 = chan_scounts1(shankPos);
    end
    if ~isempty(site2_ks)
        chan_scounts2 = readNPY(fullfile(site2_ks, 'channels.spike_count.npy'));
        channel_id2 = readNPY(fullfile(site2_ks, 'channels.id.npy'));
        load(fullfile(config_path, "bottrow50_chanmap.mat"), 'ycoords');
        load(fullfile(config_path, "bottrow50_chanmap.mat"), 'xcoords');
        ycoords = 2865 - ycoords(1:end-1,1);
        ycoords = ycoords(channel_id2); xcoords = xcoords(channel_id2);
        thisShank = [(shank - 1)*250 (shank - 1)*250 + 32];
        shankPos = logical(sum(xcoords == thisShank, 2));
        ycoords2 = ycoords(shankPos); 
        chan_scounts2 = chan_scounts2(shankPos);
    end
    if ~isempty(site1_ks) && ~isempty(site2_ks)
        chan_scounts = [chan_scounts1; chan_scounts2];
        scount_depths = [ycoords1; ycoords2];
    elseif ~isempty(site1_ks)
        chan_scounts = chan_scounts1;
        scount_depths = ycoords1;
    elseif ~isempty(site2_ks)
        chan_scounts = chan_scounts2;
        scount_depths = ycoords2;
    end
end

% Sum up channels equal in depth
eqdepth = [-15; diff(scount_depths)]; scount_depths(eqdepth==0, 1) = NaN;
for i = 2:length(chan_scounts)
    if isnan(scount_depths(i))
        chan_scounts(i) = chan_scounts(i) + chan_scounts(i-1);
        chan_scounts(i-1) = NaN; 
    end
end
chan_scounts = chan_scounts(~isnan(chan_scounts));
scount_depths = scount_depths(~isnan(scount_depths));

end