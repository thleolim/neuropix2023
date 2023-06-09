%% Extract channel spike count and their associated depths
function proc_chanspikes(kspath, site, neuropix_path)

config_path = fullfile(neuropix_path, 'files');
if strcmp(site, 'site1')
    chanMap = fullfile(config_path, 'bottrow0_chanmap.mat');
elseif strcmp(site, 'site2')
    chanMap = fullfile(config_path, 'bottrow50_chanmap.mat');
elseif strcmp(site, 'site1-shank0')
    chanMap = fullfile(config_path, 'shank0_chanmap.mat');
end
run(fullfile(config_path, 'configFile384.m'));

ops.trange = [0, Inf];
ops.NchanTOT = 385;
ops.chanMap = chanMap;
ops.minFR = 0;

fs = [dir(fullfile(kspath, '*.bin'))];
ops.fbinary = fullfile(kspath, fs.name);

[~, ich] = bd_preprocessDataSub(ops);
[channel_spikeCounts, channel_id] = groupcounts(ich);
channel_spikeCounts = gather(channel_spikeCounts);
channel_id = gather(channel_id);

writeNPY(channel_spikeCounts, fullfile(kspath, 'channels.spike_count.npy'));
writeNPY(channel_id, fullfile(kspath, 'channels.id.npy'));

end