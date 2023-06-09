%% 3. Electrophysiology to Histology Alignment Pipeline %%
% Because dyed probe tracks in histology slices do not accurately locate
% the coordinates of each channel, channel coordinates are instead
% estimated by aligning electrophysiology signatures to brain regions
% identified along the probe trajectory interpolated from Braindraw.
clear; clc; close all

%% Initialization
% This project utilizes three recording configurations that is not always
% all present for each probe insertion session: (1) site1 (bottom 50 rows
% of channels), (2) site2 (next bottom 50 rows of channels), and (3)
% site1-shank0 (all channels of shank 1).
% Specify the paths to each recording below and set to [] if unavailable.
site1_bc = 'C:\Users\taeho\Desktop\output\bombcell_output\JF090_2022-11-16_site1';
site1_ks = 'D:\recordings\JF090_2022-11-16_site1';
site2_bc = [];
site2_ks = [];
shank0 = [];
% Path to neuropix2023 package
neuropix_path = 'C:\Users\taeho\Desktop\neuropix2023';
% Specify which shank to align (from 1 to 4)
shank = 4;
% Load Braindraw output, i.e. probe_ccf.mat
animal = 'JF090'; day = '2022-11-16';
probe_ccf_path = fullfile('C:\Users\taeho\Desktop\output\braindraw_output\', ...
    animal, day, 'probe_ccf.mat');
load(probe_ccf_path)

%% Get the normalized good unit firing rate along shank depth
% Returns all good unit firing rate (normalized to maximum firing rate of
% all units) and their corresponding depth relative to the top-most
% channel's depth. 
[norm_fr, depths] = unitfr(site1_bc, site1_ks, site2_bc, site2_ks, shank);
% This helps user align depths to regions with neuron bodies, i.e. not
% ventricles or white matter, and high-firing units to corresponding
% regions, e.g. SNr, and vice versa. 

%% Get depth-binned unit firing rate correlation matrix
% Finds good and multi- unit spike density along the shank depth and
% discretizes them into bins. A correlation matrix is then computed by
% finding the Pearson's correlation coefficient between each bins.
[mua_corr, mua_depths] = get_mua(site1_bc, site1_ks, site2_bc, site2_ks, shank);
% "Clusters" of high correlation along the diagonal suggests coherent
% regions of similar spiking density. 

%% Get channel spike counts along respective channel depths
% Using the preprocessing steps of Kilosort2, the spike counts for each
% channel is computed. This step requires GPU and stores its output in the
% braindraw output folder.
if ~isempty(site1_ks) && isempty(dir(fullfile(site1_ks, 'channels.id.npy')))
    proc_chanspikes(site1_ks, 'site1', neuropix_path);
end
if ~isempty(site2_ks) && isempty(dir(fullfile(site2_ks, 'channels.id.npy')))
    proc_chanspikes(site2_ks, 'site2', neuropix_path);
end
if ~isempty(shank0) && isempty(dir(fullfile(shank0, 'channels.id.npy')))
    proc_chanspikes(shank0, 'site1-shank0', neuropix_path);
end

[chan_scounts, scount_depths] = count_spikes(site1_ks, site2_ks, shank0, shank, neuropix_path);

%% Electrophysiology-to-Histology Alignment GUI
% Load the relevant files for displaying the regions along the probe trajectory
[~, ~, st, ~] = bd_loadAllenAtlas(fullfile(neuropix_path, 'files\allen_mouse_25um_v1.2'));
load(fullfile(neuropix_path, 'files\allen_ccf_colormap_2017.mat'));
% Use the gathered electrophysiology markers and align them to the probe
% trajectory to then interpolate each good unit region.
ephys2hist_gui(shank, site1_ks, site2_ks, shank0, norm_fr, depths, mua_corr, ...
    mua_depths, chan_scounts, scount_depths, st, cmap, probe_ccf, probe_ccf_path);