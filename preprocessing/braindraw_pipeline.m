%% 2. Braindraw Probe Tracking Pipeling %%
% Interpolates the probe trajectories in the Allen Brain Atlas from
% manually-traced probe tracks of atlas-registered histology slices. This
% requires two 25um-downsampled .tif stacks: (1) a green fluorescence for
% atlas registration and (2) a red fluorescence channel for probe tracking.
clear; clc; close all

%% Initialization
% Path to registratation channel .tif stack
regimg_dir = dir(['D:\histology\JF090\downsampled_stacks\025_micron\*','green','*.tif*']);
regchan = [regimg_dir.folder, filesep, regimg_dir.name];
% Path to probe tracking channel .tif stack
traimg_dir = dir(['D:\histology\JF090\downsampled_stacks\025_micron\*','red','*.tif*']);
trachan = [traimg_dir.folder, filesep, traimg_dir.name];
% Load Allen mouse atlas (25um) for registration
neuropix_path = 'C:\Users\taeho\Desktop\neuropix2023';
[tv, av, st, bregma] = bd_loadAllenAtlas(fullfile(neuropix_path, 'files', 'allen_mouse_25um_v1.2'));
% Braindraw output folder
animal = 'JF090'; day = '2022-11-16';
savePath = fullfile('C:\Users\taeho\Desktop\output\braindraw_output', animal, day);
if isempty(dir(savePath))
    mkdir(savePath);
end

%% Brainreg: Automatic Atlas Registration
% This section requires (1) anaconda installation (https://www.anaconda.com/) 
% and (2) brainreg venv setup; run the following in anaconda prompt:
%   > conda create -n brainreg python
%   > conda activate brainreg
%   > pip install brainreg
orientation = 'psr'; atlas = 'allen_mouse_25um';
CMD = sprintf('brainreg %s %s -a %s -v 25 25 25 --orientation %s --atlas %s', ...
    regchan, savePath, trachan, orientation, atlas);
clipboard('copy', CMD);
% CMD is now saved in your clipboard. Paste and run CMD in brainreg venv to
% automatically register the two channels to the Allen Brain Atlas.

%% Braindraw: Check Registration Output
check = 1; % EDIT flag to run the script below
if check
    % Load registered .tiff stack
    regimg = loadtiff([savePath, filesep, 'downsampled_standard.tiff']);
    bd_convertToAPFormat(regimg, tv, av, savePath); % convert to standard AP format
    % Check registered atlas orientation
    bd_checkAndCorrectAtlasOrientation(tv, av, st, regimg, savePath, 2);
    % Check registered atlas alignment
    bd_checkAndCorrectAtlasAlignment(tv, av, st, regimg, savePath, 2);
end

%% Braindraw: Manual Probe Tracking
% Load the colormap file used in the allenCCF repository
load(fullfile(neuropix_path, 'files', 'allen_ccf_colormap_2017.mat'));
% Load the registered .tif stack into the manual probe tracking GUI
traimg_dir = dir([savePath, filesep, 'downsampled_standard_*.tiff']);
traimg = loadtiff([traimg_dir.folder, filesep, traimg_dir.name]);
% Braindraw GUI
probetrack_gui(cmap, tv, av, st, traimg, savePath);

% GUI Guide:
%   - Specify the number of shanks used, e.g. 4 for Neuropixels 2.0
%   - Use left/right arrow keys to go through slices. Press the 'Auto
%     brightness/ contrast' button to normalize slice brightness for dye.
%   - Press 'Probe 1' or 1 on the keyboard to trace for shank 1 and etc.
%     Press 'esc' to return. Click on the start and end of the probe track
%     in a slice to trace. Do this for each slice with probe tracks of
%     interest, e.g. 8-10 slices from the shank tip.
%   - Press 'esc' to save or quit all progress.
%   - Braindraw interpolates a linear probe trajectory through the 3D Allen
%     Brain Atlas from each slice's probe track. Inspect the output figure
%     and make sure the interpolated probe trajectories look right. 