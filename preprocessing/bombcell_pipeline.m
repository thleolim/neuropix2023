%% 1. Bombcell Quality Metrics Pipeline %%
% Automatically classify spike clusters, "units", into good/ noise/ multi.
% This requires a Kilosort output folder containing:
%   (1) spike_times.npy (spikeTimes_samples)
%   (2) spike_templates.npy (spikeTemplates)
%   (3) templates.npy (templateWaveforms)
%   (4) amplitudes.npy (templateAmplitudes)
%   (5) pc_features.npy (pcFeatures)
%   (6) pc_feature_ind.npy (pcFeaturesIdx)
%   (7) channel_positions.npy (channelPositions)
% Edit bc_qualityParamValues.m to adjust classification parameters.
clear; clc; close all

%% Initialization
% Kilosort output folder
kspath = 'D:\recordings\JF088_2022-11-29_site1';
% Bombcell output folder
animal = 'JF088'; day = '2022-11-29'; site = 'site1';
savePath = fullfile('C:\Users\taeho\Desktop\output\bombcell_output', animal, day, site);
if isempty(dir(savePath))
    mkdir(savePath); 
end
% Path to SpikeGLX output files
binpath = fullfile(kspath, '2022-11-29_JF088_g0_t0.imec0.ap.bin'); 
metadir = dir(fullfile(kspath, '2022-11-29_JF088_g0_t0.imec0.ap.meta'));
% Parameters
param = bc_qualityParamValues(metadir, binpath);

%% Run Bombcell
% Load the .npy Kilosort output files as MATLAB variables
[spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, ...
    pcFeatures, pcFeaturesIdx, channelPositions] = bc_loadEphysData(kspath);

run_bombcell = 1; % set to 1 if you wish to run Bombcell
                  % set to 0 if you wish to load existing Bombcell output

if run_bombcell
    param.computeDrift = 1; param.minSpatialDecaySlope = -0.007;
    [qMetric, unitType] = bc_runAllQualityMetrics(param, spikeTimes_samples, ...
        spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
        pcFeaturesIdx, channelPositions, savePath);
else % load bombcell
    [param, qMetric] = bc_loadSavedMetrics(savePath);
    param.minSpatialDecaySlope = -0.007;
    unitType = bc_getQualityUnitType (param, qMetric);
end

% Unit type: 0 - noise units (not a neuron)
%            1 - good units (likely a single neuron)
%            2 - multi units (multiple neurons)
%            3 - non-somatic units (non-somatic spikes)

%% Visualize Bombcell output
% GUI settings
loadRawTraces = 1;
bc_loadMetricsForGUI;

% Open GUI
unitQualityGuiHandle = bc_unitQualityGUI(memMapData, ephysData, qMetric, ...
    forGUI, rawWaveforms, param, probeLocation, unitType, loadRawTraces);

% GUI Guide
% left/right arrow: toggle between units 
% g : go to next good unit 
% m : go to next multi-unit 
% n : go to next noise unit
% up/down arrow: toggle between time chunks in the raw data
% u: brings up a input dialog to enter the unit you want to go to