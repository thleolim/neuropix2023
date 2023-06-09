%% Example Histogram of Unit 4.26 %%
clear; clc; close all
load("C:\Users\taeho\Desktop\output\variables\mean_fr.mat");
load("C:\Users\taeho\Desktop\output\results.mat");
e = 4; u = 26; s = 1; i = 3;

%% Distributions
b_meanfr = mean_fr{e,s}{u,i}(:,1);
r_meanfr = mean_fr{e,s}{u,i}(:,2);
d_meanfr = r_meanfr - b_meanfr;

%% Histograms
tiledlayout(2,1,'TileSpacing', 'compact', 'Padding', 'compact')
nexttile
histogram(b_meanfr, 'BinWidth', 5); hold on;
histogram(r_meanfr, 'BinWidth', 5); 
xlim([-20, 60])
nexttile
histogram(d_meanfr, 'BinWidth', 5, 'FaceColor','y');
xline(10, 'r--'); % pseudomedian
xlim([-20, 60])
xlabel('mean firing rate (sp/s)')