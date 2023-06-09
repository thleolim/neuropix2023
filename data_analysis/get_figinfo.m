%% Helper code for figure information %%
%% Location Images
clear; clc; close all
thispath = 'C:\Users\taeho\Desktop\output\stimuli\location';
f = figure('Position', [1 1 1200 300]); 
tiledlayout(2, 3, 'TileSpacing', 'Compact')
for i = 1:6
    nexttile
    load(fullfile(thispath, ['img' num2str(i) '.mat']))
    imagesc(img)
    colormap('gray')
    pbaspect([3 1 1])
    xticks([])
    yticks([])
end

saveas(f, 'C:\Users\taeho\Desktop\output\location.tif');
close(f)

%% Natural Images
clear; clc; close all
thispath = 'C:\Users\taeho\Desktop\output\stimuli\natimg';
f = figure('Position', [1 1 1200 1500]); 
tiledlayout(10, 3, 'TileSpacing', 'Compact')
for i = 1:30
    nexttile
    load(fullfile(thispath, ['img' num2str(i) '.mat']))
    imagesc(img)
    colormap('gray')
    pbaspect([3 1 1])
    xticks([])
    yticks([])
end

saveas(f, 'C:\Users\taeho\Desktop\output\natimg.tif');
close(f)

%% Blank Images
clear; clc; close all
thispath = 'C:\Users\taeho\Desktop\output\stimuli\location';
f = figure('Position', [1 1 1200 300]); 
load(fullfile(thispath, ['img' num2str(1) '.mat']))
thismin = min(min(img)); thismax = max(max(img));
img(1:125,:) = 128;
img(265:303, 1037:1067) = 255;
imagesc(img);
caxis([thismin thismax]);
colormap('gray')
pbaspect([3 1 1])
caxis
xticks([])
yticks([])

saveas(f, 'C:\Users\taeho\Desktop\output\blank.tif');
close(f)

%% Find median stimulation and trial duration
clear; clc; close all
load("C:\Users\taeho\Desktop\output\variables\experiment.mat");
load("C:\Users\taeho\Desktop\output\variables\mean_fr.mat");


all_stimOn = [];
all_stimOff = [];
exp_length = zeros(6, 4);
for e = 1:6
    for s = 2:5
        all_stimOn = [all_stimOn; experiment(e).Data(s).stimOn_t];
        all_stimOff = [all_stimOff; experiment(e).Data(s).stimOff_t];
        exp_length(e,s-1) = experiment(e).Data(s).stimOff_t(end) - experiment(e).Data(s).stimOn_t(1);
    end
end

trial_t = diff(all_stimOn); trial_t = trial_t(trial_t > 0);
stim_t = all_stimOff - all_stimOn; stim_t = stim_t(stim_t > 0);
disp(['Median trial duration is ' num2str(median(trial_t)) 's']);
disp(['Minimum trial duration is ' num2str(min(trial_t)) 's']);
disp(['Maximum trial duration is ' num2str(max(trial_t)) 's']);
disp(['Median stimulation duration is ' num2str(median(stim_t)) 's']);
disp(['Minimum stimulation duration is ' num2str(min(stim_t)) 's']);
disp(['Maximum stimulation duration is ' num2str(max(stim_t)) 's']);

disp(['Mean gratings experiment duration was ' num2str(mean(exp_length(:,1))/60) 'min']);
disp(['Mean location experiment duration was ' num2str(mean(exp_length(:,2))/60) 'min']);
disp(['Mean natimg experiment duration was ' num2str(mean(exp_length(:,3))/60) 'min']);
disp(['Mean cwimg experiment duration was ' num2str(mean(exp_length(:,4))/60) 'min']);

%% Get number of units for each category
load("C:\Users\taeho\Desktop\output\variables\experiment.mat");
save("C:\Users\taeho\Desktop\output\unitclass.mat")

n_gunits = 0;
for e = 1:6
    n_gunits = n_gunits + length(experiment(e).UnitID);
end

n_cp = 0; n_snr = 0; n_gpe = 0;
for e = 1:6
    n_units = length(experiment(e).UnitID);
    for u = 1:n_units
        name = experiment(e).UnitROI(u);
        if strcmp(name{:}, 'CP')
            n_cp = n_cp + 1;
        elseif strcmp(name{:}, 'SNr')
            n_snr = n_snr + 1;
        elseif strcmp(name{:}, 'GPe')
            n_gpe = n_gpe + 1;
        end
    end
end

n_rcp = 0; n_rsnr = 0;
for idx = 1:size(unit_char,1)
    if strcmp(unit_char{idx,6}, 'SNr')
        n_rsnr = n_rsnr + 1;
    else
        n_rcp = n_rcp + 1;
    end
end
