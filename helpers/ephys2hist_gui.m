%% Electrophysiology-to-Histology Alignment GUI %%
% This loads a GUI that facilitates the alignment of Neuropixel 2.0 shank
% channels to brain regions. It does so by aligning electrophysiology 
% signatures (unit firing rate, MUA correlation, channel spike density) 
% to the interpolated probe trajectory generated from braindraw.
% Inputs:
%   shank: a value between 1-4 representing the shank that will be aligned
%   site1/site2/shank0: path representing their kilosort/spikeglx output.
%                       [] represents their absence. 
%   norm_fr: the normalized firing rates of all good units across the shank
%   depths: depth along the shank with index representing norm_fr unit
%   mua_corr: a correlation matrix of multiunit spike density binned across depth
%   mua_depths: the depth along the shank of each mua_corr bin centres
%   chan_scounts: spike counts of an experiment session across channels
%   scount_depths: depth of each channel with index of chan_scounts
%   st: a table of allenCCF region index, name, and acronym
%   cmap: an RGB colormap matrix for allenCCF structures
%   probe_ccf: braindraw output with rows representing shanks 1-4
%   path: path to the probe_ccf.mat file

function ephys2hist_gui(shank, site1, site2, shank0, norm_fr, depths, mua_corr, ...
    mua_depths, chan_scounts, scount_depths, st, cmap, probe_ccf, path)

gui_fig = figure('color', 'w', 'KeyPressFcn', @keypress);

% Conditional to figure out axis ranges for recorded channels (in um)
if shank == 1 && ~isempty(shank0)
    min_depths = 0; max_depths = 2865;
else
    if ~isempty(site1)
        max_depths = 2865;
        if ~isempty(site2)
            min_depths = 1440;
        else
            min_depths = 2160;
        end
    elseif ~isempty(site2)
        min_depths = 1440; max_depths = 2145;
    end
end

% Plot normalized good units firing rate across depth
unit_ax = subplot('Position', [0.1,0.1,0.1,0.8]);
scatter(norm_fr, depths, 15, 'k', 'filled');
set(unit_ax,'YDir','reverse');
ylim([min_depths, max_depths]);
xlabel('N spikes')
title('Template depth & rate')
set(unit_ax, 'FontSize', 12)
ylabel('Depth (\mum)');

% Plot the multiunit correlation matrix
multiunit_ax = subplot('Position', [0.2,0.1,0.3,0.8]);
imagesc(mua_depths, mua_depths, mua_corr);
ylim([min_depths, max_depths]);
set(multiunit_ax, 'YTick',[]);
title('MUA correlation');
set(multiunit_ax, 'FontSize',12)
xlabel(multiunit_ax, 'Multiunit depth');

% Plot the spike density across channels
sd_ax = subplot('Position', [0.5,0.1,0.3,0.8]);
plot(smoothdata(chan_scounts, 'movmedian', 2), scount_depths)
set(sd_ax, 'YTick', []);
title(['5'' shank' num2str(shank) ' spike count' ]);
set(sd_ax, 'FontSize', 12)
xlabel(sd_ax, '# threshold crossings'); 
makepretty;

% Plot the probe trajectory
probe_areas_ax = subplot('Position', [0.8,0.1,0.05,0.8]);
% Sort probe trajectory depths by the dorsal-most coordinate
[~,dv_sort_idx] = sort(probe_ccf(shank).trajectory_coords(:,2));
% Interpolate trajectory depths from coords and convert to um
probe_trajectory_depths = ...
    pdist2(probe_ccf(shank).trajectory_coords, ...
    probe_ccf(shank).trajectory_coords((dv_sort_idx == 1),:))*25;
% Find each area boundary and center in terms of actual trajectory depths
trajectory_area_boundary_idx = ...
    [1;find(diff(double(probe_ccf(shank).trajectory_areas)) ~= 0)+1];
trajectory_area_boundaries = probe_trajectory_depths(trajectory_area_boundary_idx);
trajectory_area_centers = (trajectory_area_boundaries(1:end-1) + ...
    diff(trajectory_area_boundaries)/2);
for iArea = 1:size(trajectory_area_boundary_idx, 1)
    trajectory_area_labels(iArea) = st.name(st.id == ...
                probe_ccf(shank).trajectory_areas(trajectory_area_boundary_idx(iArea)));
end
[~,area_dv_sort_idx] = sort(trajectory_area_centers);
% Get the ROI labels and their centre depth coordinates
trajectory_area_centers = trajectory_area_centers(area_dv_sort_idx);
trajectory_area_labels = trajectory_area_labels(area_dv_sort_idx);
last_label = max(probe_trajectory_depths)/100;
last_label = floor(last_label); last_label = last_label*100;
yticks = sort([[0:100:last_label]'; trajectory_area_centers]);
ylabels = cell(1, length(yticks)); t = 1;
for i = 1:length(ylabels)
    if ismember(yticks(i), trajectory_area_centers)
        ylabels{i} = trajectory_area_labels{t}; t = t + 1;
    else
        ylabels{i} = yticks(i);
    end
end
% Finally plot the probe trajectory and their ROIs
image([],probe_trajectory_depths,probe_ccf(shank).trajectory_areas);
colormap(probe_areas_ax,cmap);
caxis([1,size(cmap,1)])
set(probe_areas_ax, 'YTick', probe_trajectory_depths);
set(probe_areas_ax,'YTick',yticks,'YTickLabels',ylabels);
set(probe_areas_ax,'XTick',[]);
set(probe_areas_ax,'YAxisLocation','right')
ylim([min_depths,max_depths]);
ylabel({'Probe areas','(Arrow/shift keys to move)','(Escape: save & quit)'});
set(probe_areas_ax,'FontSize',10)

% Draw boundary lines at borders (and undo clipping to extend across all)
boundary_lines = gobjects;
for curr_boundary = 1:length(trajectory_area_boundaries)
    boundary_lines(curr_boundary,1) = line(probe_areas_ax,[-13.5,1.5], ...
        repmat(trajectory_area_boundaries(curr_boundary),1,2),'color','b','linewidth',2);
end
set(probe_areas_ax,'Clipping','off');

% Package into gui
gui_data = struct;
gui_data.probe_ccf_fn = path; 

gui_data.probe_ccf = probe_ccf;
gui_data.shank = shank;

gui_data.unit_ax = unit_ax;
gui_data.multiunit_ax = multiunit_ax;

gui_data.probe_areas_ax = probe_areas_ax;
gui_data.probe_areas_ax_ylim = ylim(probe_areas_ax);
gui_data.probe_trajectory_depths = probe_trajectory_depths;

% Upload gui data
guidata(gui_fig, gui_data);

end

%% Function to control keypress events
function keypress(gui_fig, eventdata)

% Get guidata
gui_data = guidata(gui_fig);

% Set amounts to move by with/without shift
if any(strcmp(eventdata.Modifier,'shift'))
    y_change = 100;
    s_change = 0.1;
else
    y_change = 1;
    s_change = 0.01;
end

switch eventdata.Key
    
    % up/down: move probe areas
    case 'uparrow'
        new_ylim = gui_data.probe_areas_ax_ylim - y_change;
        ylim(gui_data.probe_areas_ax,new_ylim);
        gui_data.probe_areas_ax_ylim = new_ylim;
        % Upload gui data
        guidata(gui_fig,gui_data);
    case 'downarrow'
        new_ylim = gui_data.probe_areas_ax_ylim + y_change;
        ylim(gui_data.probe_areas_ax,new_ylim);
        gui_data.probe_areas_ax_ylim = new_ylim;
        % Upload gui data
        guidata(gui_fig,gui_data);
        
    % escape: save and quit
    case 'escape'
        opts.Default = 'Yes';
        opts.Interpreter = 'tex';
        user_confirm = questdlg('\fontsize{15} Save and quit?','Confirm exit',opts);
        if strcmp(user_confirm,'Yes')
            
            probe_ccf = gui_data.probe_ccf;
            
            % Get the probe depths corresponding to the trajectory areas
            probe_depths = gui_data.probe_areas_ax_ylim;  
            probe_ccf(gui_data.shank).probe_depths = probe_depths;
            
            % Save the appended probe_ccf structure
            save(gui_data.probe_ccf_fn, 'probe_ccf');
            
            % Close the figure
            close(gui_fig);
        end

     case 'n' % stretch - narrow
           new_ylim = [gui_data.probe_areas_ax_ylim(1), gui_data.probe_areas_ax_ylim(2) .*(1-s_change)];
           ylim(gui_data.probe_areas_ax, new_ylim);
           gui_data.probe_areas_ax_ylim = new_ylim;
           % Upload gui data
           guidata(gui_fig,gui_data);

    case 'w' % stretch - widden
           new_ylim = [gui_data.probe_areas_ax_ylim(1), gui_data.probe_areas_ax_ylim(2) .*(1+s_change)];
           ylim(gui_data.probe_areas_ax, new_ylim);
           gui_data.probe_areas_ax_ylim = new_ylim;
           % Upload gui data
           guidata(gui_fig,gui_data);
end

end