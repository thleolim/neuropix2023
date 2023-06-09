%% Confirmatory Analysis %%
% Hodges-Lehmann estimate and confidence interval calculation was taken from:
% Cardillo G. (2006). Wilcoxon test: non parametric Wilcoxon test for paired samples.
% http://www.mathworks.com/matlabcentral/fileexchange/12702
clear; clc; close all
load("C:\Users\taeho\Desktop\output\variables\experiment.mat");
load("C:\Users\taeho\Desktop\output\variables\mean_fr.mat");

%% Get the index of all BG units
allbg_idx = zeros(0,2); % exp #, unit #
idx = 1;
for e = 1:6
    n_units = length(experiment(e).UnitID);
    for u = 1:n_units
        name = experiment(e).UnitROI(u);
        if strcmp(name{:}, 'CP') || strcmp(name{:}, 'SNr') || ...
                strcmp(name{:}, 'GPe')
            allbg_idx(idx,:) = [e u]; idx = idx + 1;
        end
    end
end

clearvars -except experiment mean_fr allbg_idx

%% 1. Wilcoxon sign-rank test for unit responsiveness
n_tests = 0;
% Initialize a cell array to store all WSR results of all BG units
wsr_result = cell(size(allbg_idx, 1), 6); % column 1: experiment ID
                                         % column 2: unit ID
                                         % column 3: WSR p-values
                                         % column 4: Hodges-Lehmann estimate
                                         % column 5: rank biserial correlation
                                         % column 6: H0 rejection
                                         % column 7: [W N]

% Conduct WSR across all BG unit item pairs
id = 1; by_mat = zeros(0, 6); % for Benjamini-Yekutieli correction later
for idx = 1:size(allbg_idx, 1)
    wsr_result{idx, 1} = allbg_idx(idx, 1); % experiment ID
    e = wsr_result{idx, 1};
    wsr_result{idx, 2} = allbg_idx(idx, 2); % unit ID
    u = wsr_result{idx, 2};
    wsr_result{idx, 3} = cell(1, 5); % for each stimulus type
    wsr_result{idx, 4} = cell(1, 5); % for each stimulus type
    wsr_result{idx, 5} = cell(1, 5); % for each stimulus type
    wsr_result{idx, 6} = cell(1, 5); % for each stimulus type
    wsr_result{idx, 7} = cell(1, 5); % for each stimulus type

    for s = 1:5
        if ~isnan(mean_fr{e,s}{u,1})
            n_stims = length(unique(experiment(e).Data(s).stim_vals));
            wsr_result{idx,3}{s} = zeros(n_stims, 1);
            wsr_result{idx,4}{s} = zeros(n_stims, 1);
            wsr_result{idx,5}{s} = zeros(n_stims, 1);
            wsr_result{idx,6}{s} = zeros(n_stims, 1);
            wsr_result{idx,7}{s} = zeros(n_stims, 2);
            
            for i = 1:n_stims

                % Get the p-value
                [wsr_result{idx,3}{s}(i), ~, stats] = signrank(mean_fr{e,s}{u,i}(:,2), ...
                    mean_fr{e,s}{u,i}(:,1), 'method', 'exact'); % baseline vs response period

                % Point estimate with confidence intervals
                diffs = sort(mean_fr{e,s}{u,i}(:,2) - mean_fr{e,s}{u,i}(:,1));
                diffs = diffs(diffs ~= 0); n = length(diffs);

                if n >= 35
                    [I,J] = ndgrid(diffs,diffs); 
                    d = triu(I+J)./2; % Walsh averages triang matrix
                    ld = sort(d(d~=0)); % Walsh averages vector
                    wsr_result{idx,4}{s}(i) = median(ld); % Hodges-Lehmann estimate
    
                    % Get the rank biserial correlation coefficient
                    sum_ranks = n*(n+1)/2; w1 = stats.signedrank; w2 = sum_ranks - w1;
                    wsr_result{idx,5}{s}(i) = abs(w1-w2) / sum_ranks;
                    wsr_result{idx,7}{s}(i,:) = [w1 n]; 
                else
                    wsr_result{idx,3}{s}(i) = 1;
                    wsr_result{idx,7}{s}(i,:) = [0 0]; 
                    wsr_result{idx,4}{s}(i) = 0;
                    wsr_result{idx,5}{s}(i) = 0;
                end

                n_tests = n_tests + 1;
                by_mat(id, 1) = e; by_mat(id, 2) = u;
                by_mat(id, 3) = s; by_mat(id, 4) = i; 
                by_mat(id, 5) = wsr_result{idx,3}{s}(i); 
                id = id + 1;
            end

        else % excluded by compute_timechunks
            wsr_result{idx,3}{s} = NaN;
            wsr_result{idx,4}{s} = NaN;
            wsr_result{idx,5}{s} = NaN;
            wsr_result{idx,6}{s} = NaN;
            wsr_result{idx,7}{s} = NaN;
        end
    end

end

clearvars -except experiment mean_fr wsr_result by_mat n_tests

%% Benjamini-Yekutieli correction across all WSR tests
[~, sortidx] = sort(by_mat(:,5)); 
by_mat = by_mat(sortidx, :); % order the matrix by rank
[~,~,rnk] = unique(by_mat(:,5)); % give the rank of each p-value
rnkedge = [1; diff(rnk)]; 
rnkedge = [0; rnkedge(2:end)-rnkedge(1:end-1)];
for idx = 1:length(rnk)
    if rnkedge(idx) == 1
        rnk(idx:end) = rnk(idx:end) + (idx - rnk(idx));
    end
end

% Conduct a stepwise Benjamini-Yekutieli procedure
m = sum(1./(1:length(rnk)));
for idx = 1:length(rnk)
    q = (rnk(idx)/length(rnk)) * 0.05/m;
    if by_mat(idx, 5) >= q
        by_mat(1:idx-1,6) = 1; % reject H0
        by_mat(idx:end,6) = 0; % fail to reject H0
        n_sigs = idx - 1;
        break
    end
end

% Assign the hypothesis decisions in the wsr_result structure
for idx = 1:length(rnk)
    e = by_mat(idx,1); u = by_mat(idx,2);
    s = by_mat(idx,3); i = by_mat(idx,4);
    for row = 1:size(wsr_result,1)
        if wsr_result{row,1} == e && wsr_result{row,2} == u
            this_row = row;
            break
        end
    end
    wsr_result{this_row,6}{s}(i) = by_mat(idx,6);
end

disp([num2str(n_sigs) '/' num2str(n_tests) ' were found to be significant.'])

clearvars -except experiment mean_fr wsr_result 

%% Find the number of items each unit-stimulus pair are responsive for
% Make an array that indicates the # of items each unit-stimulus is responsive for
wsr_mat = zeros(size(wsr_result,1), 5);
for idx = 1:size(wsr_mat, 1)
    for s = 1:5
        if ~isnan(wsr_result{idx,6}{s})
            wsr_mat(idx,s) = sum(wsr_result{idx,6}{s});
        else
            wsr_mat(idx,s) = NaN;
        end
    end
end

% Also store the sign-rank results of units responsive to only one item
% within a stimulus category
wsr_forone = cell(0, 7); % column 1: experiment ID
                         % column 2: unit ID
                         % column 3: stimulus ID
                         % column 4: item ID
                         % column 5: p-value
                         % column 6: Hodges-Lehmann estimate
                         % column 7: rank biserial correlation
                         % column 8: [W N]

id = 1;
for idx = 1:size(wsr_mat, 1)
    for s = 1:5
        if wsr_mat(idx,s) == 1
            wsr_forone{id,1} = wsr_result{idx,1};
            e = wsr_forone{id,1};
            wsr_forone{id,2} = wsr_result{idx,2};
            u = wsr_forone{id,2};
            wsr_forone{id,3} = s;
            wsr_forone{id,4} = find(wsr_result{idx,6}{s});
            i = wsr_forone{id,4};
            wsr_forone{id,5} = wsr_result{idx,3}{s}(i);
            wsr_forone{id,6} = wsr_result{idx,4}{s}(i);
            wsr_forone{id,7} = wsr_result{idx,5}{s}(i);
            wsr_forone{id,8} = wsr_result{idx,7}{s}(i,:);
            id = id + 1;
        end
    end
end

clearvars -except experiment mean_fr wsr_result wsr_mat wsr_forone

%% 2. Kruskal-Wallis test with post-hoc Bonferroni-corrected Dunn's test
kw_result = cell(0, 9); % column 1: experiment ID
                        % column 2: unit ID
                        % column 3: stimulus ID
                        % column 4: item IDs
                        % column 5: Hodge-Lehmann estimate
                        % column 6: [df chi2-statistic]
                        % column 7: p-value
                        % column 8: eta2 effect size
                        % column 9: stats for post-hoc testing
id = 1;
for idx = 1:size(wsr_mat, 1)
    for s = 1:5
        if wsr_mat(idx,s) > 1
            kw_result{id, 1} = wsr_result{idx, 1};
            e = kw_result{id, 1};
            kw_result{id, 2} = wsr_result{idx, 2};
            u = kw_result{id, 2};
            kw_result{id, 3} = s;
            kw_result{id, 4} = find(wsr_result{idx, 6}{s});
            i = kw_result{id, 4};
            kw_result{id, 5} = zeros(length(i), 1);

            tempmat = zeros(0, 2);
            for item = 1:length(i)
                tempmat = vertcat(tempmat, [repmat(i(item), ...
                    size(mean_fr{e,s}{u,i(item)},1),1) ...
                    mean_fr{e,s}{u,i(item)}(:,2)-mean_fr{e,s}{u,i(item)}(:,1)]);

                x = mean_fr{e,s}{u,i(item)}(:,2)-mean_fr{e,s}{u,i(item)}(:,1);
                [I,J] = ndgrid(x,x); 
                d = triu(I+J)./2; % Walsh averages triang matrix
                ld = sort(d(d~=0)); % Walsh averages vector
                kw_result{id,5}(item) = median(ld);
            end

            [kw_result{id,7}, tbl, stats] = kruskalwallis(tempmat(:,2), ...
                tempmat(:,1), 'off');
            kw_result{id,6} = [tbl{2,3} tbl{2,5}];
%             kw_result{id,8} = tbl{2,5}*(size(tempmat,1)+1)/(size(tempmat,1)^2 - 1);
            kw_result{id,8} = abs((tbl{2,5} - numel(i) + 1)/(size(tempmat,1) - numel(i)));
            kw_result{id,9} = stats;

            id = id + 1;
        end
    end
end

ph_result = cell(0, 6); % column 1: experiment ID
                        % column 2: unit ID
                        % column 3: stimulus ID
                        % column 4: sig-different item pairs
                        % column 5: Dunn's z-test statistic
                        % column 6: Bonferroni-adjusted p-value
id = 1;
for idx = 1:size(kw_result, 1)
    if kw_result{idx,7} < 0.05
        c = multcompare(kw_result{idx,9}, 'display', 'off', 'ctype', 'bonferroni');
        sigpairs = find(c(:,6) < 0.05);
        if ~isempty(sigpairs)
            ph_result{id,1} = kw_result{idx,1};
            ph_result{id,2} = kw_result{idx,2};
            ph_result{id,3} = kw_result{idx,3};
            for p = 1:numel(sigpairs)
                ph_result{id,4} = kw_result{idx,4}(c(sigpairs, 1:2));
                ph_result{id,5} = c(sigpairs, 4);
                ph_result{id,6} = c(sigpairs, 6);
            end
            id = id + 1;
        end
    end
end

clearvars -except wsr_result wsr_mat wsr_forone kw_result ph_result

%% Save variables for figures and tables in other scripts
save("C:\Users\taeho\Desktop\output\results.mat")