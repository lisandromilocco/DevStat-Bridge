clear;
close all;

% Basic color palette
colors = lines(5);
c1 = colors(1,:);
c2 = colors(5,:);

% Create figure
figure('Color', 'w', 'Position', [100, 100, 1200, 300]);
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'loose'); % 3 plots

% Files and labels
meanFiles = {'N_individuals_Mean.mat', 'Gain_p_Mean.mat', 'Gain_m_Mean.mat'}; %this file has all the datapoints of the 100 replicates
allFiles  = {'N_individuals_ALL.mat',  'Gain_p_ALL.mat',  'Gain_m_ALL.mat'}; %this one has the parameter values in the first column
xtickScales = [1, 100, 100];
xlabels = {'Number of individuals', 'Developmental noise (%)', 'Measurement noise (%)'};
ylims = {[0 320], [0 320], [0 700]};
ks_all = 1:4; % NUMBER OF CONDITIONS PER PARAMETER

% Loop through panels
for i = 1:3
    nexttile;
    
    % Load data
    load(meanFiles{i});
    load(allFiles{i});
    
    hold on;
    
    offset = 0; % counter for the position of boxplots; gets updated inside loop
    jump = 1; % how close each pair of boxplots are with each other
    thickness = 0.35; % how close the 2 boxplots are within treatment
    jitterAmount = 0.08;

    tickLabels = store_Error(ks_all, 1); % first column has the parameter value for x-axis

    for idx = 1:length(ks_all)
        k = ks_all(idx);
        data = store_ALL{k};
        
        x = [data(:,1); data(:,3)]; % Columns 1 and 3 have the static and Kalman filtered estimates; Column 2 is the naive moving average
        group = [ones(size(data,1),1)*offset; ones(size(data,1),1)*(offset + thickness)];

        boxplot(x, group, ...
            'Positions', [offset, offset + thickness], ...
            'Colors', [c1; c2], ...
            'Symbol', '', 'Widths', 0.3, 'Whisker', 1.5);

        h = findobj(gca, 'Tag', 'Box'); % searches for the boxes of the boxplot to add the transparency
        patch(get(h(1), 'XData'), get(h(1), 'YData'), c2, ...
              'FaceAlpha', 0.2, 'EdgeColor', c2, 'LineWidth', 1.5);
        patch(get(h(2), 'XData'), get(h(2), 'YData'), c1, ...
              'FaceAlpha', 0.2, 'EdgeColor', c1, 'LineWidth', 1.5);

        % Style whiskers
        allLines = findobj(gca, 'Type', 'Line');
        whiskers = allLines(arrayfun(@(h) length(h.YData) == 2, allLines));
        set(whiskers, 'LineStyle', '-', 'LineWidth', 1.5); % 
        % forces whiskers to be continous lines instead of dashed; this affects all lines!! in the plot
                                                            

        % Add jittered scatter
        scatter(offset + (rand(size(data,1),1) - 0.5)*2*jitterAmount, data(:,1), ...
                10, c1, 'filled', 'MarkerFaceAlpha', 0.2);
        scatter(offset + thickness + (rand(size(data,1),1) - 0.5)*2*jitterAmount, data(:,3), ...
                10, c2, 'filled', 'MarkerFaceAlpha', 0.2);

        offset = offset + jump; %update location of boxplots
    end

    % Axes styling
    xlim([-0.5, offset - jump + thickness + 0.5]);
    ylim(ylims{i});
    xticks(thickness/2:jump:offset);
    xticklabels(string(xtickScales(i) .* tickLabels));
    
    set(gca, 'FontSize', 14, 'LineWidth', 1.5, 'Box', 'off', ...
             'TickDir', 'out', 'FontName', 'Arial');

    xlabel(xlabels{i}, 'FontSize', 21, 'FontName', 'Arial');
    if i == 1
        ylabel('Relative error (%)', 'FontSize', 21, 'FontName', 'Arial');
    end
end

saveas(gcf, 'Panels.svg');
