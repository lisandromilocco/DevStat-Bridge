%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
format short e
beep off

load ("predictions.mat")
n_replays=1

%%%%%%
for i = 1:size(mu_values,1)
    for rep = 1:n_replays
        obs = squeeze(observed_deltaZ_store(i,rep,:));
        lande = (deltaZ_Lande_store(i,:));
        naive = (deltaZ_naive_store(i,:));

        % Relative errors
        relErr_Lande(i,rep) = 100 * norm(obs - lande) / norm(obs);
        relErr_naive(i,rep) = 100 * norm(obs - naive) / norm(obs);

        % Angle errors
        cosine_val = dot(obs, lande) / (norm(obs) * norm(lande));
        cosine_val = min(1, max(-1, cosine_val)); % clamp to avoid NaN
        angleErr_Lande(i,rep) = acos(cosine_val) * 180/pi;

        cosine_val = dot(obs, naive) / (norm(obs) * norm(naive));
        cosine_val = min(1, max(-1, cosine_val));
        angleErr_naive(i,rep) = acos(cosine_val) * 180/pi;

        % Store observed angle
        store_ObservedChange(i, rep) = atan2(obs(2), obs(1));
    end
end



figure;
figure('Color', 'w', 'Position', [100, 100, 500, 200]); % White background and set the aspect ratio

n_conditions = 3
tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');

% Define colors

colores = lines(5); % This generates a set of 5 distinct colors
colors.obs   = [0 0 0]%
colors.lande = colores(1,:);
colors.naive = colores(2,:);


index = [1 5 10]

for i = 1:n_conditions
    nexttile;
    hold on; axis equal;
    title(sprintf('\\mu = %.3f', mu_values(index(i),2)), 'FontWeight','bold','FontSize',14);

    % Get vectors for this condition
    obs_all   = squeeze(observed_deltaZ_store(index(i),:,:));   % [n_replays, 2]
    lande = (deltaZ_Lande_store(index(i),:));      % [n_replays, 2]
    naive = (deltaZ_naive_store(index(i),:));      % [n_replays, 2]

    % Normalize to unit length (avoid division by zero)
    obs_all   = obs_all ./ max(vecnorm(obs_all,2,2), eps);
    lande = lande ./ max(vecnorm(lande,2,2), eps);
    naive = naive ./ max(vecnorm(naive,2,2), eps);

    % --- Plot observed vectors  ---
    for k = 1:n_replays
        quiver(0, 0, obs_all(k,1), obs_all(k,2), 0, ...
               'Color', "k", ... % transparency
               'LineWidth', 1.5, 'MaxHeadSize', 0.6);
    end


    % --- Plot predicted vectors  ---
    quiver(0, 0, lande(1), lande(2), 0, ...
           'Color', colors.lande, 'LineWidth', 2.2, 'MaxHeadSize', 0.8);
    quiver(0, 0, naive(1), naive(2), 0, ...
           'Color', colors.naive, 'LineWidth', 2.2, 'MaxHeadSize', 0.8);

    % Axes settings
    xlim([-1.1 1.1]); ylim([-1.1 1.1]);
        xticks([-1 0 1]); yticks([-1 0 1]);
    xlabel('x','FontSize',13); ylabel('y','FontSize',13);
    set(gca,'FontSize',13,'FontName','Arial','LineWidth',1.5, ...
            'Box','off','TickDir','out');
        xlabel('Trait 1','FontSize',13); ylabel('Trait 2','FontSize',13);

         set(gca,'FontSize',13,'FontName','Arial','LineWidth',1.5, ...
            'Box','on','TickDir','out');  % <-- Box ON to get full frame


    % Median angle errors (shorter labels, stacked neatly)
    med_lande = median(angleErr_Lande(index(i),:),'omitnan');
    med_naive = median(angleErr_naive(index(i),:),'omitnan');
    text(-1.0,0.9, sprintf('Using G: %.1f°', med_lande), ...
         'Color', colors.lande, 'FontSize',12);
    text(-1.0,0.7, sprintf('Using P: %.1f°', med_naive), ...
         'Color', colors.naive, 'FontSize',12);
end


% Legend at the bottom (short, professional labels)
lgd = legend({'Observed', ...
              'Predicted (G)','Predicted (P≈G)'}, ...
             'Orientation','horizontal','Box','off','FontSize',13);
lgd.Layout.Tile = 'south';












return



figure;
n_conditions = size(mu_values, 1);

% Customize the axes
ax = gca;
ax.FontSize = 12; % Modern, consistent font size
ax.FontName = 'Arial'; % Clean, professional font
ax.LineWidth = 1.5; % Slightly thicker axes lines
ax.Box = 'off'; % Remove box for a minimal look
ax.TickDir = 'out'; % Ticks pointing outward for clarity



for i = 1:n_conditions
    subplot(ceil(n_conditions/2), 2, i); % adjust layout based on number of mu_values
    hold on; axis equal; grid on;
    title(sprintf('\\mu = %.3f', mu_values(i,2)));
    
    % Get all replays for this i
    obs_all   = squeeze(observed_deltaZ_store(i,:,:));   % size = [n_replays, 2]
    lande = (deltaZ_Lande_store(i,:));      % size = [n_replays, 2]
    naive = (deltaZ_naive_store(i,:));      % size = [n_replays, 2]
    
    % Normalize to unit norm
    obs_all   = obs_all ./ vecnorm(obs_all,2,2);
    lande = lande ./ vecnorm(lande,2,2);
    naive = naive ./ vecnorm(naive,2,2);
    
    % Plot observed vectors (blue)
    quiver(zeros(size(obs_all,1),1), zeros(size(obs_all,1),1), ...
           obs_all(:,1), obs_all(:,2), 0, 'Color', [0 0.447 0.741], ...
           'LineWidth', 1.5, 'MaxHeadSize', 0.3);
    
    % Plot Lande prediction (red)
    quiver(zeros(size(lande,1),1), zeros(size(lande,1),1), ...
           lande(:,1), lande(:,2), 0, 'Color', [0.85 0.325 0.098], ...
           'LineWidth', 1.5, 'MaxHeadSize', 0.3);
    
    % Plot naive prediction (green)
    quiver(zeros(size(naive,1),1), zeros(size(naive,1),1), ...
           naive(:,1), naive(:,2), 0, 'Color', [0.466 0.674 0.188], ...
           'LineWidth', 1.5, 'MaxHeadSize', 0.3);
    
    xlim([-1.1 1.1]); ylim([-1.1 1.1]);
    xlabel('x'); ylabel('y');
    %legend({'Observed','Lande','Naive'}, 'Location', 'bestoutside');
end


%%%%%


% Save to file
%save('matrices.mat', 'P_store', 'G_store');


figure;
% Example: using 25th and 75th percentiles
lower_q = 0.25;   % Lower quantile (25%)
upper_q = 0.75;   % Upper quantile (75%)

% Compute median and quantiles for angleErr_Lande
med_Lande = median(angleErr_Lande, 2);
qL_Lande = quantile(angleErr_Lande, lower_q, 2);
qU_Lande = quantile(angleErr_Lande, upper_q, 2);

% Convert quantiles to asymmetric error bars (distance from median)
err_lower_Lande = med_Lande - qL_Lande;
err_upper_Lande = qU_Lande - med_Lande;

% Compute median and quantiles for angleErr_naive
med_naive = median(angleErr_naive, 2);
qL_naive = quantile(angleErr_naive, lower_q, 2);
qU_naive = quantile(angleErr_naive, upper_q, 2);

err_lower_naive = med_naive - qL_naive;
err_upper_naive = qU_naive - med_naive;

% Plot with asymmetric error bars
errorbar(mu_values(:,2), med_Lande, err_lower_Lande, err_upper_Lande, ...
    'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName','Lande'); hold on;

errorbar(mu_values(:,2), med_naive, err_lower_naive, err_upper_naive, ...
    's-', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName','Naive');
xlabel('\mu_{\theta2}');
ylabel('Angle Error (degrees)');
legend('Location','best');
title(sprintf('Directional Error (mean ± SD) over %d replays', n_replays));

ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.LineWidth = 1.5;
ax.Box = 'off';
ax.TickDir = 'out';
set(gca, 'XScale', 'log');

% %% CHANGE: add a little horizontal padding
xlim([min(mu_values(:,2))*0.9, max(mu_values(:,2))*1.1]);  




figure;
plot(mu_values(:,2), angleErr_Lande(:,:), ...
    'o-', 'LineWidth', 2, 'MarkerSize', 8, ...   %% CHANGE: thicker line, bigger markers
    'DisplayName','Lande'); hold on;
plot(mu_values(:,2), angleErr_naive(:,:), ...
    's-', 'LineWidth', 2, 'MarkerSize', 8, ...   %% CHANGE
    'DisplayName','Naive');
xlabel('\mu_{\theta2}');
ylabel('Angle Error (degrees)');
legend('Location','best');
title(sprintf('Directional Error (mean ± SD) over %d replays', n_replays));

ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.LineWidth = 1.5;
ax.Box = 'off';
ax.TickDir = 'out';
set(gca, 'XScale', 'log');

% %% CHANGE: add a little horizontal padding
xlim([min(mu_values(:,2))*0.9, max(mu_values(:,2))*1.1]);  



figure;
errorbar(mu_values(:,2), mean(relErr_Lande,2), std(relErr_Lande,0,2), ...
    'o-', 'LineWidth', 2, 'MarkerSize', 8, ...   %% CHANGE
    'DisplayName','Lande'); hold on;
errorbar(mu_values(:,2), mean(relErr_naive,2), std(relErr_naive,0,2), ...
    's-', 'LineWidth', 2, 'MarkerSize', 8, ...   %% CHANGE
    'DisplayName','Naive');
xlabel('\mu_{\theta2}');
ylabel('Relative Error (%)');
legend('Location','best');
title(sprintf('Relative Error (mean ± SD) over %d replays', n_replays));

ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.LineWidth = 1.5;
ax.Box = 'off';
ax.TickDir = 'out';
set(gca, 'XScale', 'log');

% %% CHANGE: add a little horizontal padding
xlim([min(mu_values(:,2))*0.9, max(mu_values(:,2))*1.1]);  