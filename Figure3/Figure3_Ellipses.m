clear
close all
format short e
beep off

load('matrices_10.mat');  % This loads P_all and G_all back into the workspace

angle_store = zeros(size(G_store,3),1);

colors = lines(5); % This generates a set of 5 distinct colors


for i=1:size(G_store,3)

    P=P_store(:,:,i);
    G=G_store(:,:,i);

    norm(P)
end

figure;
for i=1:size(G_store,3)

    P=P_store(:,:,i);
    G=G_store(:,:,i);

    P = P ./ norm(P);
    G = G ./ norm(G);

    subplot(2, 5, i); % 2 rows and 4 columns, for the i-th subplot

    plotEllipse(P, 0, 0, colors(2,:));
    hold on;
    plotEllipse(G, 0, 0, colors(1,:));


    %Compute the eigenvalues and eigenvectors for G and P
    [eigVecG, eigValG] = eig(G);
    [eigVecP, eigValP] = eig(P);

    %Extract the main eigenvectors (the ones corresponding to the largest eigenvalues)
    [~, idxG] = max(diag(eigValG)); % Index of largest eigenvalue for G
    [~, idxP] = max(diag(eigValP)); % Index of largest eigenvalue for P

    mainEigenVecG = eigVecG(:, idxG); % Main eigenvector of G
    mainEigenVecP = eigVecP(:, idxP); % Main eigenvector of P

    %Calculate the angle between the two eigenvectors using the dot product
    cosTheta = dot(mainEigenVecG, mainEigenVecP) / (norm(mainEigenVecG) * norm(mainEigenVecP));
    theta = acos(cosTheta); % Angle in radians
    theta_deg = rad2deg(theta); % Convert angle to degrees

    %Add legend with angle between eigenvectors
%     legendStr = sprintf('Angle: %.2f', theta_deg);
%     legend('show');
%     legend('Location', 'best');
%     legend('string', {legendStr});

    angle_store(i)=theta_deg;

    axis equal;
    grid on;
%     xlabel('X-axis');
%     ylabel('Y-axis');
%     title('Covariance Ellipse');
    xlim([-2 2])
    ylim([-2 2])
end

saveas(gca,'Figure3_ellipses.svg','svg')

%% 

col2=[9 8 7 6 5 4 3 2 1 0]';
col2= 0.5 ./ 2.^col2;


%figure(2);
figure('Color', 'w', 'Position', [100, 100, 450, 300]); % White background and set the aspect ratio

plot(col2', angle_store, 'Color', "black", 'LineWidth', 2);
hold on;
scatter(col2', angle_store, 60, 'MarkerFaceColor', "black", 'MarkerFaceAlpha', 0.6, ...
    'MarkerEdgeAlpha', 0.6, 'MarkerEdgeColor', 'black');
% Customize the axes
ax = gca;
ax.FontSize = 12; % Modern, consistent font size
ax.FontName = 'Arial'; % Clean, professional font
ax.LineWidth = 1.5; % Slightly thicker axes lines
ax.Box = 'off'; % Remove box for a minimal look
ax.TickDir = 'out'; % Ticks pointing outward for clarity

% Add labels to the axes
xlabel('Allelic frequency', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('Angle (degrees)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Extend the x and y axis limits for better visibility
xRange = max(col2) - min(col2);
yRange = max(angle_store) - min(angle_store);
% xlim([min(col2) * 0.8, max(col2) + 0.2*xRange]);
ylim([min(angle_store) - 0.1*yRange, max(angle_store) + 0.1*yRange]);

xlim([0.00065 0.6]);


% Set tight layout
outerpos = ax.OuterPosition;
ti = ax.TightInset; 


set(gca, 'XScale', 'log');
saveas(gca,'Figure3_angles.svg','svg')
