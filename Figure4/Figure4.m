clear
close all

load('sensitivities.mat')

%% -------------------------------------------------------------
%% First figure with sensitivity vectors
%% -------------------------------------------------------------
% This figure shows each parameter sensitivity vector projected
% around the fixed reference phenotype in trait space.


%% Reference phenotype
% This is the output from the tooth model. The first element is the 
% height of the central cusp, second and third and x-y position of
% posterior cusp, and fourht and fifth elements are x-y position of
% anterior cusp
ref_phenotype = [ ...
   6.3065e+00   3.7396e+00   4.0343e+00  -4.2140e+00   3.7562e+00];


figure('Color', 'w', 'Position', [100, 100, 400, 300]); 
hold on;
colors=lines(length(folders));

colors = [
    0.00, 0.45, 0.74;  
    0.85, 0.33, 0.10;  
    0.93, 0.69, 0.13;  
    0.47, 0.67, 0.19;  
    0.30, 0.75, 0.93;  
    0.64, 0.08, 0.18;  
    0.49, 0.18, 0.56;  
    1.00, 0.41, 0.71;  
    0.50, 0.50, 0.50;  
];

axis equal
hold on;

for parameter = 1:length(folders)

    s=sensitivities(2:3, 42, parameter);
    s=3*s/norm(s);
    plot([ref_phenotype(2) ref_phenotype(2)+s(1)], ...
         [ref_phenotype(3) ref_phenotype(3)+s(2)], ...
         LineWidth=2.5,Color=colors(parameter,:));

    s=sensitivities(4:5, 42, parameter);
    s=3*s/norm(s);
    plot([ref_phenotype(4) ref_phenotype(4)+s(1)], ...
         [ref_phenotype(5) ref_phenotype(5)+s(2)], ...
         LineWidth=2.5,Color=colors(parameter,:));

end

plot(ref_phenotype(2),ref_phenotype(3), '.', 'MarkerSize', 30);
plot(ref_phenotype(4),ref_phenotype(5), '.', 'MarkerSize', 30);

ax = gca;
ax.FontSize = 12; 
ax.FontName = 'Arial'; 
ax.LineWidth = 1.5; 
ax.Box = 'off'; 
ax.TickDir = 'out'; 
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
ax.Position = [outerpos(1) + ti(1), outerpos(2) + ti(2), ...
               outerpos(3) - ti(1) - ti(3), ...
               outerpos(4) - ti(2) - ti(4) - 0.05];


saveas(gcf, 'Figure4_sensitivities.svg');


%% -------------------------------------------------------------
%% 2) Compute angles between sensitivity vectors
%% -------------------------------------------------------------
% This section computes pairwise angles between parameter
% sensitivity vectors (time index 42).

angles_matrix = zeros(length(folders), length(folders));

for i = 1:length(folders)
    for j = 1:length(folders)

        a = sensitivities(2:5, 42, i);
        b = sensitivities(2:5, 42, j);

        if (a==b)
            angles_matrix(i, j)=0;
            continue
        end

        theta_rad1 = acos(dot(a, b) / (norm(a) * norm(b)));
        theta_rad2 = acos(dot(a, -b) / (norm(a) * norm(b)));

        theta_deg = rad2deg(min(theta_rad1, theta_rad2));

        angles_matrix(i, j) = theta_deg;
    end
end

%% Heatmap of angular distances
% Darker = smaller angle (more aligned sensitivities)

numColors = 256;
customGray = linspace(0.3, 1, numColors)' * [1 1 1];

figure('Color', 'w', 'Position', [100, 100, 350, 300]);
h = heatmap(round(angles_matrix), 'Colormap', customGray);
colorbar off;

saveas(gcf, 'Figure4_heatmap.svg');


%% -------------------------------------------------------------
%% 3) Compare observed angles to random expectation
%% -------------------------------------------------------------
% This section generates random normalized vectors and compares
% their angular distribution to the real sensitivities.

random_angles = [];

for i = 1:1000
    a = randn(5,1);
    b = randn(5,1);

    a=a/norm(a);
    b=b/norm(b);

    theta_rad1 = acos(dot(a, b));
    theta_rad2 = acos(dot(a, -b));

    theta_deg = rad2deg(min(theta_rad1, theta_rad2));

    random_angles=[random_angles, theta_deg];
end

%% Histogram comparison
% Overlapping histograms:
% - Real sensitivity angles
% - Random vector angles

figure('Color', 'w', 'Position', [100, 100, 400, 200]);

hold on;
histogram(angles_matrix(:), 'BinWidth', 10, ...
    'Normalization', 'probability', 'FaceAlpha', 0.5);

histogram(random_angles, 'BinWidth', 10, ...
    'Normalization', 'probability', 'FaceAlpha', 0.5);

hold off;

title('Histogram of Angles');
xlabel('Angle (degrees)');
ylabel('Frequency');
legend('Angles between sensitivities', 'Random Angles','Box', 'off');

ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial';
ax.LineWidth = 1.5;
ax.Box = 'off';
ax.TickDir = 'out';

saveas(gcf, 'Figure4_hist.svg');