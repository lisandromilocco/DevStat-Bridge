clear
close all
format short e
beep off
global theta1 theta2 IntegrationStep devNoise_x1 devNoise_x2

%------------------ Simulation Parameters ------------------%
TotalDevelopmentalTime = 50;          
IntegrationStep = 0.01;                       
IntegrationSamplingTime = 0:IntegrationStep:TotalDevelopmentalTime;
NumberOfIntegrationSamples = length(IntegrationSamplingTime);

N_individuals = 20;        % Number of individuals
ReSampling = 100;          % Downsampling factor

% Developmental noise (can be toggled on/off)
std_devNoise_x1 = 0 * 1e-2;     
std_devNoise_x2 = 0 * 1e-2;     
std_devNoise_x3 = 0 * 1e-2;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Simulate trait trajectories under genotype perturbation %%%%%%%%%%%%%%

rng(9999)  % Reproducibility

%--------------- Genotype: theta1 ----------------%
nLoci_theta1 = 1;                    
sigma_theta1 = 1e-4;                 
gScore_theta1 = randsample([-1, 0, 1], nLoci_theta1 * N_individuals, true, [0.25, 0.5, 0.25]);
gScore_theta1 = reshape(gScore_theta1, nLoci_theta1, N_individuals);
gamma_theta1 = randn(nLoci_theta1, 1) * sqrt(sigma_theta1);
loci_Theta1 = gScore_theta1 .* gamma_theta1;

%--------------- Genotype: theta2 ----------------%
nLoci_theta2 = 0;                    
sigma_theta2 = 1e-4;                 
gScore_theta2 = randsample([-1, 0, 1], nLoci_theta2 * N_individuals, true, [0.25, 0.5, 0.25]);
gScore_theta2 = reshape(gScore_theta2, nLoci_theta2, N_individuals);
gamma_theta2 = randn(nLoci_theta2, 1) * sqrt(sigma_theta2);
loci_Theta2 = gScore_theta2 .* gamma_theta2;

%---------------- Additive contributions ----------------%
sum_loci_Theta1 = sum(loci_Theta1, 1);
sum_loci_Theta2 = sum(loci_Theta2, 1);

% Initialize storage
x1 = []; x2 = [];
Delta_x1 = []; Delta_x2 = [];

options = odeset('RelTol', 1e-6);
Xio = [0.6 0.8];  % Initial state [g1, g2]

%--------------- Simulate developmental dynamics for each individual ---------------%
for individual = 1:N_individuals
    theta1 = ones(NumberOfIntegrationSamples, 1) * sum_loci_Theta1(individual);
    theta2 = ones(NumberOfIntegrationSamples, 1) * sum_loci_Theta2(individual);

    devNoise_x1 = std_devNoise_x1 * randn(1, NumberOfIntegrationSamples);
    devNoise_x2 = std_devNoise_x2 * randn(1, NumberOfIntegrationSamples);

    [t, X] = ode45(@sistemaggT, IntegrationSamplingTime, Xio, options);

    x1 = [x1; X(:, 1)'];  % Trait 1
    x2 = [x2; X(:, 2)'];  % Trait 2
end

%------------------- Population Reference -------------------%
x1_Ref = mean(x1)';
x2_Ref = mean(x2)';

Delta_x1 = x1 - x1_Ref';
Delta_x2 = x2 - x2_Ref';

%------------------- Sensitivity Calculation -------------------%
load ABr
Ar = eye(2) + IntegrationStep * Ar;
Br = IntegrationStep * Br;

[TotalSimulationTime, ~] = size(t);

TrueSensitivity_Theta1 = zeros(2, TotalSimulationTime);
TrueSensitivity_Theta2 = zeros(2, TotalSimulationTime);

theta1 = zeros(NumberOfIntegrationSamples, 1);
theta2 = zeros(NumberOfIntegrationSamples, 1);
u = zeros(NumberOfIntegrationSamples, 1);

devNoise_x1 = 0 * randn(1, NumberOfIntegrationSamples);
devNoise_x2 = 0 * randn(1, NumberOfIntegrationSamples);

[t, X] = ode45(@sistemaggT, IntegrationSamplingTime, Xio, options);

x1_N = X(:, 1)';
x2_N = X(:, 2)';

for i = 2:TotalSimulationTime
    g1 = x1_N(i); 
    g2 = x2_N(i); 

    theta1 = 0; 
    theta2 = 0;

    Arn = eval(Ar); 
    Brn = eval(Br);

    TrueSensitivity_Theta1(:, i) = Arn * TrueSensitivity_Theta1(:, i-1) + Brn(:, 1);
    TrueSensitivity_Theta2(:, i) = Arn * TrueSensitivity_Theta2(:, i-1) + Brn(:, 2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------ Resample Data ------------------------%
IntegrationStep = ReSampling * IntegrationStep;
IntegrationSamplingTime = IntegrationSamplingTime(1:ReSampling:end);
t = t(1:ReSampling:end);
[TotalSimulationTime, ~] = size(t);

x1 = x1(:, 1:ReSampling:end);
x2 = x2(:, 1:ReSampling:end);
x1_Ref = x1_Ref(1:ReSampling:end);
x2_Ref = x2_Ref(1:ReSampling:end);
Delta_x1 = Delta_x1(:, 1:ReSampling:end);
Delta_x2 = Delta_x2(:, 1:ReSampling:end);

TrueSensitivity_Theta1 = TrueSensitivity_Theta1(:, 1:ReSampling:end);
TrueSensitivity_Theta2 = TrueSensitivity_Theta2(:, 1:ReSampling:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------ Average Locus Effect ------------------%
[rows, columns] = size(x1);
E1 = []; 
E2 = [];

for i = 1:TotalSimulationTime
    Z = [gScore_theta1' gScore_theta2' ones(rows, 1)];  % Genotype matrix

    average_effects_x1 = Z \ x1(:, i);
    average_effects_x2 = Z \ x2(:, i);

    E1 = [E1, average_effects_x1];
    E2 = [E2, average_effects_x2];
end

%------------------ Plot: Sensitivity vs Additive Effect ------------------%
colors = lines(5);  % Color palette

figure('Color', 'w', 'Position', [100, 100, 600, 220]);

plot(TrueSensitivity_Theta1(1, :) * gamma_theta1(1), ...
    'Color', colors(5, :), 'LineWidth', 3.5, ...
    'DisplayName', '$s_\lambda(t)\gamma$', 'MarkerSize', 6);
hold on;

scatter(1:length(E1(1, :)), E1(1, :), 90, ...
    'MarkerEdgeColor', colors(1, :), 'MarkerFaceColor', colors(1, :), ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.6, ...
    'DisplayName', '$\alpha(t)$');

ax = gca;
ax.FontSize = 18;
ax.FontName = 'Arial';
ax.LineWidth = 2;
ax.Box = 'off';
ax.TickDir = 'out';

xlabel('Developmental time', 'FontSize', 18, 'FontName', 'Arial');
ylabel('Phenotypic effect', 'FontSize', 18, 'FontName', 'Arial');
xlim([0 50]);

legend('Interpreter', 'latex', 'Location', 'southeast', ...
    'FontSize', 18, 'Box', 'off');

ax.YAxis.Exponent = -2;
ytickformat('%.1f')

% Tight layout
outerpos = ax.OuterPosition;
ti = ax.TightInset;
ax.Position = [outerpos(1) + ti(1), outerpos(2) + ti(2), ...
    outerpos(3) - ti(1) - ti(3), outerpos(4) - ti(2) - ti(4) - 0.05];

