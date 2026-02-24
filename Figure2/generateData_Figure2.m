clear
close all
format short e
beep off

%------------- Scenario Selection (Choose only one active) -------------%
scenario = "Gain_p";  % Options: "N_individuals", "Gain_p", "Gain_m"



switch scenario
    case "N_individuals"
        X_values = [512,1024,2048,4096];
        Gain_p = 1e-2;
        Gain_m = 1e-2;
        x_label = 'N_individuals';

    case "Gain_p"
        X_values = [1e-2,5e-2,1e-1,2e-1];
        N_individuals = 4096;
        Gain_m = 1e-2;
        x_label = 'Gain_p';

    case "Gain_m"
        X_values = [1e-2,2e-2,4e-2,8e-2];
        N_individuals = 4096;
        Gain_p = 1e-2;
        x_label = 'Gain_m';

    otherwise
        error('Unknown scenario selected.');
end

%------------------ Simulation Parameters ------------------%
ReSampling = 130;
TotalDevelopmentalTime = 50;
IntegrationStep = 0.01;
IntegrationSamplingTime = 0:IntegrationStep:TotalDevelopmentalTime;
NumberOfIntegrationSamples = length(IntegrationSamplingTime);
NNrepeticiones = 100;

%------------------ Initialize Storage ------------------%
store_ALL = {};
store_Error = [];

%%------------------ Calculate Nominal Trajectory ------------------%%
options = odeset('RelTol', 1e-6);
Xio = [0.6 0.8];

theta1 = zeros(NumberOfIntegrationSamples, 1);
theta2 = zeros(NumberOfIntegrationSamples, 1);
u = zeros(NumberOfIntegrationSamples, 1);

devNoise_x1 = zeros(NumberOfIntegrationSamples, 1);
devNoise_x2 = zeros(NumberOfIntegrationSamples, 1);

odeFunction = @(t, X) parsistemaggT(t, X, theta1, theta2, ...
    IntegrationStep, devNoise_x1, devNoise_x2, u);

[t, X] = ode45(odeFunction, IntegrationSamplingTime, Xio, options);
x1_N = X(:, 1)';
x2_N = X(:, 2)';

%%------------------ Calculate True Sensitivity ------------------%%
data = load('ABr.mat', 'Ar'); Ar = data.Ar;
data = load('ABr.mat', 'Br'); Br = data.Br;

Ar = eye(2) + IntegrationStep * Ar;
Br = IntegrationStep * Br;

TrueSensitivity_Theta1_out = zeros(2, NumberOfIntegrationSamples);
TrueSensitivity_Theta2_out = zeros(2, NumberOfIntegrationSamples);
TrueSensitivity_u_out = zeros(2, NumberOfIntegrationSamples);

for i = 2:NumberOfIntegrationSamples
    g1=x1_N(i);g2=x2_N(i);
    theta1=0;
    theta2=0;
    u=0;
    Arn=eval(Ar);Brn=eval(Br);

    TrueSensitivity_Theta1_out(:, i) = Arn * TrueSensitivity_Theta1_out(:, i - 1) + Brn(:, 1);
    TrueSensitivity_Theta2_out(:, i) = Arn * TrueSensitivity_Theta2_out(:, i - 1) + Brn(:, 2);
    TrueSensitivity_u_out(:, i)      = Arn * TrueSensitivity_u_out(:, i - 1)      + Brn(:, 3);
end

%%------------------ Main Simulation Loop ------------------%%
for value = X_values

    % Set scenario-specific parameters
    switch scenario
        case "N_individuals"
            N_individuals = value;
        case "Gain_p"
            Gain_p = value;
        case "Gain_m"
            Gain_m = value;
    end

    % Initialize error accumulator
    aa = [];

    parfor rep = 1:NNrepeticiones


        % Define developmental noise levels
        std_devNoise_x1 = Gain_p * x1_N(end);
        std_devNoise_x2 = Gain_p * x2_N(end);
        measurementNoise = Gain_m;

        % Genotype generation
        nLoci = 10;
        sigma = 1e-4;

        gScore_theta1 = randsample([-1, 0, 1], nLoci * N_individuals, true, [0.25, 0.5, 0.25]);
        gScore_theta1 = reshape(gScore_theta1, nLoci, N_individuals);
        gamma_theta1 = randn(nLoci, 1) * sqrt(sigma);
        loci_Theta1 = gScore_theta1 .* gamma_theta1;

        gScore_theta2 = randsample([-1, 0, 1], nLoci * N_individuals, true, [0.25, 0.5, 0.25]);
        gScore_theta2 = reshape(gScore_theta2, nLoci, N_individuals);
        gamma_theta2 = randn(nLoci, 1) * sqrt(sigma);
        loci_Theta2 = gScore_theta2 .* gamma_theta2;

        sum_loci_Theta1 = sum(loci_Theta1, 1);
        sum_loci_Theta2 = sum(loci_Theta2, 1);

        x1 = zeros(N_individuals, NumberOfIntegrationSamples);
        x2 = zeros(N_individuals, NumberOfIntegrationSamples);

        % Simulate each individual
        for individual = 1:N_individuals
            theta1 = ones(NumberOfIntegrationSamples, 1) * sum_loci_Theta1(individual);
            theta2 = ones(NumberOfIntegrationSamples, 1) * sum_loci_Theta2(individual);
            u = zeros(NumberOfIntegrationSamples, 1);

            devNoise_x1 = std_devNoise_x1 * randn(1, NumberOfIntegrationSamples);
            devNoise_x2 = std_devNoise_x2 * randn(1, NumberOfIntegrationSamples);

            odeFunction = @(t, X) parsistemaggT(t, X, theta1, theta2, ...
                IntegrationStep, devNoise_x1, devNoise_x2, u);

            [~, X] = ode45(odeFunction, IntegrationSamplingTime, Xio, options);
            x1(individual, :) = X(:, 1)';
            x2(individual, :) = X(:, 2)';
        end

        % Resample and apply measurement noise
        x1 = x1(:, 1:ReSampling:end);
        x2 = x2(:, 1:ReSampling:end);

        x1_measured = x1 .* (1 + measurementNoise * randn(size(x1)));
        x2_measured = x2 .* (1 + measurementNoise * randn(size(x2)));

        t = IntegrationSamplingTime(1:ReSampling:end);
        TotalSimulationTime = length(t);

        TrueSensitivity_Theta1 = TrueSensitivity_Theta1_out(:, 1:ReSampling:end);
        alpha_true = TrueSensitivity_Theta1(2, :) * gamma_theta1(1);

        % Estimate average effects via regression
        Z = [gScore_theta1' gScore_theta2' ones(N_individuals, 1)];

        alpha_static = zeros(size(Z, 2), TotalSimulationTime);
        for i = 1:TotalSimulationTime
            alpha_static(:, i) = Z \ x2_measured(:, i);
        end

        % Smooth estimates
        alpha_dynamic_linear = linearFilter(alpha_static(1, :), 10);
        alpha_dynamic_KF = adaptiveKF(alpha_static(1, :), 1, 10);

        % Compute % relative errors
        aa2 = [
            100 * norm(alpha_static(1, :)     - alpha_true) / norm(alpha_true), ...
            100 * norm(alpha_dynamic_linear       - alpha_true) / norm(alpha_true), ...
            100 * norm(alpha_dynamic_KF    - alpha_true) / norm(alpha_true)
        ];

        aa = [aa; aa2];
    end

    % Store aggregated statistics
    store_Error    = [store_Error;   value, median(aa, 1)]
    store_ALL{end+1} = aa;

end

%%------------------ Save Results Based on Scenario ------------------%%
switch scenario
    case "N_individuals"
        save N_individuals_Mean store_Error
        save('N_individuals_ALL.mat', 'store_ALL')

    case "Gain_p"
        save Gain_p_Mean store_Error
        save('Gain_p_ALL.mat', 'store_ALL')

    case "Gain_m"
        save Gain_m_Mean store_Error
        save('Gain_m_ALL.mat', 'store_ALL')
end
