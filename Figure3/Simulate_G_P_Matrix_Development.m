clear
close all
format short e
beep off

%%%%%% Adjustment parameters
TotalDevelopmentalTime = 50;                      % total simulation time
IntegrationStep = .01;                     % integration sampling period
IntegrationSamplingTime=[0:IntegrationStep:TotalDevelopmentalTime];
NumberOfIntegrationSamples=length(IntegrationSamplingTime);  % simulation time


N_individuals = 5000;         % Number of simulated individuals

measurementNoise = 0*.02*1e-2;                 % measurement noise
ReSampling=100;                    % resampling

std_devNoise_x1 = 0*1*1e-2;              % noise in states
std_devNoise_x2 = 0*1*1e-2;             % noise in states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Simulation of temporal responses to genotype changes %%%%%%%%
rng(9)

col2=[ 0]';
col2= 0.5 ./ 2.^col2;

round(col2, 2, "significant")
round(col2, 1, "significant")

col1=zeros(length(col2),1)+0.5

mu_values = [ col1 col2 ] 


P_store = zeros(2, 2, size(mu_values, 1));
G_store = zeros(2, 2, size(mu_values, 1));


% Prepare the figure for 2x4 grid
figure;

% Loop through the 6 (mu_p_theta1, mu_p_theta2) pairs
for i = 1:size(mu_values, 1)
%for i = 1:length(mu_values)

    % Extract mu_p_theta1 and mu_p_theta2 for the current iteration
    mu_p_theta1 = mu_values(i, 1);
    mu_p_theta2 = mu_values(i, 2);

    % Parametro theta1
    nLoci_theta1 = 10; % Bumber of loci
    
    % Initialize gScore_theta1 matrix
    gScore_theta1 = zeros(nLoci_theta1, N_individuals);

    % Loop through loci
    for locus = 1:nLoci_theta1
        p_i = mu_p_theta1;
        q_i = 1-p_i;

        % Sample genotypes for all individuals at this locus
        gScore_theta1(locus, :) = randsample([-1, 0, 1], ...
            N_individuals, true, [q_i^2, 2*p_i*q_i, p_i^2]);
    end

    % Now we define the gamma values which determine the developmental
    % parameters

    sigma_theta1 = 1e-4; % varianza del coeficiente
    gamma_theta1 = randn(nLoci_theta1, 1)*(sigma_theta1);% * 5; %I suggest multiplying by 4 because
    loci_Theta1 = gScore_theta1 .* gamma_theta1;

    % Second parameter
    nLoci_theta2 = 10; % Número de locis

    % Initialize gScore_theta2 matrix
    gScore_theta2 = zeros(nLoci_theta2, N_individuals);

    % Loop through loci
    for locus = 1:nLoci_theta2
        p_i = mu_p_theta2;
        q_i = 1-p_i;

        % Sample genotypes for all individuals at this locus
        gScore_theta2(locus, :) = randsample([-1, 0, 1], ...
            N_individuals, true, [q_i^2, 2*p_i*q_i, p_i^2]);
    end

    % Now we define the gamma values which determine the developmental
    % parameters
    sigma_theta2 = 1e-4; % varianza del coeficiente
    gamma_theta2 = randn(nLoci_theta2, 1)*sigma_theta2;
    loci_Theta2 = gScore_theta2 .* gamma_theta2;


    % Sum the (additive) contributions of all loci to the parameter
    % values of each individual
    sum_loci_Theta1 = sum(loci_Theta1,1);
    sum_loci_Theta2 = sum(loci_Theta2,1);


    % Introduce the environmental variable
    sigma_e = 1.5e-3;
    environment = normrnd(0, sigma_e, [1,N_individuals]);

    x1=[];x2=[];
    Delta_x1=[];Delta_x2=[];

    options = odeset('RelTol',1e-6);
    Xio=[0.6 0.8];   % Initial condition for the simulations

    % Loop to simulate developmental trajectories

    % Ensure Parallel Computing Toolbox is available
    %parpool;

    x1 = zeros(N_individuals, NumberOfIntegrationSamples); % Preallocate memory
    x2 = zeros(N_individuals, NumberOfIntegrationSamples);


    parfor individual = 1:N_individuals
        %for individual = 1:N_individuals

        theta1 = ones(NumberOfIntegrationSamples, 1) * sum_loci_Theta1(individual);
        theta2 = ones(NumberOfIntegrationSamples, 1) * sum_loci_Theta2(individual);
        u = ones(NumberOfIntegrationSamples, 1) * environment(individual);

        % Simulate results
        devNoise_x1 = std_devNoise_x1 * randn(1, NumberOfIntegrationSamples);
        devNoise_x2 = std_devNoise_x2 * randn(1, NumberOfIntegrationSamples);


        odeFunction = @(t, X) parsistemaggT(t, X, theta1, theta2, ...
            IntegrationStep, devNoise_x1, devNoise_x2, u);

        [t, X] = ode45(odeFunction, IntegrationSamplingTime, Xio, options);
        % Store results
        x1(individual, :) = X(:, 1)'; % Store trait 1
        x2(individual, :) = X(:, 2)'; % Store trait 2
    end


    % Keep the phenotypes at SOME POINT IN DEVELOPMENT
    x1_ss = x1(:,1000);
    x2_ss = x2(:,1000);
    [rows,columns]=size(x1_ss);


    % Add measurement noise, optional
    x1_ss = x1_ss + measurementNoise * mean(x1_ss) .* randn(size(x1_ss));
    x2_ss = x2_ss + measurementNoise * mean(x2_ss) .* randn(size(x2_ss));

    x1_ss = x1_ss - mean(x1_ss);
    x2_ss = x2_ss - mean(x2_ss);

    % Now we obtain the average effect of each loci at steady state
    Z=[gScore_theta1' gScore_theta2' ones(rows,1)]; % rows = n_individuals
    average_effects_x1 = Z \ x1_ss;
    average_effects_x2 = Z \ x2_ss;

    % Calculate the G matrix as sum of the variance of each locus

    average_effects = [ average_effects_x1 average_effects_x2]; % Last effect is the mean

    genotype= [ gScore_theta1' gScore_theta2' ];

    G=zeros(2,2);
    for loci = 1:(length(average_effects_x1)-1) %Last one is the mean!

        n_A2 = sum(genotype(:,loci) == 1) * 2 + sum(genotype (:,loci)==0 );
        n_A1 = sum(genotype(:,loci) == -1) * 2 + sum(genotype (:,loci)==0 );
        p_i=n_A2/(n_A1+n_A2);
        q_i=n_A1/(n_A1+n_A2);

        %p_i=0.5;
        %q_i=0.5;

        G_i = average_effects(loci,:)' * average_effects(loci,:) * 2 * p_i * q_i;
        G=G+G_i;
    end

    BV=zeros(2,N_individuals);
    for ind = 1:N_individuals
        addup=0;
        for loci = 1:length(average_effects_x1)-1
            addup=addup+genotype(ind,loci)*average_effects(loci,:);
        end
        BV(:,ind)=addup+average_effects(end,:);
    end


    x_ss = [x1_ss x2_ss];
    P = cov(x_ss);
 
    % Store in 3D arrays
    P_store(:,:,i) = P;
    G_store(:,:,i) = G;

norm(P)

end

% Save to file
save('matrices_10.mat', 'P_store', 'G_store');
