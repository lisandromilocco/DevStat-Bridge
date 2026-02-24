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

col2=[9 8 7 6 5 4 3 2 1 0]';
col2= 0.5 ./ 2.^col2;


col1=zeros(length(col2),1)+0.5

mu_values = [ col1 col2 ]


P_store = zeros(2, 2, size(mu_values, 1));
G_store = zeros(2, 2, size(mu_values, 1));


%%% NEW CODE: store relative errors
n_replays = 50;   % how many times to replay reproduction per parameter set

observed_deltaZ_store = zeros(size(mu_values, 1), n_replays, 2);
deltaZ_Lande_store    = zeros(size(mu_values, 1), 2);
deltaZ_naive_store    = zeros(size(mu_values, 1), 2);



% Prepare the figure for 2x4 grid
figure;

% Loop through the 6 (mu_p_theta1, mu_p_theta2) pairs
for i = 1:size(mu_values, 1)
    %for i = 1:length(mu_values)

    % Extract mu_p_theta1 and mu_p_theta2 for the current iteration
    mu_p_theta1 = mu_values(i, 1);
    mu_p_theta2 = mu_values(i, 2);

    % Parametro theta1
    nLoci_theta1 = 10; % Number of loci

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


    % Keep the phenotypes at the end of development
    x1_ss = x1(:,end);
    x2_ss = x2(:,end);
    [rows,columns]=size(x1_ss);

    dx1_ss = x1_ss - mean(x1_ss);
    dx2_ss = x2_ss - mean(x2_ss);

    % Now we obtain the average effect of each loci at steady state
    Z=[gScore_theta1' gScore_theta2' ones(rows,1)]; % rows = n_individuals
    average_effects_x1 = Z \ dx1_ss;
    average_effects_x2 = Z \ dx2_ss;

    % Calculate the G matrix as sum of the variance of each locus

    average_effects = [ average_effects_x1 average_effects_x2]; % Last effect is the mean

    genotype= [ gScore_theta1' gScore_theta2' ];

    G=zeros(2,2);
    for loci = 1:(length(average_effects_x1)-1) %Last one is the mean!

        n_A2 = sum(genotype(:,loci) == 1) * 2 + sum(genotype (:,loci)==0 );
        n_A1 = sum(genotype(:,loci) == -1) * 2 + sum(genotype (:,loci)==0 );
        p_i=n_A2/(n_A1+n_A2);
        q_i=n_A1/(n_A1+n_A2);

        G_i = average_effects(loci,:)' * average_effects(loci,:) * 2 * p_i * q_i;
        G=G+G_i;
    end

    BV=zeros(2,N_individuals);
    for indi = 1:N_individuals
        addup=0;
        for loci = 1:length(average_effects_x1)-1
            addup=addup+genotype(indi,loci)*average_effects(loci,:);
        end
        BV(:,indi)=addup+average_effects(end,:);
    end


    dx_ss = [dx1_ss dx2_ss];
    P = cov(dx_ss);

    % Store in 3D arrays
    P_store(:,:,i) = P;
    G_store(:,:,i) = G;


    %%% One round of selection and reproduction %%%

    % Define optimum
    optimum = [4, 4];

    % Calculate Euclidean distance to optimum
    distances = sqrt((x1_ss - optimum(1)).^2 + (x2_ss - optimum(2)).^2);

    % Select the best 50% individuals
    [~, idx_sorted] = sort(distances);
    n_selected = floor(N_individuals/2);
    selected_idx = idx_sorted(1:n_selected);

    % Phenotypes of selected and full population
    mean_all = mean([x1_ss x2_ss]);
    mean_selected = mean([x1_ss(selected_idx) x2_ss(selected_idx)]);

    % Selection differential
    s = mean_selected - mean_all;

    % Lande’s prediction
    deltaZ_Lande = (G / P) * s';   % equivalent to G*inv(P)*s
    % Naive prediction (G=P)
    deltaZ_naive = s;


    %STORE PREDICTIONS

    deltaZ_Lande_store(i, :)    = deltaZ_Lande(:);
    deltaZ_naive_store(i, :)    = deltaZ_naive(:);



    %%% Reproduction with recombination %%%
    % Collect parental genotypes
    parents_theta1 = gScore_theta1(:,selected_idx);
    parents_theta2 = gScore_theta2(:,selected_idx);

    for rep = 1:n_replays

        %%% Random pairing of parents
        perm = randperm(n_selected);
        pairs = reshape(perm,2,[]);  % each column = a couple

        n_offspring = n_selected * 2; % 4 offspring per pair
        offspring_theta1 = zeros(nLoci_theta1, n_offspring);
        offspring_theta2 = zeros(nLoci_theta2, n_offspring);

        for c = 1:size(pairs,2)
            parentA1 = parents_theta1(:,pairs(1,c));
            parentB1 = parents_theta1(:,pairs(2,c));
            parentA2 = parents_theta2(:,pairs(1,c));
            parentB2 = parents_theta2(:,pairs(2,c));

            for kid = 1:4
                child_index = (c-1)*4 + kid;

                % Recombine loci for theta1
                for locus = 1:nLoci_theta1
                    offspring_theta1(locus, child_index) = recombine(parentA1(locus), parentB1(locus));
                end
                % Recombine loci for theta2
                for locus = 1:nLoci_theta2
                    offspring_theta2(locus, child_index) = recombine(parentA2(locus), parentB2(locus));
                end
            end
        end

        %%% Develop offspring phenotypes
        offspring_loci_theta1 = offspring_theta1 .* gamma_theta1;
        offspring_loci_theta2 = offspring_theta2 .* gamma_theta2;

        sum_off_Theta1 = sum(offspring_loci_theta1,1);
        sum_off_Theta2 = sum(offspring_loci_theta2,1);

        x1_off = zeros(n_offspring, NumberOfIntegrationSamples);
        x2_off = zeros(n_offspring, NumberOfIntegrationSamples);

        parfor ind = 1:n_offspring
            theta1 = ones(NumberOfIntegrationSamples, 1) * sum_off_Theta1(ind);
            theta2 = ones(NumberOfIntegrationSamples, 1) * sum_off_Theta2(ind);
            u = ones(NumberOfIntegrationSamples, 1) * normrnd(0, sigma_e);

            devNoise_x1 = std_devNoise_x1 * randn(1, NumberOfIntegrationSamples);
            devNoise_x2 = std_devNoise_x2 * randn(1, NumberOfIntegrationSamples);

            odeFunction = @(t, X) parsistemaggT(t, X, theta1, theta2, ...
                IntegrationStep, devNoise_x1, devNoise_x2, u);

            [t, X] = ode45(odeFunction, IntegrationSamplingTime, Xio, options);
            x1_off(ind,:) = X(:,1)';
            x2_off(ind,:) = X(:,2)';
        end

        x1_off_ss = x1_off(:,end);
        x2_off_ss = x2_off(:,end);
        mean_offspring = mean([x1_off_ss x2_off_ss]);

        %%% Compare observed vs predicted response
        observed_deltaZ = mean_offspring - mean_all;


        % Store the vectors (transpose if needed to be row/column consistent)
        observed_deltaZ_store(i, rep, :) = observed_deltaZ(:);


    end


    
    fprintf('Iteration %d:\n', i);
    fprintf('Observed ΔZ = [%f, %f]\n', observed_deltaZ);
    fprintf('Predicted (Lande) ΔZ = [%f, %f]\n', deltaZ_Lande);
    fprintf('Predicted (naive) ΔZ = [%f, %f]\n\n', deltaZ_naive);

end


% Optionally save
save('predictions.mat','observed_deltaZ_store','deltaZ_Lande_store','deltaZ_naive_store', ...
    'mu_values', 'n_replays');



%%% Recombination function for reproduction %%%
function child = recombine(allele1, allele2)
if allele1==1 && allele2==1
    child = 1;
elseif allele1==-1 && allele2==-1
    child = -1;
elseif (allele1==0 && allele2==-1) || (allele1==-1 && allele2==0)
    child = randsample([-1,0],1);  % 50% aa, 50% aA
elseif (allele1==0 && allele2==1) || (allele1==1 && allele2==0)
    child = randsample([0,1],1);   % 50% aA, 50% AA
elseif (allele1==-1 && allele2==1) || (allele1==1 && allele2==-1)
    child = 0;   % 100% aA
elseif allele1==0 && allele2==0
    child = randsample([-1,0,1],1,true,[0.25,0.5,0.25]);
else
    error('Unexpected allelic combination');
end
end

