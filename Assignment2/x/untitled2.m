%% Perception Assignment 2
% Name: Rachel Chen 

%% Question 1
% Simulate a two-alternative, forced-choice experiment. Each trial consists of one
% interval with stimulus A or stimulus B. In each 100 ms, the observer collects evidence for
% whether the stimulus is A or B. The evidence in each 100 ms time step consists of a
% single number. The expected value of that number is +2 * c for an A stimulus, and -2 * c
% for a B stimulus, where c is the stimulus contrast. However, that number is noisy,
% perturbed by independent random draws each time step from a standard Gaussian
% distribution (mean zero, SD = 1).

%% (a.1) Simulate this experiment for stimulus durations ranging from 100 ms to 2 s and a
% contrast of 0.2. Assume A and B are equally likely, payoffs are symmetric and an
% optimal criterion. Thus, for example, for a 500 ms stimulus there are 5 time steps.
% And, if the stimulus has contrast 0.2, that means each time step is a random draw
% from a Gaussian with mean 0.4 and variance 1.0. The observer will optimally
% combine those bits of evidence by summing the 5 numbers from the 5 time steps.
% Thus, for an A stimulus, the total evidence is expected to be 5*.4 = 2, on average,
% and -2 for a B stimulus. What is the optimal criterion? Plot your results as
% percentage correct as a function of stimulus duration. Replot the results as d’ as a
% function of duration. Compare your simulation results with theoretical predictions.

% Setting up parameters
c = 0.2; %contrast level
durations = 100:100:2000;
ntrials = 10000;
noise_level = rand(0,1); % Change this to vary noise

%Empty space for storing
pcorrect = zeros(length(durations), 1); %percentage corrctness
d_prime = zeros(length(durations), 1);
response=zeros(:,2);

% Loop over different durations
for i = 1:length(durations)
    duration = durations(i);
    num_steps = duration / 100;
    expectA = 2 * c; %output value
    ecpectB = -2 * c; %output value
    %responseA=sum(normrnd(expect A, variance, [1:ntrials], timesteps))+noise
    %noise:rand(0,1)
    % Simulate the trials
    decisions = zeros(ntrials, 1);
    actual = rand(ntrials, 1) > 0.5;

    for trial = 1:ntrials %randomly generate A or B
        evidence = 0;
        for step = 1:num_steps
            if actual(trial) == 1
                evidence = evidence + normrnd(evidence_A, noise_level);
            else
                evidence = evidence + normrnd(evidence_B, noise_level);
            end
        end
        decisions(trial) = evidence > 0;
    end
    
    correct = sum(decisions == actual);
    pcorrect(i) = (correct / ntrials) * 100;
    
    hits = sum(decisions(actual == 1) == 1) / sum(actual == 1);
    fas = sum(decisions(actual == 0) == 1) / sum(actual == 0);
    d_prime(i) = norminv(hits) - norminv(fas);
end

figure;
subplot(2,1,1);
plot(durations, pcorrect);
xlabel('Stimulus Duration (ms)');
ylabel('Percent Correct');
title('Percentage correct as a function of stimulus duration');

subplot(2,1,2);
plot(durations, d_prime);
xlabel('Stimulus Duration (ms)');
ylabel('d-prime');
title('dprime as a function of stimulus duration');

% For theoretical prediction, the longer a target is presented, the easier
% people can detect it with more details. When there are two targets, the
% longer they are presnted, the easier people can distinguish between them.
%Comparing to the results, percentage of correctness and d-prime increase
%as the duration of stimulus increase. Thus, our simulation supported the
%theoretical prediction. 

%responseA(meanA)=2xcotrastctimestep
%[mean A-meanB]/根号sigma
%% (a.2) Now, repeat the above (determine optimal criterion, simulate, plot results, 
% compare to theoretical predictions) for an experiment in which the duration is fixed 
% at 0.2 s and contrast ranges from 0.05 to 1.

% Parameters
contrast_levels = 0.05:0.05:1; % varying contrasts from 0.05 to 1
duration = 200; % fixed duration in ms
num_trials = 1000;
noise_level = 0.7; % Noise level (can be varied as needed)

% Initialize arrays to store results
percent_correct = zeros(length(contrast_levels), 1);
d_prime = zeros(length(contrast_levels), 1);

% Loop over different contrast levels
for i = 1:length(contrast_levels)
    contrast = contrast_levels(i);
    num_steps = duration / 100; % number of 100 ms steps
    evidence_A = 2 * contrast;  % mean evidence for A based on current contrast
    evidence_B = -2 * contrast; % mean evidence for B based on current contrast

    % Simulate the trials
    decisions = zeros(num_trials, 1);
    actual = rand(num_trials, 1) > 0.5; % actual conditions, 1 for A, 0 for B
    for trial = 1:num_trials
        evidence = 0;
        for step = 1:num_steps
            if actual(trial) == 1 % stimulus A
                evidence = evidence + normrnd(evidence_A, noise_level);
            else % stimulus B
                evidence = evidence + normrnd(evidence_B, noise_level);
            end
        end
        decisions(trial) = evidence > 0;
    end
    
    % Analyze performance
    correct = sum(decisions == actual);
    percent_correct(i) = (correct / num_trials) * 100;
    
    % Calculate d'
    hits = sum(decisions(actual == 1) == 1) / sum(actual == 1);
    fas = sum(decisions(actual == 0) == 1) / sum(actual == 0);
    d_prime(i) = norminv(hits) - norminv(fas);
end

% Plot percent correct
figure;
subplot(1,2,1)
plot(contrast_levels, percent_correct);
xlabel('Stimulus Contrast');
ylabel('Percent Correct');
title('Performance Across Different Contrasts');

% Plot d'
subplot(1,2,2)
plot(contrast_levels, d_prime);
xlabel('Stimulus Contrast');
ylabel('d''');
title('d'' Across Different Contrasts');

%% (b) Change the probability of an A stimulus to 0.75. Use a single fixed duration of
% 0.4 s and contrast of 0.1. What is the optimal criterion now? Simulate performance
% with that optimal criterion, and then simulate it with the criterion you used in part (a),
% comparing the performance across the two for 100 trials each. You should be able to
% compute the optimal criterion. However, you can also determine/estimate that
% criterion by doing large-scale simulations for a range of criteria (I’d prefer the closedform
% computation if you can manage it).

% Parameters
contrast = 0.1;
duration = 400; % 400ms
num_trials = 100;
noise_level = 1;
probability_A = 0.75;
probability_B = 1 - probability_A;
evidence_A = 2 * contrast;
evidence_B = -2 * contrast;

% Compute the optimal criterion
mu_diff = evidence_A - evidence_B;
optimal_criterion = log(probability_A / probability_B) * (noise_level^2 / mu_diff) + (evidence_A^2 - evidence_B^2) / (2 * mu_diff);

% Criterion from part (a) (we assumed it to be 0)
criterion_a = 0;

% Simulate trials for both criteria
results_optimal = simulate_trials(num_trials, evidence_A, evidence_B, noise_level, optimal_criterion);
results_a = simulate_trials(num_trials, evidence_A, evidence_B, noise_level, criterion_a);

% Compute performance for both criteria
percent_correct_optimal = sum(results_optimal == (rand(num_trials, 1) < probability_A)) / num_trials * 100;
percent_correct_a = sum(results_a == (rand(num_trials, 1) < probability_A)) / num_trials * 100;

% Display results
fprintf('Performance with optimal criterion: %.2f%%\n', percent_correct_optimal);
fprintf('Performance with criterion from part (a): %.2f%%\n', percent_correct_a);

%% (2) With the same setup as in (1), switch to a reaction-time experiment using the
% accumulated evidence values in a drift-diffusion framework. That is, accumulate the sum
% across time steps until the sum hits a bound representing a decision to respond “A”
% (with bound value +b) and another for response “B” (with value -b). Use simulations to 
% characterize the reaction-time distributions for correct vs. error trials. Do this for
% contrasts of 0.1 and 0.5. For each contrast, run simulations for a near decision
% boundary (relatively small value of b), a distant decision boundary (large value of b),
% and an asymmetric pair of decision boundaries (+b and -d). What happens to the hit and
% false-alarm rates (treating stimulus B as “noise” and A as “signal”) and RT distributions
% with each manipulation of the model?

% Define parameters
contrasts = [0.1, 0.5];
boundaryConfigs = {[1, 1], [5, 5], [1, 3]};  % [b, d] for each configuration
nTrials = 1000;
maxSteps = 100;  % Maximum number of time steps to simulate

% Loop over contrasts and boundary configurations
for c = 1:length(contrasts)
    contrast = contrasts(c);
    for b = 1:length(boundaryConfigs)
        boundaries = boundaryConfigs{b};
        hitRate = 0;
        falseAlarmRate = 0;
        RTsCorrect = [];
        RTsError = [];
        for t = 1:nTrials
            stimulus = randi([0, 1]);  % 0 for A, 1 for B
            evidence = 0;
            step = 0;
            while abs(evidence) < max(boundaries) && step < maxSteps
                step = step + 1;
                meanValue = (1-2*stimulus) * 2 * contrast;  % Mean value of evidence
                evidence = evidence + normrnd(meanValue, 1);  % Accumulate evidence

                % Check if decision boundaries are reached
                if evidence >= boundaries(1)
                    if stimulus == 0
                        hitRate = hitRate + 1;
                        RTsCorrect = [RTsCorrect, step];
                    else
                        RTsError = [RTsError, step];
                    end
                    break;
                elseif evidence <= -boundaries(2)
                    if stimulus == 1
                        falseAlarmRate = falseAlarmRate + 1;
                        RTsCorrect = [RTsCorrect, step];
                    else
                        RTsError = [RTsError, step];
                    end
                    break;
                end
            end
        end
        hitRate = hitRate / (nTrials / 2) * 100;  % Adjust for number of A trials
        falseAlarmRate = falseAlarmRate / (nTrials / 2) * 100;  % Adjust for number of B trials

        % Plot RT distributions
        figure;
        histogram(RTsCorrect, 'Normalization', 'probability');
        hold on;
        histogram(RTsError, 'Normalization', 'probability');
        title(['Contrast: ', num2str(contrast), ', Boundaries: ', num2str(boundaries)]);
        legend('Correct', 'Error');
        xlabel('Reaction Time (steps)');
        ylabel('Probability');
        hold off;
    end
end


%% functions
function decisions = simulate_trials(num_trials, evidence_A, evidence_B, noise_level, criterion)
    decisions = zeros(num_trials, 1);
    for trial = 1:num_trials
        actual = rand > 0.75; % 0 for A, 1 for B
        if actual == 0 % stimulus A
            evidence = normrnd(evidence_A, noise_level);
        else % stimulus B
            evidence = normrnd(evidence_B, noise_level);
        end
        decisions(trial) = evidence > criterion; 
    end
end
