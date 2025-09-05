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
response=zeros(length(ntrials),1);

% Loop over different durations
for i = 1:length(durations)
    duration = durations(i);
    num_steps = duration / 100;
    expectA = 2 * c; %output value
    expectB = -2 * c; %output value
    %responseA=sum(normrnd(expect A, variance, [1:ntrials], timesteps))+noise
    %noise:rand(0,1)
    % Simulate the trials
    decisions = zeros(ntrials, 1);
    actual = rand(ntrials, 1);

    for trial = 1:ntrials %randomly generate A or B
            if actual(trial) >0.5
                evidence=sum(normrnd(expectA, var(expectA), [1:ntrials]))+noise; %responseA
            else
                evidence=sum(normrnd(expectB, var(expectB), [1:ntrials]))+noise; %responseB
            end
        if evidence > c
            response(trial)=1;
        else
            response(trial)=-1;
        end
    end
    
    hits = sum(response == 1) / ntrials;
    pcorrect=hits*100;
    fas = sum(response == -1) / ntrials;
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