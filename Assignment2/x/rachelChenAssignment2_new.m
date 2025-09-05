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

% (a.1) Simulate this experiment for stimulus durations ranging from 100 ms to 2 s and a
% contrast of 0.2. Assume A and B are equally likely, payoffs are symmetric and an
% optimal criterion. Thus, for example, for a 500 ms stimulus there are 5 time steps.
% And, if the stimulus has contrast 0.2, that means each time step is a random draw
% from a Gaussian with mean 0.4 and variance 1.0. The observer will optimally
% combine those bits of evidence by summing the 5 numbers from the 5 time steps.
% Thus, for an A stimulus, the total evidence is expected to be 5*.4 = 2, on average,
% and -2 for a B stimulus. What is the optimal criterion? Plot your results as
% percentage correct as a function of stimulus duration. Replot the results as d’ as a
% function of duration. Compare your simulation results with theoretical predictions.

% Given parameters
contrast = 0.2;
durations = 100:100:2000; % From 100 ms to 2 s

% Preallocate arrays for storing results
pcorrect = zeros(size(durations));
dprime = zeros(size(durations));

% Loop over different stimulus durations
for i = 1:length(durations)
    duration = durations(i);
    nSteps = duration / 100; % Number of time steps
    nTrials = 1000; % Number of trials
    
    % Simulate evidence collection for each trial
    totalEvidence = zeros(nTrials, 1);
    for t = 1:nTrials
        for step = 1:nSteps
            % Generate evidence based on stimulus and noise
            noise = randn; % Random draw from standard Gaussian
            % Assuming stimulus A is presented in all trials
            evidence = 2 * contrast + noise; 
            totalEvidence(t) = totalEvidence(t) + evidence;
        end
    end
    
    % Determine optimal criterion (assuming 0 as optimal)
    criterion = 0;
    
    % Determine responses based on total evidence and criterion
    % Assume choosing A if total evidence > criterion, B otherwise
    responses = totalEvidence > criterion;
    
    % Simulate the actual stimulus labels (0 for B, 1 for A)
    stimulusLabels = rand(nTrials, 1) > 0.5;
    
    % Calculate percentage correct
    pcorrect(i) = sum(responses == stimulusLabels) / nTrials;
    
    % Calculate hit rate and false alarm rate for d' calculation
    hits = sum(stimulusLabels == 1 & responses == 1);
    falseAlarms = sum(stimulusLabels == 0 & responses == 1);
    hitRate = hits / sum(stimulusLabels == 1);
    falseAlarmRate = falseAlarms / sum(stimulusLabels == 0);
    
    % Avoid division by zero or log(0) in d' calculation
    hitRate = max(hitRate, 1/nTrials);
    falseAlarmRate = max(falseAlarmRate, 1/nTrials);
    
    % Calculate d'
    dprime(i) = norminv(hitRate) - norminv(falseAlarmRate);
end

% Plot percentage correct
figure;
plot(durations, pcorrect * 100);
xlabel('Stimulus Duration (ms)');
ylabel('Percentage Correct');
title('Percentage Correct vs Stimulus Duration');

% Plot d'
figure;
plot(durations, dprime);
xlabel('Stimulus Duration (ms)');
ylabel('d''');
title('d'' vs Stimulus Duration');



























































%% Setting up the experiment

% Define parameters
ntrials=1000;
stimDur=100:100:2000; %100ms to 2s with 100ms for each step
c=0.2; % stimulus contrast, manually adjustable

%Initializatioin
correct=zeros(size(stimDur)); %making a blank space
dprime=zeros(size(stimDur)); %making a blank space

%making the experiment for varying duration
for i=1:length(stimDur);
    numTimesteps=round(stimDur(i)/100);
    decision=zeros(ntrials,1);  %making a blank space
    dispStimulusA=zeros(ntrials,1); %making a blank space to record
    dispStimulusB=zeros(ntrials,1);

    for trial=1:ntrials
        %rand generate stimulus A or B
        if rand()<0.5
            stimuls='A stimulus';
        else
            stimuls='B stimulus';
        end

        for step=1:numTimesteps
            noise=randn; %randomly assign noise 
            if stimuls=='A stimulus';
                evidence=2*c+noise;
            else 
                evidence=-2*c+noise;
            end 
            if stimulus=='A stimulus'
                dispStimulusA(trial)=dispStimulusA(trial)+evidence;
            else
                dispStimulusB(trial)=dispStimulusB(trial)+evidence;
            end 
        end
    end 

%Determine optimal criterion (setting 0 as optimal in this case)
criterion=0;

%Calculate percentage correctness
correctA=sum(dispStimulusA>criterion);
correctB = sum(dispStimulusB < criterion);
correct(i)=(correctA+correctB)/(2*ntrials);

%Calculate d-prime
meanA=mean (dispStimulusA);
varA=var(dispStimulusA);
meanB=mean(dispStimulusB);
varB=var(dispStimulusB);
dprime(i)=(meanA-meanB)/sqrt((varA+varB)/2);









end

%Plotting
figure1=figure;
subplot(2,1,1);
plot(stimDur, correct*100,'-o');
xlabel('Stimulus Duration (s)');
ylabel('Percentage of Correctness (%)');
title('Percentage of correctness as a function of stimulus duration');

subplot(2,1,2);
plot(stimDur, dprime,'-o');
xlabel('Stimulus Duration (s)');
ylabel('d-prime');
title('d-prime as a function of stimulus duration');

% For theoretical prediction, the longer a target is presented, the easier
% people can detect it with more details. When there are two targets, the
% longer they are presnted, the easier people can distinguish between them.
%Comparing to the results, percentage of correctness and d-prime increase
%as the duration of stimulus increase. Thus, our simulation supported the
%theoretical prediction. 

%% (a.2)Now, repeat the above (determine optimal criterion, simulate, plot results, compare
% to theoretical predictions) for an experiment in which the duration is fixed at 0.2 s
% and contrast ranges from 0.05 to 1

stimDur=0.2;
contrast=0.05:0.01:1;
ntrials=1000;

%Initializatioin
correct2=zeros(size(contrast)); %making a blank space
dprime2=zeros(size(contrast)); %making a blank space


%making the experiment for varying duration
for i=1:length(contrast);
    c=contrast(i);
    numTimesteps=round(stimDur/0.1);
    decision=zeros(ntrials,1);  %making a blank space
    dispStimulus=zeros(ntrials,1); %making a blank space to record
    for trial=1:ntrials
        if rand()<0.5
            stimuls='A stimulus';
            dispStimulus(trial)=1; %set A as 1
            outputValue=2*c;
        else
            stimuls='B stimulus';
            dispStimulus(trial)=-1; %set B as -1
            outputValue=-2*c;
        end
        evidence=sum(outputValue+randn(1,numTimesteps)); %calculate evidence

        if evidence>0
            decision(trial) = 1; % Decide A
        else
            decision(trial) = -1; % Decide B
        end
    end
    correct2(i)=sum(decision==dispStimulus)/ntrials;
    dprime2(i)=norminv(correct2(i));
end

%Plotting
figure1=figure;
subplot(2,1,1);
plot(contrast, correct2*100,'-o');
xlabel('Stimulus Contrast Level');
ylabel('Percentage of Correctness (%)');
title('Percentage of correctness as a function of stimulus contrast level');

subplot(2,1,2);
plot(contrast, dprime2,'-o');
xlabel('Stimulus Contrast Level');
ylabel('d-prime');
title('d-prime as a function of stimulus contrast level');

% For theoretical prediction, the higher contrast level that a target is presented, the easier
% people can detect it with more details. When there are two targets, the
% higher contrast they are presnted, the easier people can distinguish between them.
%Comparing to the results, percentage of correctness and d-prime increase
%as the contrast of stimuli. Thus, our simulation supported the
%theoretical prediction. 

%% (b) Change the probability of an A stimulus to 0.75. Use a single fixed duration of
% 0.4 s and contrast of 0.1. What is the optimal criterion now? Simulate performance
% with that optimal criterion, and then simulate it with the criterion you used in part (a),
% comparing the performance across the two for 100 trials each. You should be able to
% compute the optimal criterion. However, you can also determine/estimate that
% criterion by doing large-scale simulations for a range of criteria (I’d prefer the closedform
% computation if you can manage it)


stimDur=0.4;
contrast=0.1;
ntrials=1000;

%Initializatioin
correct2=zeros(size(contrast)); %making a blank space
dprime2=zeros(size(contrast)); %making a blank space


%making the experiment for varying duration
for i=1:length(contrast);
    c=contrast(i);
    numTimesteps=round(stimDur/0.1);
    decision=zeros(ntrials,1);  %making a blank space
    dispStimulus=zeros(ntrials,1); %making a blank space to record
    for trial=1:ntrials
        if rand()==0.75
            stimuls='A stimulus';
            dispStimulus(trial)=1; %set A as 1
            outputValue=2*c;
        else
            stimuls='B stimulus';
            dispStimulus(trial)=-1; %set B as -1
            outputValue=-2*c;
        end
        evidence=sum(outputValue+randn(1,numTimesteps)); %calculate evidence

        if evidence>0
            decision(trial) = 1; % Decide A
        else
            decision(trial) = -1; % Decide B
        end
    end
    correct2(i)=sum(decision==dispStimulus)/ntrials;
    dprime2(i)=norminv(correct2(i));
end

%Plotting
figure1=figure;
subplot(2,1,1);
plot(contrast, correct2*100,'-o');
xlabel('Stimulus Contrast Level');
ylabel('Percentage of Correctness (%)');
title('Percentage of correctness as a function of stimulus contrast level');

subplot(2,1,2);
plot(contrast, dprime2,'-o');
xlabel('Stimulus Contrast Level');
ylabel('d-prime');
title('d-prime as a function of stimulus contrast level');

% For theoretical prediction, the higher contrast level that a target is presented, the easier
% people can detect it with more details. When there are two targets, the
% higher contrast they are presnted, the easier people can distinguish between them.
%Comparing to the results, percentage of correctness and d-prime increase
%as the contrast of stimuli. Thus, our simulation supported the
%theoretical prediction. 