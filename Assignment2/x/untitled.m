% Parameters
contrast = 0.2;
durations = 100:100:2000;  % from 100 ms to 2 s

% Pre-allocate arrays for results
pcorrect = zeros(size(durations));
dprime = zeros(size(durations));

% Loop over stimulus durations
for i = 1:length(durations)
    duration = durations(i);
    nSteps = duration / 100;  % number of time steps
    nTrials = 1000;  % number of trials (arbitrary, adjust as needed)
    
    % Simulate evidence collection
    evidence = zeros(nTrials, 2);  % columns for stimulus A and B
    for t = 1:nTrials
        for step = 1:nSteps
            noise = randn;  % random draw from standard Gaussian
            evidence(t, 1) = evidence(t, 1) + (2 * contrast + noise);  % evidence for stimulus A
            evidence(t, 2) = evidence(t, 2) + (-2 * contrast + noise);  % evidence for stimulus B
        end
    end
    
    % Determine responses based on collected evidence
    responses = evidence(:, 1) > evidence(:, 2);  % 1 if stimulus A chosen, 0 if stimulus B chosen
    
    % Stimulus labels (0 for B and 1 for A, equally likely)
    stimulusLabels = rand(nTrials, 1) > 0.5;
    
    % Calculate percentage correct
    pcorrect(i) = sum(responses == stimulusLabels) / nTrials;
    
    % Calculate hit rate and false alarm rate
    hits = sum(responses == 1 & stimulusLabels == 1);
    falseAlarms = sum(responses == 1 & stimulusLabels == 0);
    hitRate = hits / sum(stimulusLabels == 1);
    falseAlarmRate = falseAlarms / sum(stimulusLabels == 0);
    
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
