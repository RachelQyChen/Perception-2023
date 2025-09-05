%% Question 4b
%Right-wards grating
clear all; close all; clc;

deltaT = 1; %ms
deltaX=1/120; %spatial sampling rate
duration = 1000; %ms
t = 0:deltaT:duration-deltaT;
xXarray = -2:deltaX:2;
xYarray = -2:deltaX:2;
tau = 25; % ms
sigma = 0.1; %Gaussian sd
sf = 30; 
iniPhase = 0; 
contrasts_right = [1, 5, 10, 15, 25, 45]; % contrasts levels
%contrasts_right = [1, 10,  45];
contrast_up = 50; %50%
[evenFilt, oddFilt] =generate_gabor(xXarray, sigma, 4);
Phase_shift_R = -2*pi/125; % rightwards shifts
Phase_shift = 2*pi/125; % rightwards shifts

%Whether there is upwards signals or whether there is not upwards signals
energy_container_no = zeros(length(contrasts_right), 4);
energy_container_yes = zeros(length(contrasts_right), 4);

%When there is no upwards
for cc = 1:length (contrasts_right)
    contrast = contrasts_right(cc);
    amplitude_right = contrast*(1/100);

    [evenLeftq, oddLeftq, evenRightq, oddRightq, energyA_q, energyB_q] = Q3d_horizontal_h(xXarray, xYarray, t, deltaT, tau, amplitude_right, iniPhase, ...
        Phase_shift_R, sf, oddFilt, evenFilt);
    [evenUpw, oddUpw, evenDownw, oddDownw, energyA_w, energyB_w] = Q3d_vertical_h(xXarray, xYarray, t, deltaT, tau, amplitude_right, iniPhase, ...
        Phase_shift_R, sf, oddFilt, evenFilt);

    timeseries = 241; time = 500;
    [leftEnergyNorm, rightEnergyNorm, upEnergyNorm, downEnergyNorm] = generate_normalization(energyA_q, energyB_q, energyA_w, energyB_w);
    energy_container_no(cc, :) = [mean(leftEnergyNorm(timeseries, timeseries, time:end)), ...
        mean(rightEnergyNorm(timeseries, timeseries, time:end)), ...
        mean(upEnergyNorm(timeseries, timeseries, time:end)), ...
        mean(downEnergyNorm(timeseries, timeseries, time:end))];
end

%When there is upwards
for kk = 1:length (contrasts_right)
    contrast = contrasts_right(kk);
    amplitude_right = contrast*(1/100);

    [leftEven, leftOdd, rightEven, rightOdd, leftEnergy, rightEnergy, upEven, upOdd, downEven, downOdd, upEnergy, downEnergy] = ...
        generate_superimpose(xXarray, xYarray, t, deltaT, tau, amplitude_right, iniPhase, Phase_shift, sf, oddFilt, evenFilt);

    timeseries = 241; time = 500;
    [leftEnergyNorm, rightEnergyNorm, upEnergyNorm, downEnergyNorm] = generate_normalization(leftEnergy, rightEnergy, upEnergy, downEnergy);
    energy_container_yes(kk, :) = [mean(leftEnergyNorm(timeseries, timeseries, time:end)), ...
        mean(rightEnergyNorm(timeseries, timeseries, time:end)), ...
        mean(upEnergyNorm(timeseries, timeseries, time:end)), ...
        mean(downEnergyNorm(timeseries, timeseries, time:end))];
end

figure();
subplot(2,2,1)
plot (contrasts_right, energy_container_no(:,1), 'o-', 'DisplayName', 'Single', 'LineWidth', 2); hold on;
plot (contrasts_right, energy_container_yes(:,1), 'o-', 'DisplayName', 'Superimpose', 'LineWidth', 2); hold on;
legend ('Location','best'); title('Mean of Normalized Up Energy'); hold on;

subplot(2,2,2)
plot (contrasts_right, energy_container_no(:,2), 'o-', 'DisplayName', 'Single', 'LineWidth', 2); hold on;
plot (contrasts_right, energy_container_yes(:,2), 'o-', 'DisplayName', 'Superimpose', 'LineWidth', 2); hold on;
legend ('Location','best'); title('Mean of Normalized Down Energy'); hold on;

subplot(2,2,3)
plot (contrasts_right, energy_container_no(:,3), 'o-', 'DisplayName', 'Single', 'LineWidth', 2); hold on;
plot (contrasts_right, energy_container_yes(:,3), 'o-', 'DisplayName', 'Superimpose', 'LineWidth', 2); hold on;
legend ('Location','best'); title('Mean of Normalized Left Energy'); hold on;

subplot(2,2,4)
plot (contrasts_right, energy_container_no(:,4), 'o-', 'DisplayName', 'Single', 'LineWidth', 2); hold on;
plot (contrasts_right, energy_container_yes(:,4), 'o-', 'DisplayName', 'Superimpose', 'LineWidth', 2); hold on;
legend ('Location','best'); title('Mean of Normalized Right Energy'); hold on;

xlabel('Contrast Levels'); ylabel('Normalized Energy'); legend('location', 'best'); 


