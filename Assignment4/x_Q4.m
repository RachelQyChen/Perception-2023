%% Question 4a
%Right-wards grating
clear; close all; clc;

deltaT = 1; %ms
deltaX=1/120; %spatial sampling rate
duration = 1000; %ms
t = 0:deltaT:duration-deltaT;
xXarray = -2:deltaX:2;
xYarray = -2:deltaX:2;
tau = 25; % ms
sigma = 0.1; %Gaussian sd
sf = 32; %4 cyc/deg for the preferred sinusoid wave
iniPhase = 0; 
n = 10; %10 contrast steps
contrasts = logspace(0.01, 1, n); % contrasts levels
[evenFilt, oddFilt] =generate_gabor(xXarray, sigma, 4);
energy_container = zeros(length(contrasts), 4);
Phase_shift = -2*pi/125; % rightwards shifts


for cc = 1:length (contrasts)
    contrast = contrasts(cc);
    amplitude = logContrastToAmplitude(contrast);

    [evenLeftq, oddLeftq, evenRightq, oddRightq, energyA_q, energyB_q] = Q3d_horizontal_h(xXarray, xYarray, t, deltaT, tau, amplitude, iniPhase, ...
            Phase_shift, sf, oddFilt, evenFilt);
    [evenUpw, oddUpw, evenDownw, oddDownw, energyA_w, energyB_w] = Q3d_vertical_h(xXarray, xYarray, t, deltaT, tau, amplitude, iniPhase, ...
            Phase_shift, sf, oddFilt, evenFilt);

    timeseries = 241; time = 400; 
    [leftEnergyNorm, rightEnergyNorm, upEnergyNorm, downEnergyNorm] = generate_normalization(energyA_q, energyB_q, energyA_w, energyB_w);
    energy_container(cc, :) = [mean(leftEnergyNorm(timeseries, timeseries, time:end)), ...
        mean(rightEnergyNorm(timeseries, timeseries, time:end)), ...
        mean(upEnergyNorm(timeseries, timeseries, time:end)), ...
        mean(downEnergyNorm(timeseries, timeseries, time:end))];
end 

figure();
plot(exp(contrasts), energy_container(:, 1), 'o-', 'DisplayName','Mean of Normalized Left Energy', 'LineWidth', 2);  hold on;
plot(exp(contrasts), energy_container(:, 2), 'o-', 'DisplayName','Mean of Normalized Right Energy', 'LineWidth', 2);  hold on;
plot(exp(contrasts), energy_container(:, 3), 'o-', 'DisplayName','Mean of Normalized Up Energy', 'LineWidth', 2);  hold on;
plot(exp(contrasts), energy_container(:, 4), 'o-', 'DisplayName','Mean of Normalized Down Energy', 'LineWidth', 2);  hold on;
xlabel('Contrast of grating'); ylabel('Normalized Energy'); legend('location', 'best'); 

%% Function

function amplitude = logContrastToAmplitude(logContrast)
linearContrast = exp(logContrast);
amplitude = linearContrast;
end 
%Time filter
function [f1, f2] = Q3_filters(x, t, deltaT, tau)
    
    [x_size, y_size, t_length] = size(x);
    y = zeros(x_size, y_size, t_length, 7);
    f1 = zeros(x_size, y_size, t_length);
    f2 = zeros(x_size, y_size, t_length);

    for tt = 1:t_length - 1
        for type = 1:7
            if type == 1
                deltaY = (deltaT / tau) * (-y(:, :, tt, type) + x(:, :, tt));
                y(:, :, tt + 1, type) = y(:, :, tt, type) + deltaY;
            else
                deltaY = (deltaT / tau) * (-y(:, :, tt, type) + y(:, :, tt, type - 1));
                y(:, :, tt + 1, type) = y(:, :, tt, type) + deltaY;
            end
        end
        % filter 1
        f1(:, :, tt) = y(:, :, tt, 3) - y(:, :, tt, 5); 
         % filter 2
        f2(:, :, tt) = y(:, :, tt, 5) - y(:, :, tt, 7);
    end
end

%Generate Gabor
function [evenFilt, oddFilt] = generate_gabor(x, sigma, sf)
    evenFilt = exp(-(x.^2)./(2*sigma^2)) .* cos(2*pi*sf*x);
    oddFilt = exp(-(x.^2)./(2*sigma^2)) .* sin(2*pi*sf*x);
    integral = sum(evenFilt.^2 + oddFilt.^2);
    evenFilt = evenFilt / integral;
    oddFilt = oddFilt / integral;
end 

%Convolve filter
function [oddFast, evenFast, oddSlow, evenSlow] = temp_gabor(f1, f2, oddFilt, evenFilt)
    [lx, ly, lt] = size(f1);
    oddFast = zeros(lx, ly, lt);
    oddSlow = zeros(lx, ly, lt);
    evenSlow = zeros(lx, ly, lt);
    evenFast = zeros(lx, ly, lt);
    
    for tt = 1:lt
        oddFast(:, :, tt) = conv2(f1(:, :, tt), oddFilt, 'same');
        evenFast(:, :, tt) = conv2(f1(:, :, tt), evenFilt, 'same');
        oddSlow(:, :, tt) = conv2(f2(:, :, tt), oddFilt, 'same');
        evenSlow(:, :, tt) = conv2(f2(:, :, tt), evenFilt, 'same');
    end
end

%Side-selective filter
function [evenLeft, oddLeft, evenRight, oddRight] = selective_filter(oddFast, oddSlow, evenFast, evenSlow)  
    evenLeft = oddFast + evenSlow;
    oddLeft = -oddSlow + evenFast;
    evenRight = -oddFast + evenSlow;
    oddRight = oddSlow + evenFast; 
end

%Energy calculation
function [energyA, energyB] = generate_energy(oddFast, oddSlow, evenFast, evenSlow)
evenLeft = oddFast + evenSlow;
oddLeft = -oddSlow + evenFast;
evenRight = -oddFast + evenSlow;
oddRight = oddSlow + evenFast;
energyA = evenLeft.^2 + oddLeft.^2;
energyB = evenRight.^2. + oddRight.^2;
end



%Horizontal drifting
function [evenLeftq, oddLeftq, evenRightq, oddRightq, energyA_q, energyB_q] = Q3d_horizontal_h(xXarray, xYarray, t, deltaT, tau, amplitude, phase, ...
            Phase_shift, sf, oddFilt, evenFilt)

sinusoid_input = zeros(length(xXarray), length(xYarray), length(t)); 
    
    %generating a sinusoid input
    for tt = 1:1000
        phase = phase + Phase_shift;
        [Horizontal_vale, Vertical_value] = meshgrid(xXarray .* sf);
        sinusoid_input(:, :, tt) = amplitude * sin(Horizontal_vale + phase);
    end 
    
    %Applying temporal filter first
    [f1, f2] = Q3_filters(sinusoid_input, t, deltaT, tau);
    
    %Convolving filters by honrizontal filters
    [oddFastq, evenFastq, oddSlowq, evenSlowq] = temp_gabor(f1, f2, oddFilt, evenFilt);

    %Side-selective
    [evenLeftq, oddLeftq, evenRightq, oddRightq] = selective_filter(oddFastq, oddSlowq, evenFastq, evenSlowq);

    %Compute energie for horizontal
    [energyA_q, energyB_q] = generate_energy(oddFastq, oddSlowq, evenFastq, evenSlowq);

end

%Vertical drifting
function [evenUpw, oddUpw, evenDownw, oddDownw, energyA_w, energyB_w] = Q3d_vertical_h(xXarray, xYarray, t, deltaT, tau, amplitude, phase, ...
            Phase_shift, sf, oddFilt, evenFilt)

sinusoid_input2 = zeros(length(xXarray), length(xYarray), length(t)); 
    
    %generating a sinusoid input
    for tt = 1:1000
        phase = phase + Phase_shift;
        [Horizontal_vale, Vertical_value] = meshgrid(xXarray .* sf);
        sinusoid_input2(:, :, tt) = amplitude * sin(Horizontal_vale + phase);
    end 
    
    %Applying temporal filter first
    [f1, f2] = Q3_filters(sinusoid_input2, t, deltaT, tau);
    
    %Convolving filters by honrizontal filters
    [oddFastw, evenFastw, oddSloww, evenSloww] = temp_gabor(f1, f2, oddFilt, evenFilt);

    %Side-selective
    [evenUpw, oddUpw, evenDownw, oddDownw] = selective_filter(oddFastw, evenFastw, oddSloww, evenSloww);

    %Compute energie for horizontal
    [energyA_w, energyB_w] = generate_energy(oddFastw, oddSloww, evenFastw, evenSloww);

end

%Normalization
function [leftEnergyNorm, rightEnergyNorm, upEnergyNorm, downEnergyNorm] = generate_normalization(leftEnergy, rightEnergy, upEnergy, downEnergy)
sdN = 0.02;
sumEnergy = leftEnergy + rightEnergy + upEnergy + downEnergy;
leftEnergyNorm = leftEnergy ./ (sumEnergy + sdN^2);
rightEnergyNorm = rightEnergy ./ (sumEnergy + sdN^2);
upEnergyNorm = upEnergy ./ (sumEnergy + sdN^2);
downEnergyNorm = downEnergy ./ (sumEnergy + sdN^2);
end







