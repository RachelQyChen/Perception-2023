%% Perception HW4
%Name: Rachel Chen
clear all; clc; close all;

%% Question 1.a
%Ground setting up
deltaT = 1; %ms
duration = 1000; % ms
t = 0:deltaT:duration-deltaT;
x = zeros(size(t));
x(100) = 1;
tau = 25; % ms

y1 = iirFilter(x, t, deltaT, tau);
figure;
plot(t, y1, 'b-');
xlabel('Time (ms)');
ylabel('Output');
title('Impulse response through IIR Filter')

% The impulse at time-step = 100 ms
x = zeros(size(t));
x(100) = 1;
tauVals = [50, 75, 10, 5]; % different value of tau

figure;
for a = 1:length(tauVals)
    tau = tauVals(a); 
    y1 = iirFilter(x, t, deltaT, tau);

    subplot(length(tauVals),1,  a);
    plot(t, y1, 'bo-'); hold on;
    
    % Lotting the exponential function
    exponential_function = exp(-(t./tau)); 
    plot(exponential_function, 'r-');

    xlabel('Time (ms)')
    ylabel('Output')
    ylim([0, max(2*y1)])
    title(['Impulse response through IIR Filter when tau = ', num2str(tau)])
    legend('y-save', 'expoential y');
end

%nomailize ysave (y1), y1./max(y1)

%%The overall tendency of impulse response output have a same tendency as
%%the exponential functioin (can be fitted by the exponential function).
%%However, the output can be recalibrated to 

%% Question 1.b
x = zeros(size(t));
x(100:1000) = 1;
tauVals = [50, 75, 10]; % different value of tau

figure;
for a = 1:length(tauVals)
    tau = tauVals(a); 
    y1 = iirFilter(x, t, deltaT, tau);

    subplot(1,length(tauVals),  a);
    plot(t, y1, 'bo-'); hold on;
    
    % Lotting the exponential function
    exponential_function = 1 - exp(-(t./tau)); 
    plot(exponential_function, 'r-');

    xlabel('Time (ms)')
    ylabel('Output')
    ylim([0, 1.3])
    legend('y-save in step response', 'exponential y')
    
end

%% Question 1.c
%According to handout, the differential equation is: τ dy(t)/dt = -y(t) + x(t);
%Taking the Fourier transform of both sides, the equation in fourier domain
%wil be: d-hat(w)y-hat(w) = -y-hat(w)+x-hat(w).
%The effective frequency response is: h-hat(w) = 1/[1 + d-hat(w)];
%Thus, the equation now will be: y-hat(w) = x-hat(w)/[1+d-hat(w)];
frequencys = [2, 10, 20];
amplitude_theory = 1;
phase_theory = 1.5;

for uu = 1:length(frequencys)

frequency = frequencys(uu);
%generate sinusodio wave
x_sinusoid =sin(2*pi* frequency/1000*t+phase_theory);
%simulated output
y1_sinusoid = iirFilter(x_sinusoid, t, deltaT, tau); 

amplitude_response = amplitude_theory / (1+2*pi*frequency/1000); 
phase_response = phase_theory / (pi / 2 + phase_theory); 
theory_prediction = amplitude_response*sin(2*pi*frequency/1000*t+phase_response); 

subplot(3,1,uu);
plot (t, theory_prediction, 'DisplayName', 'Theoretical prediction'); hold on
plot(t, y1_sinusoid , 'r','DisplayName', 'Simulated results'); hold on;
legend(); 
end


%% Question 2
deltaT = 1; %ms
duration = 1000; % ms
t = [0:deltaT:duration-deltaT];
x = zeros(size(t));
x(100) = 1;
tau = 25; % ms

[f1, f2] = lpFilter(x, t, deltaT, tau); 
figure;
plot(t, f1, 'DisplayName', 'Filter 1'); hold on;
plot(t, f2, 'DisplayName', 'Filter 2')
xlabel('time (ms)')
title('Temporal filters')
legend()

%% Functions
%Impulse response filter
%This function compute the response obtained through applying IIR filter
function y1 = iirFilter(x, t, deltaT, tau)
y1 = zeros(length(t),1);
for tt = 1:length(t) - 1
    deltaY1 = (deltaT/tau) * (-y1(tt) + x(tt));
    y1(tt + 1) = y1(tt) + deltaY1;
end
end

function [f1, f2] = lpFilter(x, t, deltaT, tau)
y = zeros(length(t),7);
for tt = 1:length(t) - 1
    for type = 1:7
        if type == 1 || type == 2 
            deltaY = (deltaT/tau) * (-y(tt, type) + x(tt));
            y(tt + 1, type) = y(tt, type) + deltaY;
        else
            deltaY = (deltaT / tau) * (-y(tt, type) + y(tt, type - 1));
            y(tt + 1, type) = y(tt, type) + deltaY;
        end
    end
end
%filter 1
f1 = y(:, 3) - y(:, 5);
%filter 2
f2 = y(:, 5) - y(:, 7); % Slow filter
end

