%% Question 3a
deltaT = 1; %ms
deltaX=1/120; %spatial sampling rate
duration = 1000; %ms
t = 0:deltaT:duration-deltaT;
xXarray = -2:deltaX:2;
xYarray = -2:deltaX:2;
tau = 25; % ms
sigma = 0.1; %Gaussian sd
sf = 4; %4 cyc/deg for the preferred sinusoid wave

%signal
x = zeros(length(xXarray), length(xYarray), length(t)); 
x(241, 241, 1) = 1; %input of 480 units, and 240 each onleft and right side of zero, +1so started from oringin

%Filtering
[f1, f2] = Q3_filters(x, t, deltaT, tau); % computes temporal filters

%Generate gabor
[evenFilt, oddFilt] = generate_gabor(xXarray, sigma, sf); 

%Temoral and Horizontal spatial filters
[oddFast_h, evenFast_h, oddSlow_h, evenSlow_h] = temp_gabor(f1, f2, oddFilt, evenFilt); 
upLeft = squeeze(oddFast_h(241, :, :))';
upRight = squeeze(oddSlow_h(241, :, :))';
downLeft = squeeze(evenFast_h(241, :, :))';
downRight = squeeze(evenSlow_h(241, :, :))';

 figure1 = figure; 
 subplot(2,2,1); imagesc(upLeft); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('oddFast (Horizontal)'); 

 subplot(2,2,2); imagesc(upRight); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('oddSlow (Horizontal)'); 

 subplot(2,2,3); imagesc(downLeft); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('evenFast (Horizontal)'); 

 subplot(2,2,4); imagesc(downRight); colormap(gray);
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('evenSlow (Horizontal)'); 

%Temoral and Vertical spatial filters
[oddFast_v, evenFast_v, oddSlow_v, evenSlow_v] = temp_gabor(f1, f2, oddFilt', evenFilt'); 
upLeftv = squeeze(oddFast_v(:, 241, :))';
upRightv = squeeze(oddSlow_v(:, 241, :))';
downLeftv = squeeze(evenFast_v(:, 241, :))';
downRightv = squeeze(evenSlow_v(:, 241, :))';

 figure2 = figure; 
 subplot(2,2,1); imagesc(upLeftv); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('oddFast (Vertical)'); 

 subplot(2,2,2); imagesc(upRightv); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('oddSlow (Vertical)'); 

 subplot(2,2,3); imagesc(downLeftv); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('evenFast (Vertical)'); 

 subplot(2,2,4); imagesc(downRightv); colormap(gray);
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('evenSlow (Vertical)'); 

 %% Question 3b

 %Horizontal 
 [leftEven, leftOdd, rightEven, rightOdd] = selective_filter(oddFast_h, oddSlow_h, evenFast_h, evenSlow_h); 
 upLeft1 = squeeze(leftEven(241, :, :))';
 upRight1 = squeeze(leftOdd(241, :, :))';
 downLeft1 = squeeze(rightEven(241, :, :))';
 downRight1 = squeeze(rightOdd(241, :, :))';

 figure3 = figure; 
 subplot(2,2,1); imagesc(upLeft1); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Left Even (Horizontal)'); 

 subplot(2,2,2); imagesc(upRight1); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Left Odd (Horizontal)'); 

 subplot(2,2,3); imagesc(downLeft1); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Right Even (Horizontal)'); 

 subplot(2,2,4); imagesc(downRight1); colormap(gray);
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Right Odd (Horizontal)');


 %Vertical
[upEven, upOdd, downEven, downOdd] = selective_filter(oddFast_v, oddSlow_v, evenFast_v, evenSlow_v); 
upLeftv1 = squeeze(upEven(:, 241, :))';
upRightv1 = squeeze(upOdd(:, 241, :))';
downLeftv1 = squeeze(downEven(:, 241, :))';
downRightv1 = squeeze(downOdd(:, 241, :))';

 figure4 = figure; 
 subplot(2,2,1); imagesc(upLeftv1); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Left Even (Vertical)'); 

 subplot(2,2,2); imagesc(upRightv1); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Left Odd (Vertical)'); 

 subplot(2,2,3); imagesc(downLeftv1); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Right Even (Vertical)'); 

 subplot(2,2,4); imagesc(downRightv1); colormap(gray);
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Right Odd (Vertical)');

  %% Question 3c
[energyA_h, energyB_h] = generate_energy(oddFast_h, oddSlow_h, evenFast_h, evenSlow_h);
upLeft_energy = squeeze(energyA_h(241, :, :))';
upRighty_energy = squeeze(energyB_h(241, :, :))';

[energyA_v, energyB_v] = generate_energy(oddFast_v, oddSlow_v, evenFast_v, evenSlow_v);
downLeft_energy = squeeze(energyA_v(:, 241, :))';
downRighty_energy = squeeze(energyB_v(:, 241, :))';

figure5 = figure; 
 subplot(2,2,1); imagesc(upLeft_energy); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Left Even (Vertical)'); 

 subplot(2,2,2); imagesc(upRighty_energy); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Left Odd (Vertical)'); 

 subplot(2,2,3); imagesc(downLeft_energy); colormap(gray); 
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Right Even (Vertical)'); 

 subplot(2,2,4); imagesc(downRighty_energy); colormap(gray);
 xticks([1, 121, 241, 361, 481]); xticklabels([-2, -1, 0, 1, 2]); yticks(0:100:1000)
 xlabel('Visual angle (deg)'); ylabel('time (ms)'); title('Right Odd (Vertical)');

 %% Question 3d
iniPhase = 0;
shiftPhases = [2*pi/125, -2*pi/125];
sf = 32;
amplitude = 1; 

%Horizontally
for Phase_shift = shiftPhases

[evenLeftq, oddLeftq, evenRightq, oddRightq, energyA_q, energyB_q] = Q3d_horizontal_h(xXarray, xYarray, t, deltaT, tau, amplitude, iniPhase, ...
            Phase_shift, sf, oddFilt, evenFilt);
[evenUpw, oddUpw, evenDownw, oddDownw, energyA_w, energyB_w] = Q3d_vertical_h(xXarray, xYarray, t, deltaT, tau, amplitude, iniPhase, ...
            Phase_shift, sf, oddFilt, evenFilt);

if  Phase_shift > 0
    direction_h = 'Left';
else
    direction_h = 'Right';
end

%Plotting
figure6 = figure;
sgtitle(direction_h)

subplot (2, 2, 1);
plot(t, squeeze(evenLeftq(241, 241, :)), 'DisplayName', 'Even Left'); hold on;
plot(t, squeeze(oddLeftq(241, 241, :)), 'DisplayName', 'Odd Left'); hold on;
plot(t, squeeze(energyA_q(241, 241, :)), 'DisplayName', 'Energy Left'); hold on;
xlim([0 500]); xlabel('Time (ms)'); ylabel('Response'); legend()

subplot (2, 2, 2);
plot(t, squeeze(evenRightq(241, 241, :)), 'DisplayName', 'Even Right'); hold on;
plot(t, squeeze(oddRightq(241, 241, :)), 'DisplayName', 'Odd Right'); hold on;
plot(t, squeeze(energyB_q(241, 241, :)), 'DisplayName', 'Energy Right'); hold on;
xlim([0 500]); xlabel('Time (ms)'); ylabel('Response'); legend()

subplot (2, 2, 3);
plot(t, squeeze(evenUpw(241, 241, :)), 'DisplayName', 'Even Up'); hold on;
plot(t, squeeze(oddUpw(241, 241, :)), 'DisplayName', 'Odd Up'); hold on;
plot(t, squeeze(energyA_w(241, 241, :)), 'DisplayName', 'Energy Up'); hold on;
xlim([0 500]); xlabel('Time (ms)'); ylabel('Response'); legend()

subplot (2, 2, 4);
plot(t, squeeze(evenDownw(241, 241, :)), 'DisplayName', 'Even Down'); hold on;
plot(t, squeeze(oddDownw(241, 241, :)), 'DisplayName', 'Odd Down'); hold on;
plot(t, squeeze(energyB_w(241, 241, :)), 'DisplayName', 'Energy Down'); hold on;
xlim([0 500]); xlabel('Time (ms)'); ylabel('Response'); legend()
end

%% Vertically
for Phase_shift = shiftPhases

[evenLeftq, oddLeftq, evenRightq, oddRightq, energyA_q, energyB_q] = Q3d_horizontal_v(xXarray, xYarray, t, deltaT, tau, amplitude, iniPhase, ...
            Phase_shift, sf, oddFilt, evenFilt);
[evenUpw, oddUpw, evenDownw, oddDownw, energyA_w, energyB_w] = Q3d_vertical_v(xXarray, xYarray, t, deltaT, tau, amplitude, iniPhase, ...
            Phase_shift, sf, oddFilt, evenFilt);

if  Phase_shift > 0
    direction_v = 'Up';
else
    direction_v = 'Down';
end

%Plotting
figure7 = figure;
sgtitle(direction_v)

subplot (2, 2, 1);
plot(t, squeeze(evenLeftq(241, 241, :)), 'DisplayName', 'Even Left'); hold on;
plot(t, squeeze(oddLeftq(241, 241, :)), 'DisplayName', 'Odd Left'); hold on;
plot(t, squeeze(energyA_q(241, 241, :)), 'DisplayName', 'Energy Left'); hold on;
xlim([0 500]); xlabel('Time (ms)'); ylabel('Response'); legend()

subplot (2, 2, 2);
plot(t, squeeze(evenRightq(241, 241, :)), 'DisplayName', 'Even Right'); hold on;
plot(t, squeeze(oddRightq(241, 241, :)), 'DisplayName', 'Odd Right'); hold on;
plot(t, squeeze(energyB_q(241, 241, :)), 'DisplayName', 'Energy Right'); hold on;
xlim([0 500]); xlabel('Time (ms)'); ylabel('Response'); legend()

subplot (2, 2, 3);
plot(t, squeeze(evenUpw(241, 241, :)), 'DisplayName', 'Even Up'); hold on;
plot(t, squeeze(oddUpw(241, 241, :)), 'DisplayName', 'Odd Up'); hold on;
plot(t, squeeze(energyA_w(241, 241, :)), 'DisplayName', 'Energy Up'); hold on;
xlim([0 500]); xlabel('Time (ms)'); ylabel('Response'); legend()

subplot (2, 2, 4);
plot(t, squeeze(evenDownw(241, 241, :)), 'DisplayName', 'Even Down'); hold on;
plot(t, squeeze(oddDownw(241, 241, :)), 'DisplayName', 'Odd Down'); hold on;
plot(t, squeeze(energyB_w(241, 241, :)), 'DisplayName', 'Energy Down'); hold on;
xlim([0 500]); xlabel('Time (ms)'); ylabel('Response'); legend()
end


%% Function

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

%Horizontal drifting
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

%Vertical drifting
function [evenLeftq, oddLeftq, evenRightq, oddRightq, energyA_q, energyB_q] = Q3d_horizontal_v(xXarray, xYarray, t, deltaT, tau, amplitude, phase, ...
            Phase_shift, sf, oddFilt, evenFilt)

sinusoid_input = zeros(length(xXarray), length(xYarray), length(t)); 
    
    %generating a sinusoid input
    for tt = 1:1000
        phase = phase + Phase_shift;
        [Horizontal_vale, Vertical_value] = meshgrid(xXarray .* sf);
        sinusoid_input(:, :, tt) = amplitude * sin(Vertical_value + phase);
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
function [evenUpw, oddUpw, evenDownw, oddDownw, energyA_w, energyB_w] = Q3d_vertical_v(xXarray, xYarray, t, deltaT, tau, amplitude, phase, ...
            Phase_shift, sf, oddFilt, evenFilt)

sinusoid_input2 = zeros(length(xXarray), length(xYarray), length(t)); 
    
    %generating a sinusoid input
    for tt = 1:1000
        phase = phase + Phase_shift;
        [Horizontal_vale, Vertical_value] = meshgrid(xXarray .* sf);
        sinusoid_input2(:, :, tt) = amplitude * sin(Vertical_value + phase);
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
