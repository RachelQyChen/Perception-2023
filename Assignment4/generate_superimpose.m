%Superimposing
function [leftEven, leftOdd, rightEven, rightOdd, leftEnergy, rightEnergy, upEven, upOdd, downEven, downOdd, upEnergy, downEnergy] = ...
            generate_superimpose(xXarray, xYarray, t, deltaT, tau, amplitude, phase, Phase_shift, sf, oddFilt, evenFilt)

superimpose_input= zeros(length(xXarray), length(xYarray), length(t)); 

for tt = 1:length(t)
    [Horizontal_vale, Vertical_value] = meshgrid(xXarray .* sf);
        phase_right = phase - Phase_shift;
        rightwardsinput = amplitude * sin(Horizontal_vale + phase_right);
        phase_up = phase + Phase_shift;
        upwardsinput = 0.5 * sin(Vertical_value + phase_up);
        superimpose_input(:, :, tt) = rightwardsinput + upwardsinput;
end 
[f1, f2] = Q3_filters(superimpose_input, t, deltaT, tau);

 %Convolving filters by honrizontal filters
    [oddFastq, evenFastq, oddSlowq, evenSlowq] = temp_gabor(f1, f2, oddFilt, evenFilt);
[oddFastw, evenFastw, oddSloww, evenSloww] = temp_gabor(f1, f2, oddFilt, evenFilt);
    %Side-selective
    [leftEven, leftOdd, rightEven, rightOdd] = selective_filter(oddFastq, oddSlowq, evenFastq, evenSlowq);
[upEven, upOdd, downEven, downOdd] = selective_filter(oddFastw, evenFastw, oddSloww, evenSloww);
    
%Compute energie for horizontal
    [leftEnergy, rightEnergy] = generate_energy(oddFastq, oddSlowq, evenFastq, evenSlowq);
    [upEnergy, downEnergy] = generate_energy(oddFastw, oddSloww, evenFastw, evenSloww);
end

