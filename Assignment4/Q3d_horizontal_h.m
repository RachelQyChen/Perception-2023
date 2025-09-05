
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
