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