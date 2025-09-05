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