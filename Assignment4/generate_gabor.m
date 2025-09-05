%Generate Gabor
function [evenFilt, oddFilt] = generate_gabor(x, sigma, sf)
    evenFilt = exp(-(x.^2)./(2*sigma^2)) .* cos(2*pi*sf*x);
    oddFilt = exp(-(x.^2)./(2*sigma^2)) .* sin(2*pi*sf*x);
    integral = sum(evenFilt.^2 + oddFilt.^2);
    evenFilt = evenFilt / integral;
    oddFilt = oddFilt / integral;
end 