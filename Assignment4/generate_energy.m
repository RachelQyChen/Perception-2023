%Energy calculation
function [energyA, energyB] = generate_energy(oddFast, oddSlow, evenFast, evenSlow)
evenLeft = oddFast + evenSlow;
oddLeft = -oddSlow + evenFast;
evenRight = -oddFast + evenSlow;
oddRight = oddSlow + evenFast;
energyA = evenLeft.^2 + oddLeft.^2;
energyB = evenRight.^2. + oddRight.^2;
end