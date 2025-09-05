
%Normalization
function [leftEnergyNorm, rightEnergyNorm, upEnergyNorm, downEnergyNorm] = generate_normalization(leftEnergy, rightEnergy, upEnergy, downEnergy)
sdN = 0.02;
sumEnergy = leftEnergy + rightEnergy + upEnergy + downEnergy;
leftEnergyNorm = leftEnergy ./ (sumEnergy + sdN^2);
rightEnergyNorm = rightEnergy ./ (sumEnergy + sdN^2);
upEnergyNorm = upEnergy ./ (sumEnergy + sdN^2);
downEnergyNorm = downEnergy ./ (sumEnergy + sdN^2);
end