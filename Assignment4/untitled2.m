
for ii = 1:length(contrasts)
    contrast = contrasts(ii);
    [leftEven, leftOdd, rightEven, rightOdd, leftEnergy, rightEnergy, ...
                upEven, upOdd, downEven, downOdd, upEnergy, downEnergy] = ...
                neuron_responses(x_x, x_y, t, deltaT, tau, contrast, phase, ...
                phase_shift, sf, ori, oddFilt, evenFilt);
    [leftEnergyNorm, rightEnergyNorm, upEnergyNorm, downEnergyNorm] = ...
        energy_norm(leftEnergy, rightEnergy, upEnergy, downEnergy);
    
    x_y_dim = 241;
    st_tm = 400; % ms
    mean_energies(ii, :) = [mean(leftEnergyNorm(x_y_dim, x_y_dim, st_tm:end)), ...
        mean(rightEnergyNorm(x_y_dim, x_y_dim, st_tm:end)), ...
        mean(upEnergyNorm(x_y_dim, x_y_dim, st_tm:end)), ...
        mean(downEnergyNorm(x_y_dim, x_y_dim, st_tm:end))];
end

plt_titles = ["leftEnergyNorm mean", "rightEnergyNorm mean", ...
    "upEnergyNorm mean", "downEnergyNorm mean"];

figure()
for ss = 1:size(mean_energies, 2)
    plot(contrasts, mean_energies(:, ss), 'o-', ...
        'DisplayName', plt_titles(ss), 'LineWidth', 2);
    xlabel('Contrast of grating')
    ylabel('Normalized Energy')
    legend('location', 'northwest')
    hold on;
end



function [leftEnergyNorm, rightEnergyNorm, upEnergyNorm, downEnergyNorm] = ...
    energy_norm(leftEnergy, rightEnergy, upEnergy, downEnergy)
    sumEnergy = leftEnergy + rightEnergy + upEnergy + downEnergy;
    stdev = 0.02;
    leftEnergyNorm = leftEnergy ./ (sumEnergy + stdev^2);
    rightEnergyNorm = rightEnergy ./ (sumEnergy + stdev^2);
    upEnergyNorm = upEnergy ./ (sumEnergy + stdev^2);
    downEnergyNorm = downEnergy ./ (sumEnergy + stdev^2);
end