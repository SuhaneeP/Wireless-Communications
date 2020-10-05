% Sample size (L) = 1000;
% Pf=0.01

close all;
clear all;
clc;

SNR_dB = (-50:10); %SNR in db
SNR_Lin = 10.^(0.1*(SNR_dB)); % SNR Linear values

L = 1000; % Time bandwidth product
Pf = 0.01; % Prob of false alarm

% Analytical expression
threshold = (qfuncinv(Pf)./sqrt(L)) + 1; %lambda
theory_Pd = qfunc(((threshold-(SNR_Lin+1)).*L)./sqrt(2.*L.*(SNR_Lin+1).^2)); % Probability of detection



for i = 1:length(SNR_dB)
    detect = 0; % variable to check detection of error
    % Monte Carlo Simulation
    for ii = 1:10000
        noise = randn(1,L); % white gaussian noise
        signal = sqrt(SNR_Lin(i)).*randn(1,L);
        rec_sig = signal + noise; % y = hx + n; Here, h=1
        energy = abs((rec_sig).^2);% compute the received signal energy over L samples
        test = (1/L).*sum(energy); % average energy
        threshold(i) = (qfuncinv(Pf)./sqrt(L)) + 1;
        if(test >= threshold(i)) % check whether the received signal is greater than the threshold
            detect = detect + 1; % detecting of condition
        end
    end
    sim_Pd(i) = detect/ii; % compute prob of detection
end

% Plots
close all
figure
plot(SNR_dB, theory_Pd, 'r-', 'LineWidth', 1);
hold on;
plot(SNR_dB, sim_Pd, 'b+', 'LineWidth', 1);
axis([-20 5 0 1])
grid on
legend('Theory', 'Simulation');
xlabel('Signal to Noise Ratio (dB)');
ylabel('Probability of Detection');
title('Pd v/s SNR for Pf=0.01');
