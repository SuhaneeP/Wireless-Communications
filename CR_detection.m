close all;
clc;

L = 100; % u (Channel Length)
false_alarm_probability = 0.01; % Assumed probability of fasle alarm used for determining threshold value
threshold = (qfuncinv(false_alarm_probability)./sqrt(L))+1; % calculating threshold from Pf
snr_dB = -20:5; % range of SNR
snr_linear = power(10,snr_dB/10); % Converting value of SNR from dBs to linear

for i=1:length(snr_dB)
    detection = 0; % initialising with 0
    for k=1:10000 % Number of Monte-Carlo Simulations
        noise = randn(1,L); % noise generation (Real AWGN Noise)
        signal = sqrt(snr_linear(i)).*randn(1,L); % Input signal*gain
        rsignal = signal+noise; % received signal
        energy = power(abs(rsignal),2); % Energy Calculation
        test_statistic = sum(energy)/L; % Calculating the Test Statistic
        if(test_statistic >= threshold) % using the definition of probability of detection
            detection = detection + 1;
        end
    end
    Pd_simulation(i) = detection / k; % calculating detection probability
end

numerator = L.*(threshold - (snr_linear+1));
denominator = sqrt(2*L*(snr_linear+1));
Pd_analytical = qfunc(numerator./denominator); % analytical formula for finding Pd

close all;
plot(snr_dB,Pd_analytical,'r.-','LineWidth',1.5);
hold on
plot(snr_dB,Pd_simulation,'b*');
axis([-20 5 0 1])

LOC = "northwest";
legend('Analytical','Simulation','Location',LOC);
title('$P_{d} \;vs\; SNR \;for\; P_{f} = 0.01$','Interpreter','latex','FontWeight','Normal');
xlabel('$\frac{E_b}{N_0} (dBs)$','Interpreter','latex');
ylabel('Probability of Detection');
grid on