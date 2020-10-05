close all;
clear all;
clc;

N = 10^6;  % Simulation

snr_db = (0:20);
s1_inp = rand(1,N) > 0.5;   % Generating random 0 and 1
s2_inp = rand(1,N) > 0.5;
s1 = 2*s1_inp - 1;
s2 = 2*s2_inp - 1;
s_inp = s1 + 1i*s2; % QPSK modulated signal

for i = 1:length(snr_db)
    n = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); % AWGN noise
    h = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); % Rayleigh channel
    s_noise = s_inp + (power(10,-snr_db(i)/20) * n);   % Adding noise
    s_noise_wireless = h.*s_inp + (power(10,-snr_db(i)/20) * n);    % Adding noise for wireless tarnsmission
    
    % For AWGN 
    s1_out = real(s_noise) > 0; % Decoding
    s2_out = imag(s_noise) > 0;
    
    diff1 = xor(s1_inp,s1_out); % Counting error
    diff2 = xor(s2_inp,s2_out);
    err1 = sum(diff1(:)==1);
    err2 = sum(diff2(:)==1);
    ber(i) = mean([err1 err2]); % Taking mean
    
    % For wireless
    s_noise_wireless_e = s_noise_wireless./h; % Equalization 
    s1_out_wireless = real(s_noise_wireless_e) > 0; % Decoding
    s2_out_wireless = imag(s_noise_wireless_e) > 0;
    
    diff1_wireless = xor(s1_inp,s1_out_wireless);    % Counting error
    diff2_wireless = xor(s2_inp,s2_out_wireless);
    err1_wireless = sum(diff1_wireless(:)==1);
    err2_wireless = sum(diff2_wireless(:)==1);
    ber_wireless(i) = mean([err1_wireless err2_wireless]);  % Taking mean
end

ber_qpsk = ber/N;   % simulated BER for QPSK modulation
ser_qpsk = 2*ber_qpsk;  % simulated SER for QPSK modulation

ber_qpsk_wireless = ber_wireless/N; % simulated BER for QPSK modulation for wireless
ser_qpsk_wireless = 2*ber_qpsk_wireless;    % simulated SER for QPSK modulation for wireless

close all;
figure;
semilogy(snr_db, ser_qpsk, 'ro-');  % Plot for AWGN channel
hold on
semilogy(snr_db, ser_qpsk_wireless, 'm^-'); % Plot for Rayleigh Channel
hold off

axis([0 20 10^-5 1]);
grid on;

legend('QPSK over AWGN Channel','QPSK over Rayleigh Channel')
title('Simulated plot over Wireless Channel and comparision with AWGN channel.');
xlabel('Eb/No_{db}');
ylabel('Symbol Error Rate');