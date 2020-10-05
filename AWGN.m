clear all;
close all;
clc;

N = 10^6;

%BPSK

snr_db = (1:10);
s_inp = rand(1,N) > 0.5;   % Generate random 1 and 0
sb = 2*s_inp - 1;   % BPSK modulated 

for i = 1:length(snr_db)
    n = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));  % AWGN Noise
    sb_n = sb + 10^(-snr_db(i)/20)*n;    % Noise addition AWGN  
    s_out = real(sb_n) > 0;    % Decoding
    diff = xor(s_inp,s_out); % Checking error
    err(i) = sum(diff(:)==1); % Counting errors
end

bpsk_ber = err/N;    % simulated BER for BPSK modulation

%QPSK

s1_inp = rand(1,N) > 0.5;   % Generate random 1 and 0
s2_inp = rand(1,N) > 0.5;   % Generate random 1 and 0
s1 = 2*s1_inp - 1;
s2 = 2*s2_inp - 1; 
sq = s1 + j*s2;   % QPSK modulated 

for x=1:length(snr_db)
    n = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));  % AWGN Noise
    sq_n = sq + 10^(-snr_db(x)/20)*n;    % Noise addition AWGN
    s1_out = real(sq_n) > 0; % Decoding
    s2_out = imag(sq_n) > 0; % Decoding
    diff1 = xor(s1_inp,s1_out);
    diff2 = xor(s2_inp,s2_out);
    err1 = sum(diff1(:)==1);
    err2 = sum(diff2(:)==1);
    ber(x) = mean([err1,err2]);    % Mean error
end

qpsk_ber = ber/N;   % simulated BER for QPSK modulation
qpsk_ser = 2*qpsk_ber;   % simulated SER for QPSK modulation

close all;
figure;
semilogy(snr_db, bpsk_ber, 'go-');
hold on
semilogy(snr_db, qpsk_ser, 'r*-');
hold off;
grid on;
axis([1 10 10^-6 1]);
legend('BPSK','QPSK');
title('Simulated plot under AWGN channel with BPSK and QPSK modulation.');
xlabel('SNR_{db}');
ylabel('Error Rate');
