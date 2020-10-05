close all;
clc;

% Number of bits to be transmitted
N = 10^6;

% SNR Range
Eb_N0_dB = (-5:20);
% Conversion of dBs to linear scale
Eb_N0_linear = power(10,Eb_N0_dB/10);


mu = sqrt(Eb_N0_linear./(Eb_N0_linear + 1));

% ro1 is for the imperfect CSI
ro1 = 0.9;
%ro2 is for the perfect CSI
ro2 = 1;

% analytical formula to find BER for SIMO having perfect CSI (becuse ro = 1)
analytical_pscsi_simo = 1/4*(2 - 3*ro2*mu + power(ro2*mu,3));

% analytical formula to find BER for SIMO having perfect CSI (becuse ro = 0.9)
analytical_icsi_simo = 1/4*(2 - 3*ro1*mu + power(ro1*mu,3));

% Generating input signal
input = rand(1,N)>0.5;
% BPSK modulation
signal = 2*input - 1;

for i=1:length(Eb_N0_dB)
    % AWGN Noise generation
    n1 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); 
    n2 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    
    % generating rayleigh random variable which is modelled for the channel coefficient
    hh1 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    hh2 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    
    % generating the channel coefficient for the given value of ro = 1 (perfect csi) and generated h for Rx no. 1
    h1 =  ro2*hh1 + sqrt(1-ro2^2)*1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    % generating the channel coefficient for the given value of ro = 1 (perfect csi) and generated h for Rx no. 2
    h2 =  ro2*hh2 + sqrt(1-ro2^2)*1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    
    % generating the channel coefficient for the given value of ro = 0.9 (imperfect csi) and generated h for Rx no. 1
    h3 =  ro1*hh1 + sqrt(1-ro1^2)*1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    % generating the channel coefficient for the given value of ro = 0.9 (imperfect csi) and generated h for Rx no. 2
    h4 =  ro1*hh2 + sqrt(1-ro1^2)*1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    

    p = power(10,(-Eb_N0_dB(i)/20));
    
    % recieved signal at reciever 1 along with the noise for ro = 1 (perfect csi)
    y1_pcsi = h1.*signal + p*n1;
    % recieved signal at reciever 2 along with the noise for ro = 1 (perfect csi)
    y2_pcsi = h2.*signal + p*n2;
    
    % recieved signal at reciever 1 along with the noise for ro = 0.9 (imperfect csi)
    y1_icsi = h3.*signal + p*n1;
    % recieved signal at reciever 2 along with the noise for ro = 0.9 (imperfect csi)
    y2_icsi = h4.*signal + p*n2;
    
    % using the detection rule for SIMO Wireless system for both perfect & imperfect csi
    y_pcsi = conj(h1).*y1_pcsi + conj(h2).*y2_pcsi;
    y_icsi = conj(hh1).*y1_icsi + conj(hh2).*y2_icsi; % here we are taking h1 & h2 as it is imperfect csi we don't know the values of h at reciever end
    
    % Applying MLD rule
    rsig_pcsi = real(y_pcsi) >= 0;
    % counting errors
    nerr_pcsi(i) = size(find([input-rsig_pcsi]), 2);
    
    % Applying MLD rule
    rsig_icsi = real(y_icsi) >= 0;
    % counting errors
    nerr_icsi(i) =  size(find([input-rsig_icsi]), 2);
    
end

% BER Calculation
ber_pcsi = nerr_pcsi/N;
ber_icsi = nerr_icsi/N;

% Plots
close all;
figure;
semilogy(Eb_N0_dB,analytical_pscsi_simo,'m.-');
hold on
semilogy(Eb_N0_dB,ber_pcsi,'kd');
hold on
semilogy(Eb_N0_dB,analytical_icsi_simo,'r.-');
hold on
semilogy(Eb_N0_dB,ber_icsi,'bo');
hold off

axis([-5 20 10^-5 1]);
LOC = "southwest";
legend('Analytical Perfect CSI Simo (L=2 ,\rho=1)','Simulation Perfect CSI Simo (L=2 ,\rho=1)','Analytical Perfect CSI Simo (L=2 ,\rho=0.9)','Simulation Perfect CSI Simo (L=2 ,\rho=0.9)','Location',LOC);
title('Plots','FontWeight','Normal');
xlabel('$\frac{E_b}{N_0} (dBs)$','Interpreter','latex');
ylabel('BER');