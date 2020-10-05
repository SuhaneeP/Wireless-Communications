close all;
clc;

% Number of bits to be transmitted
N = 10^6;

% Range of SNR
Eb_N0_dB = 0:35;
Eb_N0_linear = power(10,Eb_N0_dB/10);

% Analytical BER for Rayleigh Wireless Channel for SISO
nRx_1_analytical = 0.5*(1-sqrt(Eb_N0_linear./(Eb_N0_linear + 1)));

mu = sqrt(Eb_N0_linear./(Eb_N0_linear + 1));

% Analytical BER for Wireless Channel for SIMO ( 1 transmitter & 2 reciever)
nRx_2_analytical = 1/4*(2 - 3*mu + mu.^3);

% Generating signal
bits = rand(1,N)>0;
signal = 2*bits - 1; % BPSK modulation

for i=1:length(Eb_N0_dB)
    
    % AWGN noice generation
    n1 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    n2 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    n3 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    n4 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    
    % channel gain generation
    h1 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N)); 
    h2 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    h3 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    h4 = 1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
    
    % normalising gain on basis of SNR value
    p = power(10,-Eb_N0_dB(i)/20);
    
    % Signal at the recieving end
    y1 = h1.*signal + p*n1;
    y2 = h2.*signal + p*n2;
    y3 = h3.*signal + p*n3;
    y4 = h4.*signal + p*n4;
    
    % Recieved signal for SISO wireless
    rnRx_1 = y1./h1;
    
    % Recieved signal for SIMO wireless ( 1 x 2 )
    rnRx_2 = conj(h1).*y1 + conj(h2).*y2;
    
    % Recieved signal for SIMO wireless ( 1 x 4 )
    rnRx_4 = conj(h1).*y1 + conj(h2).*y2 + conj(h3).*y3 + conj(h4).*y4;
    
    % The following code decodes the signal by taking the real part and
    % counting the error bits.
    dbnRx_1 = real(rnRx_1) >0;
    dfnRx_1 = xor(dbnRx_1,bits);
    BER_nRx_1(i)=sum(dfnRx_1(:)==1); % Error count for SISO wireless
    
    dbnRx_2 = real(rnRx_2) >0;
    dfnRx_2 = xor(dbnRx_2,bits);
    BER_nRx_2(i)=sum(dfnRx_2(:)==1); % Error count for SIMO wireless (1x2)
    
    dbnRx_4 = real(rnRx_4) >0;
    dfnRx_4 = xor(dbnRx_4,bits);
    BER_nRx_4(i)=sum(dfnRx_4(:)==1); % Error count for SIMO wireless (1x4)
end

% BER calculation
BER_nRx_1_simulation = BER_nRx_1/N;
BER_nRx_2_simulation = BER_nRx_2/N;
BER_nRx_4_simulation = BER_nRx_4/N;

% Plots
close all;
figure;
semilogy(Eb_N0_dB,nRx_1_analytical,'b*-');
hold on
semilogy(Eb_N0_dB,BER_nRx_1_simulation,'md-');
hold on
semilogy(Eb_N0_dB,nRx_2_analytical,'rp-');
hold on
semilogy(Eb_N0_dB,BER_nRx_2_simulation,'bo-');
hold on
semilogy(Eb_N0_dB,BER_nRx_4_simulation,'g*-');
hold off

axis([0 35 10^-5 1]);
legend('nRx=1 Analytical','nRx=1 Simulation','nRx=2 Analytical','nRx=2 Simulation','nRx=4 Simulation');
title('SISO Wireless & SIMO Wireless (1 X 2 & 1 X 4)','FontWeight','Normal');
xlabel('$\frac{E_b}{N_0} (dBs)$','Interpreter','latex');
ylabel('BER');