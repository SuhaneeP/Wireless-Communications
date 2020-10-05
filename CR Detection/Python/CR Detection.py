import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

if __name__ == "__main__":
    SNR_dB = range(-20, 5)
    L = 1000
    Pf = 0.01
    threshold_constant = (norm.isf(Pf) / np.sqrt(L)) + 1
    threshold = threshold_constant * np.ones(len(SNR_dB))

    Pd_Pr = [0 for i in SNR_dB]
    for i, snr_db in enumerate(SNR_dB):
        SNR = 10 ** (snr_db / 10)
        detect = 0
        for j in range(10000):
            noise = np.random.normal(size=L)
            signal = np.sqrt(SNR) * np.random.normal(size=L)
            received_signal = signal + noise
            received_signal_energy = np.abs(received_signal) ** 2
            test_statistic = (1 / L) * np.sum(received_signal_energy)
            if test_statistic > threshold[i]:
                detect += 1
        Pd_Pr[i] = detect / 10000

    Pd_Th = [
        norm.sf(
            ((threshold_constant - ((10 ** (snr / 10)) + 1)) * L)
            / (np.sqrt(2 * L * (((10 ** (snr / 10)) + 1) ** 2)))
        )
        for snr in SNR_dB
    ]
    plt.plot(SNR_dB, Pd_Th, "-", label="Theoretical")
    plt.plot(SNR_dB, Pd_Pr, "o", label="Simulations")
    plt.xlabel("SNR (dB)")
    plt.ylabel("Detection Probability")
    plt.legend()
    plt.grid()
    plt.savefig("Figure_1.png")
    plt.show()
