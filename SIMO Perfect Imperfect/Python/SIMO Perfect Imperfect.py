import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


def generate_BER_vs_EbN0_Th(EbN0dB_range, rho):
    """Generate the theoretical BER array for a given signal and SNR range

    The values that fall in the given SNR range are converted to linear form and
    the formula of BER for (1*2) SIMO MRC is used to compute the BER value.

    Args:
        EbN0dB_range (range): range object containing the range of SNR
                              values used in simulations
        rho: correlation coefficient of CSI and ICSI

    Returns:
        CSI_Theory_BER (list): list of CSI theoretical BER values
        ICSI_Theory_BER (list): list of ICSI theoretical BER values
    """
    CSI_Theory_BER = []
    ICSI_Theory_BER = []
    for EbN0dB in EbN0dB_range:
        EbN0 = 10.0 ** (EbN0dB / 10.0)
        mu = np.sqrt(EbN0 / (EbN0 + 1))
        CSI_Theory_Pe = 0.25 * (2 - 3 * mu + mu ** 3)
        ICSI_Theory_Pe = 0.25 * (2 - 3 * rho * mu + (rho ** 3) * (mu ** 3))
        CSI_Theory_BER.append(CSI_Theory_Pe)
        ICSI_Theory_BER.append(ICSI_Theory_Pe)

    return CSI_Theory_BER, ICSI_Theory_BER


def generate_noise(signal, EbN0):
    """Generate the noise array for a given signal and its SNR

    The array of randomly generated values is multiplied with the noise variance
    which is computed using the SNR value. The array thus obtained is the noise array.

    Args:
        signal (array): signal array
        EbN0 (float): SNR value

    Returns:
        noise (array): The noise array
    """
    noise = (1 / np.sqrt(2 * EbN0)) * (
        np.random.normal(size=(np.shape(signal)))
        + np.random.normal(size=(np.shape(signal))) * 1j
    )
    return noise


def generate_BER_vs_EbN0_Pr(EbN0dB_range, rho):
    """Generate practical BER for a BPSK signal

    The randomly generated message is converted to its corresponding BPSK signal by
    changing all the 0s to -1. The resulting signal array denotes the phase shift
    property. The same signal is duplicated nRx times and the final signal array is
    generated. For wireless channel, we multiply the signal array with the h array
    which is an array of complex numbers belonging to the normal distribution.
    Noise of the same dimension is added to the signal. On the receiver end, we
    multiply the signal received on each antenna by their corresponding channel
    factor. The signals received are summed up and by checking the real parts of
    the sum array, the received message is obtained.

    Args:
        EbN0dB_range (range): range object containing the range of SNR
                              values used in simulations
        rho: correlation coefficient of CSI and ICSI

    Returns:
        CSI_BER (list): list of CSI BER values
        ICSI_BER (list): list of ICSI BER values
    """
    msg = np.random.randint(low=0, high=2, size=int(1e6))
    signal = np.ones(len(msg))
    for i, j in enumerate(msg):
        if j == 0:
            signal[i] = -1
    x = []
    x.append(signal)
    signal = np.asarray(x)

    CSI_BER = []
    ICSI_BER = []
    for EbN0dB in EbN0dB_range:
        EbN0 = 10 ** (EbN0dB / 10.0)
        signal_prime = np.repeat(signal, 2, axis=0)
        n = generate_noise(signal_prime, EbN0)
        CSI_h = (1 / np.sqrt(2)) * (
            np.random.normal(size=(np.shape(signal_prime)))
            + 1j * np.random.normal(size=(np.shape(signal_prime)))
        )
        ICSI_h = rho * CSI_h + (np.sqrt(1 - rho * rho)) * (1 / np.sqrt(2)) * (
            np.random.normal(size=(np.shape(signal_prime)))
            + 1j * np.random.normal(size=(np.shape(signal_prime)))
        )
        CSI_y = CSI_h * x + n
        ICSI_y = ICSI_h * x + n
        CSI_d = np.sum(np.conjugate(CSI_h) * CSI_y, axis=0)
        ICSI_d = np.sum(np.conjugate(CSI_h) * ICSI_y, axis=0)
        CSI_received_msg = np.real(CSI_d) >= 0
        ICSI_received_msg = np.real(ICSI_d) >= 0
        CSI_Pb_pr = np.count_nonzero(msg != CSI_received_msg) / len(msg)
        ICSI_Pb_pr = np.count_nonzero(msg != ICSI_received_msg) / len(msg)
        CSI_BER.append(CSI_Pb_pr)
        ICSI_BER.append(ICSI_Pb_pr)

    return CSI_BER, ICSI_BER


if __name__ == "__main__":
    EbN0dB_range = range(-5, 21)
    rho = 0.9
    CSI_Theory_BER, ICSI_Theory_BER = generate_BER_vs_EbN0_Th(EbN0dB_range, rho)
    CSI_BER, ICSI_BER = generate_BER_vs_EbN0_Pr(EbN0dB_range, rho)
    plt.plot(
        EbN0dB_range,
        CSI_Theory_BER,
        fillstyle="none",
        marker="o",
        linestyle="-",
        label=r"Analytical Perfect CSI SIMO (L = 2, $\rho$ = 1)",
    )
    plt.plot(
        EbN0dB_range,
        CSI_BER,
        marker="x",
        linestyle="-",
        label=r"Simulation Perfect CSI SIMO (L = 2, $\rho$ = 1)",
    )
    plt.plot(
        EbN0dB_range,
        ICSI_Theory_BER,
        fillstyle="none",
        marker="o",
        linestyle="-",
        label=r"Analytical ICSI SIMO (L = 2, $\rho$ = 0.9)",
    )
    plt.plot(
        EbN0dB_range,
        ICSI_BER,
        marker="x",
        linestyle="-",
        label=r"Simulation ICSI SIMO (L = 2, $\rho$ = 0.9)",
    )
    plt.xscale("linear")
    plt.xlabel("SNR (dB)")
    plt.ylabel("BER")
    plt.ylim(1e-5, 1e0)
    plt.yscale("log")
    plt.legend()
    plt.grid()
    plt.show()
